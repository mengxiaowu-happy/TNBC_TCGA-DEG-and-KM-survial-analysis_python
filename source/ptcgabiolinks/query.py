"""GDC query implementation."""

from __future__ import annotations

import json
import re
import warnings
from typing import Any

import pandas as pd

from .client import get_json, is_serve_ok
from .models import GDCQuery

PAGE_SIZE = 1000


def GDCquery(
    project: str | list[str],
    data_category: str | None = None,
    data_type: str | None = None,
    workflow_type: str | None = None,
    access: str | list[str] | None = None,
    platform: str | list[str] | None = None,
    barcode: str | list[str] | None = None,
    data_format: str | None = None,
    experimental_strategy: str | list[str] | None = None,
    sample_type: str | list[str] | None = None,
    **kwargs: Any,
) -> GDCQuery:
    """Search GDC files using arguments modeled after R TCGAbiolinks.

    Python keywords use underscores, but R-style names such as
    ``"data.category"`` are accepted through ``**kwargs``.
    """

    data_category = _pop_alias(kwargs, data_category, "data.category")
    data_type = _pop_alias(kwargs, data_type, "data.type")
    workflow_type = _pop_alias(kwargs, workflow_type, "workflow.type")
    data_format = _pop_alias(kwargs, data_format, "data.format")
    experimental_strategy = _pop_alias(
        kwargs, experimental_strategy, "experimental.strategy"
    )
    sample_type = _pop_alias(kwargs, sample_type, "sample.type")
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise TypeError(f"Unknown GDCquery argument(s): {unknown}")
    if data_category is None:
        raise ValueError("Please set data_category (R: data.category)")

    is_serve_ok()
    projects = _as_list(project)
    access_values = _as_optional_list(access)
    platform_values = _as_optional_list(platform)
    barcode_values = _as_optional_list(barcode)
    experimental_values = _as_optional_list(experimental_strategy)
    sample_values = _as_optional_list(sample_type)

    records: list[dict[str, Any]] = []
    for proj in projects:
        hits = _fetch_files_for_project(
            project=proj,
            data_category=data_category,
            data_type=data_type,
            workflow_type=workflow_type,
            access=access_values,
            platform=platform_values,
            experimental_strategy=experimental_values,
            sample_type=sample_values,
        )
        records.extend(_flatten_hit(hit, proj, data_category) for hit in hits)

    results = pd.DataFrame.from_records(records)
    if results.empty:
        raise RuntimeError("Sorry, no results were found for this query")

    results = _post_filter(results, "data_category", data_category)
    if data_type:
        results = _post_filter(results, "data_type", data_type)
    if workflow_type:
        results = _post_filter(results, "analysis_workflow_type", workflow_type)
    if data_format:
        results = _post_filter(results, "data_format", data_format)
    if access_values:
        results = _post_filter(results, "access", access_values)
    if experimental_values:
        results = _post_filter(results, "experimental_strategy", experimental_values)
    if platform_values and "platform" in results:
        results = _post_filter(results, "platform", platform_values)
    if barcode_values:
        results = _filter_by_barcode(results, barcode_values)
    if sample_values:
        results = _filter_by_sample_type(results, sample_values)

    if results.empty:
        raise RuntimeError("Sorry, no results were found for this query")
    if results["cases"].duplicated().any():
        warnings.warn(
            "There are more than one file for the same case. Please verify query results.",
            RuntimeWarning,
            stacklevel=2,
        )

    return GDCQuery(
        results=results.reset_index(drop=True),
        project=projects,
        data_category=data_category,
        data_type=data_type,
        workflow_type=workflow_type,
        access=access_values,
        experimental_strategy=experimental_values,
        platform=platform_values,
        sample_type=sample_values,
        barcode=barcode_values,
        data_format=data_format,
    )


def getResults(
    query: GDCQuery,
    rows: list[int] | slice | None = None,
    cols: list[str] | str | None = None,
) -> pd.DataFrame | pd.Series:
    data = query.results
    if rows is not None:
        data = data.iloc[rows]
    if cols is not None:
        data = data.loc[:, cols]
    return data


def _fetch_files_for_project(
    project: str,
    data_category: str,
    data_type: str | None,
    workflow_type: str | None,
    access: list[str] | None,
    platform: list[str] | None,
    experimental_strategy: list[str] | None,
    sample_type: list[str] | None,
) -> list[dict[str, Any]]:
    filters = _build_filters(
        project=project,
        data_category=data_category,
        data_type=data_type,
        workflow_type=workflow_type,
        access=access,
        platform=platform,
        experimental_strategy=experimental_strategy,
        sample_type=sample_type,
    )
    expand = (
        "cases,cases.project,center,analysis"
        if data_category in {"Clinical", "Biospecimen"}
        else "cases,cases.samples.portions.analytes.aliquots,cases.project,center,analysis,cases.samples"
    )
    hits: list[dict[str, Any]] = []
    offset = 0
    while True:
        payload = get_json(
            "/files",
            params={
                "format": "JSON",
                "pretty": "true",
                "expand": expand,
                "filters": json.dumps(filters, separators=(",", ":")),
                "size": PAGE_SIZE,
                "from": offset,
            },
        )
        data = payload.get("data", {})
        page_hits = data.get("hits", [])
        hits.extend(page_hits)
        pagination = data.get("pagination", {})
        total = int(pagination.get("total", pagination.get("count", len(hits))))
        if len(hits) >= total or not page_hits:
            break
        offset += len(page_hits)
    return hits


def _build_filters(
    project: str,
    data_category: str,
    data_type: str | None,
    workflow_type: str | None,
    access: list[str] | None,
    platform: list[str] | None,
    experimental_strategy: list[str] | None,
    sample_type: list[str] | None,
) -> dict[str, Any]:
    content: list[dict[str, Any]] = [
        _in_filter("cases.project.project_id", [project]),
        _in_filter("files.data_category", [data_category]),
    ]
    if data_type:
        content.append(_in_filter("files.data_type", [data_type]))
    if workflow_type:
        content.append(_in_filter("files.analysis.workflow_type", [workflow_type]))
    if access:
        content.append(_in_filter("files.access", access))
    if platform:
        content.append(_in_filter("files.platform", platform))
    if experimental_strategy:
        content.append(
            _in_filter("files.experimental_strategy", experimental_strategy)
        )
    if sample_type:
        normalized = [
            "Primary Tumor" if value == "Primary solid Tumor" else value
            for value in sample_type
        ]
        content.append(_in_filter("cases.samples.sample_type", normalized))
    return {"op": "and", "content": content}


def _in_filter(field: str, values: list[str]) -> dict[str, Any]:
    return {"op": "in", "content": {"field": field, "value": values}}


def _flatten_hit(hit: dict[str, Any], project: str, data_category: str) -> dict[str, Any]:
    row: dict[str, Any] = {}
    for key, value in hit.items():
        if key in {"cases", "analysis", "center", "acl"}:
            continue
        if isinstance(value, (dict, list)):
            continue
        row[key] = value

    row["file_id"] = hit.get("file_id") or hit.get("id")
    row["file_name"] = hit.get("file_name") or hit.get("filename")
    row["project"] = project

    analysis = hit.get("analysis") or {}
    if isinstance(analysis, dict):
        for key, value in analysis.items():
            if key.endswith("datetime"):
                continue
            out_key = "analysis_id" if key == "analysis_id" else f"analysis_{key}"
            row[out_key] = value

    center = hit.get("center") or {}
    if isinstance(center, dict):
        for key, value in center.items():
            out_key = key if key == "center_id" else f"center_{key}"
            row[out_key] = value

    cases = hit.get("cases") or []
    case_submitters: list[str] = []
    sample_submitters: list[str] = []
    aliquot_submitters: list[str] = []
    sample_types: list[str] = []
    ffpe_values: list[bool] = []

    for case in cases:
        if not isinstance(case, dict):
            continue
        if case.get("submitter_id"):
            case_submitters.append(case["submitter_id"])
        for sample in case.get("samples") or []:
            if not isinstance(sample, dict):
                continue
            if sample.get("submitter_id"):
                sample_submitters.append(sample["submitter_id"])
            if sample.get("sample_type"):
                sample_types.append(sample["sample_type"])
            if "is_ffpe" in sample and sample["is_ffpe"] is not None:
                ffpe_values.append(bool(sample["is_ffpe"]))
            for portion in sample.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for aliquot in analyte.get("aliquots") or []:
                        if isinstance(aliquot, dict) and aliquot.get("submitter_id"):
                            aliquot_submitters.append(aliquot["submitter_id"])

    row["cases.submitter_id"] = ";".join(_unique(case_submitters))
    row["sample.submitter_id"] = ";".join(_unique(sample_submitters))
    row["sample_type"] = ";".join(_unique(sample_types))
    row["is_ffpe"] = any(ffpe_values) if ffpe_values else None

    if data_category == "Clinical":
        row["cases"] = ";".join(_unique(case_submitters))
    elif data_category == "Biospecimen":
        row["cases"] = ",".join(_unique(case_submitters))
    else:
        identifiers = aliquot_submitters or sample_submitters or case_submitters
        row["cases"] = ";".join(_unique(identifiers))
    return row


def _post_filter(
    results: pd.DataFrame, column: str, value: str | list[str]
) -> pd.DataFrame:
    if column not in results.columns:
        raise ValueError(f"GDC results do not include column {column!r}")
    values = _as_list(value)
    mask = results[column].astype(str).str.lower().isin([v.lower() for v in values])
    if not mask.any():
        available = "\n  => ".join(sorted(results[column].dropna().astype(str).unique()))
        raise ValueError(f"Please set a valid {column} argument:\n  => {available}")
    return results.loc[mask].copy()


def _filter_by_barcode(results: pd.DataFrame, barcodes: list[str]) -> pd.DataFrame:
    mask = pd.Series(False, index=results.index)
    cases = results["cases"].fillna("").astype(str)
    for barcode in barcodes:
        mask |= cases.str.contains(re.escape(barcode), case=False, regex=True)
    if not mask.any():
        available = "\n  => ".join(cases.unique())
        raise ValueError(f"None of the barcodes were matched. Available barcodes:\n  => {available}")
    return results.loc[mask].copy()


def _filter_by_sample_type(results: pd.DataFrame, sample_types: list[str]) -> pd.DataFrame:
    wanted = {value.lower() for value in sample_types}

    def matches(value: object) -> bool:
        parts = [part.strip().lower() for part in str(value).split(";")]
        return any(part in wanted for part in parts)

    mask = results["sample_type"].map(matches)
    if not mask.any():
        available = "\n  => ".join(
            sorted(results["sample_type"].dropna().astype(str).unique())
        )
        raise ValueError(f"Please set a valid sample_type argument:\n  => {available}")
    return results.loc[mask].copy()


def _pop_alias(kwargs: dict[str, Any], value: Any, *aliases: str) -> Any:
    for alias in aliases:
        if alias in kwargs:
            if value is not None:
                raise TypeError(f"Received both {alias!r} and Python alias")
            value = kwargs.pop(alias)
    return None if value is False else value


def _as_list(value: str | list[str]) -> list[str]:
    if isinstance(value, str):
        return [value]
    return list(value)


def _as_optional_list(value: str | list[str] | None) -> list[str] | None:
    if value is None or value is False:
        return None
    values = _as_list(value)
    return None if not values else values


def _unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    unique_values: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            unique_values.append(value)
    return unique_values
