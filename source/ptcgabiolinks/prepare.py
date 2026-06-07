"""Prepare downloaded GDC files into Python data objects."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .barcode import col_data_prepare
from .client import get_gdc_info, is_serve_ok
from .models import GDCQuery, SummarizedExperiment

STAR_ASSAY_COLUMNS = {
    "unstranded": "unstranded",
    "stranded_first": "stranded_first",
    "stranded_second": "stranded_second",
    "tpm_unstranded": "tpm_unstrand",
    "fpkm_unstranded": "fpkm_unstrand",
    "fpkm_uq_unstranded": "fpkm_uq_unstrand",
}


def GDCprepare(
    query: GDCQuery,
    save: bool = False,
    save_filename: str | None = None,
    directory: str | Path = "GDCdata",
    summarized_experiment: bool = True,
    remove_files_prepared: bool = False,
    **kwargs: Any,
) -> SummarizedExperiment | pd.DataFrame:
    """Read downloaded files and prepare an object like TCGAbiolinks."""

    save_filename = _pop_alias(kwargs, save_filename, "save.filename")
    summarized_experiment = _pop_alias(
        kwargs, summarized_experiment, "summarizedExperiment"
    )
    remove_files_prepared = _pop_alias(
        kwargs, remove_files_prepared, "remove.files.prepared"
    )
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise TypeError(f"Unknown GDCprepare argument(s): {unknown}")
    if query is None:
        raise ValueError("Please set query parameter")
    if remove_files_prepared and not save:
        raise ValueError("To remove files, please set save=True")

    is_serve_ok()
    files = _query_files(query, Path(directory))
    missing = [str(path) for path in files if not path.exists()]
    if missing:
        preview = "\n".join(missing[:5])
        raise FileNotFoundError(
            "I couldn't find all the files from the query. "
            "Please check if directory is right or GDCdownload downloaded the samples.\n"
            f"Missing examples:\n{preview}"
        )

    if _is_star_counts_query(query):
        prepared = _read_star_counts(
            files,
            _cases_for_prepare(query),
            summarized_experiment,
            projects=query.results["project"].astype(str).tolist(),
            sample_types=query.results.get("sample_type", pd.Series(dtype=str)).astype(str).tolist(),
        )
    else:
        raise NotImplementedError(
            "This Python port currently prepares Transcriptome Profiling / "
            "Gene Expression Quantification / STAR - Counts."
        )

    if isinstance(prepared, SummarizedExperiment):
        prepared.metadata["data_release"] = get_gdc_info().get("data_release")

    if save:
        import pickle

        if save_filename is None:
            project = "_".join(query.project)
            save_filename = f"{project}_{query.data_category.replace(' ', '_')}.pkl"
        with open(save_filename, "wb") as handle:
            pickle.dump(prepared, handle)

    if remove_files_prepared:
        for path in files:
            path.unlink(missing_ok=True)

    return prepared


def assay(data: SummarizedExperiment | pd.DataFrame, name: str | None = None) -> pd.DataFrame:
    if isinstance(data, SummarizedExperiment):
        return data.assay(name)
    if name is not None:
        raise TypeError("name can only be used with SummarizedExperiment objects")
    return data


def assays(data: SummarizedExperiment) -> dict[str, pd.DataFrame]:
    return data.assays


def _is_star_counts_query(query: GDCQuery) -> bool:
    return (
        query.data_category.lower() == "transcriptome profiling"
        and (query.data_type or "").lower() == "gene expression quantification"
        and (query.workflow_type or "").lower() == "star - counts"
    )


def _query_files(query: GDCQuery, directory: Path) -> list[Path]:
    files: list[Path] = []
    for _, row in query.results.iterrows():
        files.append(
            directory
            / str(row["project"])
            / _space_to_underscore(str(row["data_category"]))
            / _space_to_underscore(str(row["data_type"]))
            / _space_to_underscore(str(row["file_id"]))
            / _space_to_underscore(str(row["file_name"]))
        )
    return files


def _cases_for_prepare(query: GDCQuery) -> list[str]:
    rows = query.results
    if all(str(project).startswith(("TCGA", "TARGET", "CGCI-HTMCP-CC", "CPTAC-2", "BEATAML1.0-COHORT")) for project in rows["project"]):
        return rows["cases"].astype(str).tolist()
    if "sample.submitter_id" in rows:
        return rows["sample.submitter_id"].astype(str).tolist()
    return rows["cases"].astype(str).tolist()


def _read_star_counts(
    files: list[Path],
    cases: list[str],
    summarized_experiment: bool,
    projects: list[str] | None = None,
    sample_types: list[str] | None = None,
) -> SummarizedExperiment | pd.DataFrame:
    if len(files) != len(cases):
        raise ValueError("files and cases must have the same length")

    per_case: list[pd.DataFrame] = []
    for path, case in zip(files, cases, strict=True):
        data = _read_star_counts_file(path)
        data["case_barcode"] = case
        per_case.append(data)

    first = per_case[0]
    genes = first.loc[first["gene_id"].astype(str).str.startswith("ENSG"), ["gene_id", "gene_name", "gene_type"]].copy()
    genes = genes.drop_duplicates("gene_id")
    genes = genes.set_index("gene_id", drop=False)

    metric_frames = [
        frame.loc[~frame["gene_id"].astype(str).str.startswith("ENSG")].assign(case_barcode=case)
        for frame, case in zip(per_case, cases, strict=True)
    ]
    metrics = pd.concat(metric_frames, ignore_index=True)

    assays_map: dict[str, pd.DataFrame] = {}
    for source_col, assay_name in STAR_ASSAY_COLUMNS.items():
        columns: dict[str, pd.Series] = {}
        for frame, case in zip(per_case, cases, strict=True):
            gene_rows = frame.loc[frame["gene_id"].astype(str).str.startswith("ENSG")]
            series = gene_rows.set_index("gene_id")[source_col]
            columns[case] = pd.to_numeric(series, errors="coerce")
        assay_frame = pd.DataFrame(columns).reindex(genes.index)
        assays_map[assay_name] = assay_frame

    if not summarized_experiment:
        merged = genes.copy()
        for assay_name, assay_frame in assays_map.items():
            renamed = assay_frame.add_prefix(f"{assay_name}_")
            merged = merged.join(renamed)
        return merged

    col_data = col_data_prepare(cases)
    if projects is not None and len(projects) == len(cases):
        col_data["project_id"] = projects
    if sample_types is not None and len(sample_types) == len(cases):
        col_data["sample_type"] = sample_types

    return SummarizedExperiment(
        assays=assays_map,
        row_data=genes,
        col_data=col_data,
        metadata={"metrics": metrics},
    )


def _read_star_counts_file(path: Path) -> pd.DataFrame:
    skiprows = _find_star_header(path)
    data = pd.read_csv(path, sep="\t", skiprows=skiprows, comment="#")
    missing = set(STAR_ASSAY_COLUMNS) - set(data.columns)
    required = {"gene_id", "gene_name", "gene_type"}
    missing_required = required - set(data.columns)
    if missing_required:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_required)}")
    if missing:
        raise ValueError(f"{path} is missing STAR assay columns: {sorted(missing)}")
    return data


def _find_star_header(path: Path) -> int:
    with path.open("rt", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            if line.startswith("gene_id\t"):
                return idx
    return 0


def _space_to_underscore(value: str) -> str:
    return value.replace(" ", "_")


def _pop_alias(kwargs: dict[str, Any], value: Any, *aliases: str) -> Any:
    for alias in aliases:
        if alias in kwargs:
            if value is not None:
                raise TypeError(f"Received both {alias!r} and Python alias")
            value = kwargs.pop(alias)
    return value
