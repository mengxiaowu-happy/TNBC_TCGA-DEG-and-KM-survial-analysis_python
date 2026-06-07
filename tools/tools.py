"""Reusable tools for the TNBC expression-survival workflow.

The functions here intentionally mirror the data handling in the supplied R
script while making the cleaning rules explicit and testable.
"""

from __future__ import annotations

import json
import math
from io import StringIO
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

from ptcgabiolinks.client import get_json
from ptcgabiolinks.models import GDCQuery
from ptcgabiolinks.query import _flatten_hit


TNBC_ER_COL = "ER_STATUS_BY_IHC"
TNBC_PR_COL = "PR_STATUS_BY_IHC"
TNBC_HER2_COL = "IHC_HER2"
PATIENT_COL = "patientId"
CBIOPORTAL_API_ROOT = "https://www.cbioportal.org/api"
DEFAULT_CBIOPORTAL_STUDY = "brca_tcga"


def filter_tnbc_samples(
    clinical: pd.DataFrame,
    er_col: str = TNBC_ER_COL,
    pr_col: str = TNBC_PR_COL,
    her2_col: str = TNBC_HER2_COL,
    patient_col: str = PATIENT_COL,
    negative_value: str = "Negative",
) -> pd.DataFrame:
    """Keep triple-negative BRCA samples and drop rows without patient IDs."""

    _require_columns(clinical, [er_col, pr_col, her2_col, patient_col])
    mask = (
        clinical[er_col].eq(negative_value)
        & clinical[pr_col].eq(negative_value)
        & clinical[her2_col].eq(negative_value)
        & clinical[patient_col].notna()
    )
    return clinical.loc[mask].copy()


def discover_star_counts_files(
    input_path: str | Path,
    project: str = "TCGA-BRCA",
) -> pd.DataFrame:
    """Discover STAR Counts files in a GDCdownload-style directory tree."""

    root = Path(input_path)
    if not root.exists():
        raise FileNotFoundError(f"Input path does not exist: {root}")

    files = sorted(root.rglob("*.rna_seq.augmented_star_gene_counts.tsv"))
    if not files:
        files = sorted(root.rglob("*star_gene_counts.tsv"))
    records: list[dict[str, object]] = []
    for file_path in files:
        file_id = file_path.parent.name
        records.append(
            {
                "project": project,
                "data_category": "Transcriptome Profiling",
                "data_type": "Gene Expression Quantification",
                "analysis_workflow_type": "STAR - Counts",
                "file_id": file_id,
                "file_name": file_path.name,
                "file_path": file_path,
            }
        )
    if not records:
        raise FileNotFoundError(
            "No STAR Counts files were found. Expected files like "
            "<file_id>/*.rna_seq.augmented_star_gene_counts.tsv"
        )
    return pd.DataFrame(records)


def infer_gdc_download_root(input_path: str | Path, project: str = "TCGA-BRCA") -> Path:
    """Infer the directory argument expected by GDCprepare from a nested path."""

    path = Path(input_path)
    parts = list(path.parts)
    if project in parts:
        index = parts.index(project)
        if index == 0:
            return Path(".")
        return Path(*parts[:index])
    if (path / project).exists():
        return path
    return path


def build_query_from_gdc_directory(
    input_path: str | Path,
    project: str = "TCGA-BRCA",
    metadata_chunk_size: int = 100,
) -> GDCQuery:
    """Create a GDCQuery from already-downloaded GDC STAR Counts files."""

    discovered = discover_star_counts_files(input_path, project=project)
    metadata = fetch_gdc_file_metadata(
        discovered["file_id"].astype(str).tolist(),
        project=project,
        chunk_size=metadata_chunk_size,
    )
    merged = discovered.drop(columns=["file_path"]).merge(
        metadata,
        on="file_id",
        how="left",
        suffixes=("", "_metadata"),
    )
    for column in [
        "file_name",
        "project",
        "data_category",
        "data_type",
        "analysis_workflow_type",
    ]:
        metadata_column = f"{column}_metadata"
        if metadata_column in merged.columns:
            merged[column] = merged[metadata_column].combine_first(merged[column])
            merged = merged.drop(columns=metadata_column)

    if "cases" not in merged.columns or merged["cases"].isna().any():
        missing = merged.loc[merged.get("cases", pd.Series(index=merged.index)).isna(), "file_id"].tolist()
        preview = ", ".join(str(value) for value in missing[:5])
        raise RuntimeError(
            "Could not map all downloaded files to sample barcodes through GDC API. "
            f"Missing file_id examples: {preview}"
        )

    return GDCQuery(
        results=merged.reset_index(drop=True),
        project=[project],
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
    )


def fetch_gdc_file_metadata(
    file_ids: Sequence[str],
    project: str = "TCGA-BRCA",
    chunk_size: int = 100,
) -> pd.DataFrame:
    """Fetch GDC metadata for downloaded file IDs."""

    unique_ids = list(dict.fromkeys(str(file_id) for file_id in file_ids))
    records: list[dict[str, object]] = []
    for start in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[start : start + chunk_size]
        filters = {
            "op": "in",
            "content": {"field": "files.file_id", "value": chunk},
        }
        payload = get_json(
            "/files",
            params={
                "format": "JSON",
                "expand": "cases,cases.samples.portions.analytes.aliquots,cases.project,center,analysis,cases.samples",
                "filters": json.dumps(filters, separators=(",", ":")),
                "size": len(chunk),
            },
        )
        hits = payload.get("data", {}).get("hits", [])
        for hit in hits:
            cases = hit.get("cases") or [{}]
            hit_project = cases[0].get("project", {}).get("project_id", project)
            records.append(_flatten_hit(hit, hit_project, "Transcriptome Profiling"))
    if not records:
        raise RuntimeError("GDC API returned no metadata for discovered files")
    return pd.DataFrame.from_records(records)


def fetch_cbioportal_clinical(
    study_id: str = DEFAULT_CBIOPORTAL_STUDY,
    clinical_data_type: str = "PATIENT",
    api_root: str = CBIOPORTAL_API_ROOT,
    timeout: int = 120,
) -> pd.DataFrame:
    """Fetch cBioPortal clinical data and pivot it to one row per patient."""

    url = f"{api_root.rstrip('/')}/studies/{study_id}/clinical-data"
    response = requests.get(
        url,
        params={
            "clinicalDataType": clinical_data_type,
            "projection": "SUMMARY",
            "pageSize": 100000,
        },
        timeout=timeout,
    )
    response.raise_for_status()
    long_frame = pd.DataFrame(response.json())
    return pivot_cbioportal_clinical(long_frame)


def fetch_gdc_clinical_cases(
    project: str = "TCGA-BRCA",
    page_size: int = 200,
) -> list[dict[str, object]]:
    """Download all GDC case-level clinical records for a project."""

    filters = {"op": "in", "content": {"field": "project.project_id", "value": [project]}}
    expand = ",".join(
        [
            "demographic",
            "diagnoses",
            "diagnoses.treatments",
            "exposures",
            "follow_ups",
            "follow_ups.molecular_tests",
            "samples",
        ]
    )
    cases: list[dict[str, object]] = []
    offset = 0
    while True:
        payload = get_json(
            "/cases",
            params={
                "format": "JSON",
                "filters": json.dumps(filters, separators=(",", ":")),
                "expand": expand,
                "size": page_size,
                "from": offset,
            },
        )
        data = payload.get("data", {})
        hits = data.get("hits", [])
        cases.extend(hits)
        total = int(data.get("pagination", {}).get("total", len(cases)))
        if len(cases) >= total or not hits:
            break
        offset += len(hits)
    return cases


def fetch_gdc_clinical(
    project: str = "TCGA-BRCA",
    raw_output_file: str | Path | None = None,
) -> pd.DataFrame:
    """Fetch and normalize GDC clinical data for BRCA-style survival analysis."""

    cases = fetch_gdc_clinical_cases(project=project)
    if raw_output_file is not None:
        Path(raw_output_file).write_text(json.dumps(cases, ensure_ascii=False, indent=2), encoding="utf-8")
    return parse_gdc_brca_clinical(cases)


def parse_gdc_brca_clinical(cases: Sequence[dict[str, object]]) -> pd.DataFrame:
    """Parse GDC BRCA case records into the columns used by the R workflow."""

    rows: list[dict[str, object]] = []
    for case in cases:
        patient_id = str(case.get("submitter_id") or "")
        diagnoses = _as_list_of_dicts(case.get("diagnoses"))
        follow_ups = _as_list_of_dicts(case.get("follow_ups"))
        demographic = case.get("demographic") if isinstance(case.get("demographic"), dict) else {}

        receptor_status = _extract_receptor_status(follow_ups)
        os_status, os_months = _derive_os(demographic, diagnoses, follow_ups)
        dfs_status, dfs_months = _derive_dfs(diagnoses, follow_ups)

        rows.append(
            {
                "patientId": patient_id,
                "case_id": case.get("case_id"),
                "ER_STATUS_BY_IHC": receptor_status.get("ESR1"),
                "PR_STATUS_BY_IHC": receptor_status.get("PGR"),
                "IHC_HER2": receptor_status.get("ERBB2"),
                "OS_STATUS": os_status,
                "OS_MONTHS": os_months,
                "DFS_STATUS": dfs_status,
                "DFS_MONTHS": dfs_months,
                "vital_status": demographic.get("vital_status"),
                "days_to_death": demographic.get("days_to_death"),
            }
        )
    return pd.DataFrame(rows)


def pivot_cbioportal_clinical(long_frame: pd.DataFrame) -> pd.DataFrame:
    """Pivot cBioPortal clinical data from attribute rows to patient rows."""

    _require_columns(long_frame, ["patientId", "clinicalAttributeId", "value"])
    wide = (
        long_frame.pivot_table(
            index="patientId",
            columns="clinicalAttributeId",
            values="value",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    return wide


def gene_info_from_star_row_data(row_data: pd.DataFrame) -> pd.DataFrame:
    """Build an Ensembl-to-symbol table from downloaded STAR Counts row data."""

    _require_columns(row_data, ["gene_id", "gene_name"])
    gene_info = row_data.loc[:, ["gene_id", "gene_name"]].copy()
    gene_info["ensembl_gene_id"] = gene_info["gene_id"].astype(str).str.replace(
        r"\.[0-9]+$", "", regex=True
    )
    gene_info["hgnc_symbol"] = gene_info["gene_name"].astype(str)
    gene_info = gene_info.loc[
        gene_info["ensembl_gene_id"].notna()
        & gene_info["hgnc_symbol"].notna()
        & gene_info["hgnc_symbol"].ne("")
    ]
    return gene_info.loc[:, ["ensembl_gene_id", "hgnc_symbol"]].drop_duplicates()


def simplify_tcga_sample_ids(expression: pd.DataFrame, length: int = 12) -> pd.DataFrame:
    """Rename expression matrix columns to patient-level TCGA IDs."""

    simplified = expression.copy()
    simplified.columns = [str(column)[:length] for column in simplified.columns]
    return simplified


def match_expression_to_clinical(
    expression: pd.DataFrame,
    clinical: pd.DataFrame,
    patient_col: str = PATIENT_COL,
) -> pd.DataFrame:
    """Keep expression columns whose patient ID exists in the clinical table."""

    _require_columns(clinical, [patient_col])
    patients = set(clinical[patient_col].dropna().astype(str))
    return expression.loc[:, [column in patients for column in expression.columns]].copy()


def drop_duplicate_sample_columns(expression: pd.DataFrame, keep: str = "first") -> pd.DataFrame:
    """Remove duplicated sample columns after aliquot IDs are collapsed."""

    if keep != "first":
        raise ValueError("Only keep='first' is implemented to reproduce the R script")
    return expression.loc[:, ~pd.Index(expression.columns).duplicated(keep=keep)].copy()


def strip_ensembl_versions(expression: pd.DataFrame) -> pd.DataFrame:
    """Remove Ensembl version suffixes such as '.17' from row IDs."""

    stripped = expression.copy()
    stripped.index = stripped.index.astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    return stripped


def map_ensembl_to_hgnc(
    gene_ids: Sequence[str],
    host: str = "https://www.ensembl.org/biomart/martservice",
    chunk_size: int = 400,
    timeout: int = 120,
) -> pd.DataFrame:
    """Map Ensembl gene IDs to HGNC symbols using Ensembl BioMart."""

    unique_ids = list(dict.fromkeys(str(gene_id) for gene_id in gene_ids if pd.notna(gene_id)))
    records: list[pd.DataFrame] = []
    for start in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[start : start + chunk_size]
        xml = _biomart_query_xml(chunk)
        response = requests.post(host, data={"query": xml}, timeout=timeout)
        response.raise_for_status()
        if not response.text.strip():
            continue
        frame = pd.read_csv(
            StringIO(response.text),
            sep="\t",
            header=None,
            names=["ensembl_gene_id", "hgnc_symbol"],
            dtype=str,
        )
        records.append(frame)
    if not records:
        return pd.DataFrame(columns=["ensembl_gene_id", "hgnc_symbol"])
    return pd.concat(records, ignore_index=True).drop_duplicates()


def convert_ensembl_to_hgnc(
    expression: pd.DataFrame,
    gene_info: pd.DataFrame,
    duplicate_rule: str = "first",
) -> pd.DataFrame:
    """Attach HGNC symbols, remove unmapped genes, and use symbols as row names."""

    _require_columns(gene_info, ["ensembl_gene_id", "hgnc_symbol"])
    exp = expression.copy()
    exp["ensembl_gene_id"] = exp.index.astype(str)
    merged = exp.reset_index(drop=True).merge(gene_info, on="ensembl_gene_id", how="left")
    merged = merged.loc[merged["hgnc_symbol"].notna() & merged["hgnc_symbol"].ne("")].copy()

    sample_cols = [column for column in merged.columns if column not in {"ensembl_gene_id", "hgnc_symbol"}]
    if duplicate_rule == "first":
        merged = merged.loc[~merged["hgnc_symbol"].duplicated(keep="first")]
    elif duplicate_rule == "max_mean":
        numeric = merged[sample_cols].apply(pd.to_numeric, errors="coerce")
        merged["_mean_expression"] = numeric.mean(axis=1)
        merged = (
            merged.sort_values("_mean_expression", ascending=False)
            .loc[lambda frame: ~frame["hgnc_symbol"].duplicated(keep="first")]
            .drop(columns="_mean_expression")
        )
    else:
        raise ValueError("duplicate_rule must be 'first' or 'max_mean'")

    merged = merged.set_index("hgnc_symbol")
    return merged.loc[:, sample_cols].apply(pd.to_numeric, errors="coerce")


def log2_transform(expression: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """Apply log2(x + pseudocount) to expression values."""

    numeric = expression.apply(pd.to_numeric, errors="coerce")
    return np.log2(numeric + pseudocount)


def subset_expression_genes(expression: pd.DataFrame, gene_list: Sequence[str]) -> pd.DataFrame:
    """Keep genes in the supplied gene list, preserving expression row order."""

    wanted = set(gene_list)
    return expression.loc[[gene for gene in expression.index if gene in wanted]].copy()


def build_survival_expression_data(
    clinical: pd.DataFrame,
    log_expression: pd.DataFrame,
    gene_list: Sequence[str],
    patient_col: str = PATIENT_COL,
    status_col: str = "DFS_STATUS",
    time_col: str = "DFS_MONTHS",
) -> pd.DataFrame:
    """Merge selected gene expression with clinical survival fields."""

    _require_columns(clinical, [patient_col, status_col, time_col])
    expression = subset_expression_genes(log_expression, gene_list).T
    expression[patient_col] = expression.index.astype(str)
    clinical_subset = clinical.loc[:, [patient_col, status_col, time_col]].copy()
    merged = clinical_subset.merge(expression, on=patient_col, how="inner")
    merged = merged.loc[merged[status_col].notna()].copy()
    return merged


def analyze_gene_survival(
    merged_data: pd.DataFrame,
    gene_list: Sequence[str],
    output_dir: str | Path | None = None,
    patient_col: str = PATIENT_COL,
    status_col: str = "DFS_STATUS",
    time_col: str = "DFS_MONTHS",
    event_value: str = "1:Recurred/Progressed",
    high_label: str = "High",
    low_label: str = "Low",
    plot: bool = True,
    risk_table: bool = True,
) -> pd.DataFrame:
    """Median-split each gene, run log-rank tests, and optionally save KM plots."""

    _require_columns(merged_data, [patient_col, status_col, time_col])
    output_path = Path(output_dir) if output_dir is not None else None
    if output_path is not None:
        output_path.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, object]] = []
    for gene in gene_list:
        if gene not in merged_data.columns:
            results.append({"gene": gene, "status": "missing_gene", "p_value": np.nan})
            continue

        frame = merged_data.loc[:, [patient_col, status_col, time_col, gene]].copy()
        frame[time_col] = pd.to_numeric(frame[time_col], errors="coerce")
        frame[gene] = pd.to_numeric(frame[gene], errors="coerce")
        frame = frame.dropna(subset=[status_col, time_col, gene])
        if frame.empty:
            results.append({"gene": gene, "status": "empty_after_cleaning", "p_value": np.nan})
            continue

        event = frame[status_col].eq(event_value).astype(int)
        median_value = frame[gene].median(skipna=True)
        group = np.where(frame[gene] > median_value, high_label, low_label)
        if pd.Series(group).nunique() < 2:
            results.append(
                {
                    "gene": gene,
                    "status": "single_group",
                    "median": median_value,
                    "p_value": np.nan,
                }
            )
            continue

        p_value = logrank_p_value(
            time=frame[time_col].to_numpy(dtype=float),
            event=event.to_numpy(dtype=int),
            group=group,
            high_label=high_label,
        )

        output_file = None
        if plot and output_path is not None:
            output_file = output_path / f"survival_curve_{gene}.png"
            plot_survival_curve(
                time=frame[time_col],
                event=event,
                group=group,
                title=f"Survival Curve for {gene}",
                p_value=p_value,
                output_file=output_file,
                high_label=high_label,
                low_label=low_label,
                risk_table=risk_table,
            )

        results.append(
            {
                "gene": gene,
                "status": "ok",
                "median": median_value,
                "p_value": p_value,
                "n_high": int((group == high_label).sum()),
                "n_low": int((group == low_label).sum()),
                "output_file": str(output_file) if output_file else None,
            }
        )

    return pd.DataFrame(results)


def kaplan_meier_curve(time: Sequence[float], event: Sequence[int]) -> pd.DataFrame:
    """Return Kaplan-Meier step points for a single group."""

    frame = pd.DataFrame({"time": pd.to_numeric(pd.Series(time), errors="coerce"), "event": event})
    frame = frame.dropna(subset=["time"])
    frame["event"] = frame["event"].astype(int)
    event_times = np.sort(frame.loc[frame["event"].eq(1), "time"].unique())
    survival = 1.0
    rows = [{"time": 0.0, "survival": survival, "at_risk": len(frame), "events": 0}]
    for event_time in event_times:
        at_risk = int((frame["time"] >= event_time).sum())
        events = int(((frame["time"] == event_time) & frame["event"].eq(1)).sum())
        if at_risk > 0:
            survival *= 1 - events / at_risk
        rows.append(
            {
                "time": float(event_time),
                "survival": survival,
                "at_risk": at_risk,
                "events": events,
            }
        )
    return pd.DataFrame(rows)


def logrank_p_value(
    time: Sequence[float],
    event: Sequence[int],
    group: Sequence[str],
    high_label: str = "High",
) -> float:
    """Two-group log-rank p-value with a 1-df chi-square survival function."""

    frame = pd.DataFrame({"time": time, "event": event, "group": group}).dropna(subset=["time"])
    frame["event"] = frame["event"].astype(int)
    event_times = np.sort(frame.loc[frame["event"].eq(1), "time"].unique())
    observed = 0.0
    expected = 0.0
    variance = 0.0
    for event_time in event_times:
        at_risk = frame["time"] >= event_time
        events_at_time = frame["time"].eq(event_time) & frame["event"].eq(1)
        n = int(at_risk.sum())
        n1 = int((at_risk & frame["group"].eq(high_label)).sum())
        n2 = n - n1
        d = int(events_at_time.sum())
        d1 = int((events_at_time & frame["group"].eq(high_label)).sum())
        if n <= 1 or d == 0:
            continue
        observed += d1
        expected += d * n1 / n
        variance += (n1 * n2 * d * (n - d)) / (n * n * (n - 1))
    if variance <= 0:
        return float("nan")
    chi_square = (observed - expected) ** 2 / variance
    return float(math.erfc(math.sqrt(chi_square / 2)))


def plot_survival_curve(
    time: Sequence[float],
    event: Sequence[int],
    group: Sequence[str],
    title: str,
    p_value: float,
    output_file: str | Path,
    high_label: str = "High",
    low_label: str = "Low",
    risk_table: bool = True,
) -> None:
    """Save a Kaplan-Meier plot, optionally with a compact risk table."""

    frame = pd.DataFrame({"time": time, "event": event, "group": group}).dropna(subset=["time"])
    groups = [high_label, low_label]
    if risk_table:
        fig, (ax, risk_ax) = plt.subplots(
            2,
            1,
            figsize=(8, 6),
            gridspec_kw={"height_ratios": [4, 1]},
            sharex=True,
        )
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        risk_ax = None

    for label in groups:
        subset = frame.loc[frame["group"].eq(label)]
        km = kaplan_meier_curve(subset["time"], subset["event"])
        ax.step(km["time"], km["survival"], where="post", label=f"{label} Expression")

    ax.set_title(title)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.05)
    ax.legend()
    p_text = "p = NA" if pd.isna(p_value) else f"p = {p_value:.4g}"
    ax.text(0.03, 0.08, p_text, transform=ax.transAxes)

    if risk_ax is not None:
        ticks = _risk_table_ticks(frame["time"])
        risk_ax.set_yticks(range(len(groups)))
        risk_ax.set_yticklabels([f"{label} at risk" for label in groups])
        risk_ax.set_xlabel("Time (months)")
        for row, label in enumerate(groups):
            subset = frame.loc[frame["group"].eq(label)]
            counts = [int((subset["time"] >= tick).sum()) for tick in ticks]
            for tick, count in zip(ticks, counts, strict=True):
                risk_ax.text(tick, row, str(count), ha="center", va="center")
        risk_ax.set_xticks(ticks)
        risk_ax.set_ylim(-0.5, len(groups) - 0.5)
        risk_ax.grid(False)
        for spine in risk_ax.spines.values():
            spine.set_visible(False)

    fig.tight_layout()
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def _biomart_query_xml(gene_ids: Sequence[str]) -> str:
    values = ",".join(gene_ids)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="ensembl_gene_id" value="{values}"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="hgnc_symbol"/>
  </Dataset>
</Query>"""


def _risk_table_ticks(time: pd.Series, ticks: int = 5) -> list[float]:
    max_time = float(pd.to_numeric(time, errors="coerce").max())
    if not math.isfinite(max_time) or max_time <= 0:
        return [0.0]
    return [float(value) for value in np.linspace(0, max_time, ticks)]


def _extract_receptor_status(follow_ups: Sequence[dict[str, object]]) -> dict[str, str | None]:
    status: dict[str, str | None] = {"ESR1": None, "PGR": None, "ERBB2": None}
    for follow_up in follow_ups:
        for test in _as_list_of_dicts(follow_up.get("molecular_tests")):
            if str(test.get("molecular_analysis_method") or "").upper() != "IHC":
                continue
            gene_symbol = str(test.get("gene_symbol") or "")
            if gene_symbol not in status or status[gene_symbol] is not None:
                continue
            result = _normalize_status(test.get("test_result"))
            if result is not None:
                status[gene_symbol] = result
    return status


def _derive_os(
    demographic: dict[str, object],
    diagnoses: Sequence[dict[str, object]],
    follow_ups: Sequence[dict[str, object]],
) -> tuple[str | None, float | None]:
    vital_status = str(demographic.get("vital_status") or "")
    if vital_status.lower() == "dead":
        status = "1:DECEASED"
        days = _to_float(demographic.get("days_to_death"))
    elif vital_status.lower() == "alive":
        status = "0:LIVING"
        days = None
    else:
        status = None
        days = None

    if days is None:
        days = _max_not_none(
            [
                *(_to_float(diagnosis.get("days_to_last_follow_up")) for diagnosis in diagnoses),
                *(_to_float(follow_up.get("days_to_follow_up")) for follow_up in follow_ups),
            ]
        )
    return status, _days_to_months(days)


def _derive_dfs(
    diagnoses: Sequence[dict[str, object]],
    follow_ups: Sequence[dict[str, object]],
) -> tuple[str | None, float | None]:
    event_days: list[float] = []
    for follow_up in follow_ups:
        progression = str(follow_up.get("progression_or_recurrence") or "").lower()
        event_candidate = progression == "yes"
        recurrence_day = _to_float(follow_up.get("days_to_recurrence"))
        progression_day = _to_float(follow_up.get("days_to_progression"))
        if recurrence_day is not None or progression_day is not None:
            event_candidate = True
        if event_candidate:
            event_day = _min_not_none(
                [
                    recurrence_day,
                    progression_day,
                    _to_float(follow_up.get("days_to_follow_up")),
                ]
            )
            if event_day is not None:
                event_days.append(event_day)

    if event_days:
        return "1:Recurred/Progressed", _days_to_months(min(event_days))

    censor_day = _max_not_none(
        [
            *(_to_float(diagnosis.get("days_to_last_follow_up")) for diagnosis in diagnoses),
            *(_to_float(follow_up.get("days_to_follow_up")) for follow_up in follow_ups),
        ]
    )
    if censor_day is None:
        return None, None
    return "0:DiseaseFree", _days_to_months(censor_day)


def _normalize_status(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in {"not reported", "unknown", "", "nan"}:
        return None
    if text.lower() == "positive":
        return "Positive"
    if text.lower() == "negative":
        return "Negative"
    if text.lower() == "equivocal":
        return "Equivocal"
    return text


def _to_float(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _days_to_months(days: float | None) -> float | None:
    if days is None:
        return None
    return days / 30.4375


def _min_not_none(values: Iterable[float | None]) -> float | None:
    clean = [value for value in values if value is not None]
    return min(clean) if clean else None


def _max_not_none(values: Iterable[float | None]) -> float | None:
    clean = [value for value in values if value is not None]
    return max(clean) if clean else None


def _as_list_of_dicts(value: object) -> list[dict[str, object]]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, dict)]


def _require_columns(frame: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(f"Missing required column(s): {', '.join(missing)}")
