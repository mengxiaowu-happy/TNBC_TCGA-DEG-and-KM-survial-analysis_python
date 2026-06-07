from __future__ import annotations

import pandas as pd

from tools.tools import (
    analyze_gene_survival,
    build_query_from_gdc_directory,
    build_survival_expression_data,
    convert_ensembl_to_hgnc,
    discover_star_counts_files,
    drop_duplicate_sample_columns,
    filter_tnbc_samples,
    gene_info_from_star_row_data,
    infer_gdc_download_root,
    log2_transform,
    match_expression_to_clinical,
    parse_gdc_brca_clinical,
    pivot_cbioportal_clinical,
    simplify_tcga_sample_ids,
    strip_ensembl_versions,
)


def test_tnbc_expression_cleaning_pipeline():
    clinical = pd.DataFrame(
        {
            "patientId": ["TCGA-01-0001", "TCGA-01-0002", None, "TCGA-01-0004"],
            "ER_STATUS_BY_IHC": ["Negative", "Positive", "Negative", "Negative"],
            "PR_STATUS_BY_IHC": ["Negative", "Negative", "Negative", "Negative"],
            "IHC_HER2": ["Negative", "Negative", "Negative", "Positive"],
            "DFS_STATUS": ["1:Recurred/Progressed", "0:DiseaseFree", "0:DiseaseFree", "0:DiseaseFree"],
            "DFS_MONTHS": [10, 12, 20, 30],
        }
    )
    expression = pd.DataFrame(
        {
            "TCGA-01-0001-01A-11R": [4, 8],
            "TCGA-01-0001-02A-11R": [5, 9],
            "TCGA-01-0002-01A-11R": [6, 10],
        },
        index=["ENSG000001.1", "ENSG000002.2"],
    )

    clin_data = filter_tnbc_samples(clinical)
    assert clin_data["patientId"].tolist() == ["TCGA-01-0001"]

    simplified = simplify_tcga_sample_ids(expression)
    matched = match_expression_to_clinical(simplified, clin_data)
    deduped = drop_duplicate_sample_columns(matched)
    stripped = strip_ensembl_versions(deduped)

    assert deduped.shape == (2, 1)
    assert stripped.index.tolist() == ["ENSG000001", "ENSG000002"]


def test_gene_mapping_log_and_survival_analysis(tmp_path):
    expression = pd.DataFrame(
        {
            "TCGA-01-0001": [0, 8, 2],
            "TCGA-01-0002": [4, 8, 5],
            "TCGA-01-0003": [16, 8, 9],
            "TCGA-01-0004": [32, 8, 12],
        },
        index=["ENSG1", "ENSG2", "ENSG3"],
    )
    gene_info = pd.DataFrame(
        {
            "ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "hgnc_symbol": ["GENE1", "GENE1", "GENE3"],
        }
    )
    converted = convert_ensembl_to_hgnc(expression, gene_info)
    assert converted.index.tolist() == ["GENE1", "GENE3"]

    logged = log2_transform(converted)
    assert logged.loc["GENE1", "TCGA-01-0001"] == 0

    clinical = pd.DataFrame(
        {
            "patientId": ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003", "TCGA-01-0004"],
            "DFS_STATUS": [
                "0:DiseaseFree",
                "1:Recurred/Progressed",
                "0:DiseaseFree",
                "1:Recurred/Progressed",
            ],
            "DFS_MONTHS": [5, 10, 15, 20],
        }
    )
    merged = build_survival_expression_data(clinical, logged, ["GENE1", "GENE3"])
    results = analyze_gene_survival(merged, ["GENE1", "MISSING"], output_dir=tmp_path)

    ok = results.loc[results["gene"].eq("GENE1")].iloc[0]
    missing = results.loc[results["gene"].eq("MISSING")].iloc[0]
    assert ok["status"] == "ok"
    assert pd.notna(ok["p_value"])
    assert (tmp_path / "survival_curve_GENE1.png").exists()
    assert missing["status"] == "missing_gene"


def test_discover_downloaded_gdc_format_and_build_query(tmp_path, monkeypatch):
    file_dir = (
        tmp_path
        / "GDCdata"
        / "TCGA-BRCA"
        / "Transcriptome_Profiling"
        / "Gene_Expression_Quantification"
        / "file-1"
    )
    file_dir.mkdir(parents=True)
    counts_file = file_dir / "sample.rna_seq.augmented_star_gene_counts.tsv"
    counts_file.write_text("gene_id\tgene_name\tgene_type\tunstranded\n", encoding="utf-8")

    discovered = discover_star_counts_files(tmp_path / "GDCdata")
    assert discovered.loc[0, "file_id"] == "file-1"
    assert discovered.loc[0, "file_name"] == counts_file.name
    assert infer_gdc_download_root(file_dir, project="TCGA-BRCA") == tmp_path / "GDCdata"

    def fake_metadata(file_ids, project="TCGA-BRCA", chunk_size=100):
        return pd.DataFrame(
            {
                "file_id": ["file-1"],
                "file_name": [counts_file.name],
                "md5sum": ["abc"],
                "file_size": [10],
                "state": ["released"],
                "project": ["TCGA-BRCA"],
                "data_category": ["Transcriptome Profiling"],
                "data_type": ["Gene Expression Quantification"],
                "analysis_workflow_type": ["STAR - Counts"],
                "cases": ["TCGA-AA-0001-01A-11R-0000-07"],
                "sample_type": ["Primary Tumor"],
            }
        )

    monkeypatch.setattr("tools.tools.fetch_gdc_file_metadata", fake_metadata)
    query = build_query_from_gdc_directory(tmp_path / "GDCdata")
    assert query.results.loc[0, "cases"] == "TCGA-AA-0001-01A-11R-0000-07"
    assert query.workflow_type == "STAR - Counts"


def test_pivot_cbioportal_and_star_gene_map():
    long_frame = pd.DataFrame(
        {
            "patientId": ["P1", "P1", "P2"],
            "clinicalAttributeId": ["DFS_STATUS", "DFS_MONTHS", "DFS_STATUS"],
            "value": ["0:DiseaseFree", "12", "1:Recurred/Progressed"],
        }
    )
    wide = pivot_cbioportal_clinical(long_frame)
    assert wide.loc[wide["patientId"].eq("P1"), "DFS_MONTHS"].iloc[0] == "12"

    row_data = pd.DataFrame(
        {
            "gene_id": ["ENSG1.1", "ENSG2.2"],
            "gene_name": ["GENE1", "GENE2"],
        }
    )
    gene_info = gene_info_from_star_row_data(row_data)
    assert gene_info.to_dict("records")[0] == {
        "ensembl_gene_id": "ENSG1",
        "hgnc_symbol": "GENE1",
    }


def test_parse_gdc_brca_clinical_derives_receptors_and_survival():
    cases = [
        {
            "submitter_id": "TCGA-AA-0001",
            "case_id": "case-1",
            "demographic": {"vital_status": "Dead", "days_to_death": 365},
            "diagnoses": [{"days_to_last_follow_up": 200}],
            "follow_ups": [
                {
                    "days_to_follow_up": 100,
                    "molecular_tests": [
                        {"molecular_analysis_method": "IHC", "gene_symbol": "ESR1", "test_result": "Negative"},
                        {"molecular_analysis_method": "IHC", "gene_symbol": "PGR", "test_result": "Negative"},
                        {"molecular_analysis_method": "IHC", "gene_symbol": "ERBB2", "test_result": "Negative"},
                    ],
                },
                {
                    "progression_or_recurrence": "Yes",
                    "days_to_progression": 120,
                    "days_to_follow_up": 150,
                },
            ],
        }
    ]

    clinical = parse_gdc_brca_clinical(cases)
    row = clinical.iloc[0]
    assert row["patientId"] == "TCGA-AA-0001"
    assert row["ER_STATUS_BY_IHC"] == "Negative"
    assert row["PR_STATUS_BY_IHC"] == "Negative"
    assert row["IHC_HER2"] == "Negative"
    assert row["OS_STATUS"] == "1:DECEASED"
    assert row["DFS_STATUS"] == "1:Recurred/Progressed"
    assert round(row["DFS_MONTHS"], 2) == round(120 / 30.4375, 2)
