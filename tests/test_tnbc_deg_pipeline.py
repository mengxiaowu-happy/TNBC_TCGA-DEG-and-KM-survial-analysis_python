from __future__ import annotations

import pandas as pd

from experiments.tnbc_survival_workflow import (
    build_deg_config,
    build_deg_metadata,
    filter_downregulated_genes,
    filter_upregulated_genes,
    match_tnbc_expression_and_gene_list,
    parse_tcga_sample_type,
    run_differential_expression,
    select_tnbc_patients,
)


def test_metadata_only_contains_tnbc_and_normal():
    expression = pd.DataFrame(
        {
            "TCGA-XX-0001-01A": [100, 200],
            "TCGA-XX-0002-11A": [20, 30],
            "TCGA-XX-0003-01A": [50, 60],
        },
        index=["GeneA", "GeneB"],
    )
    clinical = pd.DataFrame(
        [
            {
                "patient_id": "TCGA-XX-0001",
                "er_status": "Negative",
                "pr_status": "Negative",
                "her2_status": "Negative",
            },
            {
                "patient_id": "TCGA-XX-0003",
                "er_status": "Positive",
                "pr_status": "Negative",
                "her2_status": "Negative",
            },
        ]
    )

    metadata = build_deg_metadata(expression, clinical)

    assert set(metadata["condition"].unique()) == {"TNBC", "Normal"}


def test_non_tnbc_tumor_excluded_from_deg():
    expression = pd.DataFrame(
        {
            "TCGA-XX-0001-01A": [100, 200],
            "TCGA-XX-0002-11A": [20, 30],
            "TCGA-XX-0003-01A": [50, 60],
        },
        index=["GeneA", "GeneB"],
    )
    clinical = pd.DataFrame(
        [
            {
                "patient_id": "TCGA-XX-0001",
                "er_status": "Negative",
                "pr_status": "Negative",
                "her2_status": "Negative",
            },
            {
                "patient_id": "TCGA-XX-0003",
                "er_status": "Positive",
                "pr_status": "Negative",
                "her2_status": "Negative",
            },
        ]
    )

    metadata = build_deg_metadata(expression, clinical)
    non_tnbc_tumor_patients = clinical[
        ~(
            clinical["er_status"].eq("Negative")
            & clinical["pr_status"].eq("Negative")
            & clinical["her2_status"].eq("Negative")
        )
    ]["patient_id"]

    assert not metadata["patient_id"].isin(non_tnbc_tumor_patients).any()


def test_parse_primary_tumor_sample():
    sample_type = parse_tcga_sample_type("TCGA-A1-A0SB-01A")
    assert sample_type == "tumor"


def test_parse_normal_sample():
    sample_type = parse_tcga_sample_type("TCGA-A1-A0SB-11A")
    assert sample_type == "normal"


def test_select_tnbc_patient():
    clinical = pd.DataFrame(
        [
            {
                "patient_id": "TCGA-XX-0001",
                "er_status": "Negative",
                "pr_status": "Negative",
                "her2_status": "Negative",
            }
        ]
    )

    result = select_tnbc_patients(clinical)

    assert result == {"TCGA-XX-0001"}


def test_exclude_er_positive_patient():
    clinical = pd.DataFrame(
        [
            {
                "patient_id": "TCGA-XX-0002",
                "er_status": "Positive",
                "pr_status": "Negative",
                "her2_status": "Negative",
            }
        ]
    )

    result = select_tnbc_patients(clinical)

    assert "TCGA-XX-0002" not in result


def test_exclude_unknown_receptor_status():
    clinical = pd.DataFrame(
        [
            {
                "patient_id": "TCGA-XX-0003",
                "er_status": "Negative",
                "pr_status": "Negative",
                "her2_status": "Unknown",
            }
        ]
    )

    result = select_tnbc_patients(clinical)

    assert "TCGA-XX-0003" not in result


def test_deg_contrast_is_tnbc_vs_normal():
    deg_config = build_deg_config(case_group="TNBC", control_group="Normal")

    assert deg_config.case_group == "TNBC"
    assert deg_config.control_group == "Normal"
    assert deg_config.contrast_name == "TNBC_vs_Normal"


def test_deg_metadata_condition_order():
    expression = pd.DataFrame(
        {
            "TCGA-XX-0001-01A": [100, 200],
            "TCGA-XX-0002-11A": [20, 30],
        },
        index=["GeneA", "GeneB"],
    )
    clinical = pd.DataFrame(
        [
            {
                "patient_id": "TCGA-XX-0001",
                "er_status": "Negative",
                "pr_status": "Negative",
                "her2_status": "Negative",
            }
        ]
    )

    metadata = build_deg_metadata(expression, clinical)

    assert "TNBC" in metadata["condition"].values
    assert "Normal" in metadata["condition"].values


def test_log2fc_direction_tnbc_over_normal():
    deg_result = pd.DataFrame(
        [
            {
                "gene_id": "GeneA",
                "log2FoldChange": 2.0,
                "padj": 0.001,
            },
            {
                "gene_id": "GeneB",
                "log2FoldChange": -1.5,
                "padj": 0.002,
            },
        ]
    )

    up_genes = filter_upregulated_genes(deg_result)
    down_genes = filter_downregulated_genes(deg_result)

    assert "GeneA" in set(up_genes["gene_id"])
    assert "GeneB" in set(down_genes["gene_id"])


def test_run_differential_expression_outputs_tnbc_over_normal_direction():
    counts = pd.DataFrame(
        {
            "TCGA-XX-0001-01A": [100, 20],
            "TCGA-XX-0004-01A": [120, 25],
            "TCGA-XX-0002-11A": [10, 90],
            "TCGA-XX-0005-11A": [15, 95],
        },
        index=["ENSG1", "ENSG2"],
    )
    metadata = pd.DataFrame(
        {
            "sample_id": [
                "TCGA-XX-0001-01A",
                "TCGA-XX-0004-01A",
                "TCGA-XX-0002-11A",
                "TCGA-XX-0005-11A",
            ],
            "patient_id": [
                "TCGA-XX-0001",
                "TCGA-XX-0004",
                "TCGA-XX-0002",
                "TCGA-XX-0005",
            ],
            "sample_type": ["tumor", "tumor", "normal", "normal"],
            "condition": ["TNBC", "TNBC", "Normal", "Normal"],
        }
    )
    gene_info = pd.DataFrame(
        {
            "ensembl_gene_id": ["ENSG1", "ENSG2"],
            "hgnc_symbol": ["GeneA", "GeneB"],
        }
    )

    result = run_differential_expression(counts, metadata, gene_info, build_deg_config())

    gene_a = result.loc[result["gene_id"].eq("ENSG1")].iloc[0]
    gene_b = result.loc[result["gene_id"].eq("ENSG2")].iloc[0]
    assert gene_a["contrast"] == "TNBC_vs_Normal"
    assert gene_a["log2FoldChange"] > 0
    assert gene_b["log2FoldChange"] < 0


def test_survival_gene_list_comes_from_matched_tnbc_expression():
    log_expression = pd.DataFrame(
        {
            "TCGA-XX-0001": [1.0, 2.0, None],
            "TCGA-XX-9999": [3.0, 4.0, None],
        },
        index=["GeneA", "GeneB", "GeneWithoutMatchedValues"],
    )
    tnbc_clinical = pd.DataFrame(
        [
            {
                "patientId": "TCGA-XX-0001",
                "DFS_STATUS": "0:DiseaseFree",
                "DFS_MONTHS": 10,
            }
        ]
    )

    matched_expression, gene_list = match_tnbc_expression_and_gene_list(log_expression, tnbc_clinical)

    assert list(matched_expression.columns) == ["TCGA-XX-0001"]
    assert gene_list == ["GeneA", "GeneB"]
