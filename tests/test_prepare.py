from __future__ import annotations

import pandas as pd

from ptcgabiolinks import GDCQuery, GDCprepare, assay, assays


STAR_CONTENT = """gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded
N_unmapped\tN_unmapped\tN_unmapped\t1\t2\t3\t0\t0\t0
ENSG00000000001.1\tGENE1\tprotein_coding\t10\t11\t12\t1.5\t2.5\t3.5
ENSG00000000002.1\tGENE2\tprotein_coding\t20\t21\t22\t4.5\t5.5\t6.5
"""


def test_gdcprepare_star_counts(tmp_path, monkeypatch):
    monkeypatch.setattr("ptcgabiolinks.prepare.is_serve_ok", lambda: True)
    monkeypatch.setattr(
        "ptcgabiolinks.prepare.get_gdc_info",
        lambda: {"data_release": "Data Release Test"},
    )

    rows = []
    for idx, case in enumerate(
        ["TCGA-AA-0001-01A-11R-0000-07", "TCGA-AA-0002-01A-11R-0000-07"],
        start=1,
    ):
        file_id = f"file-{idx}"
        file_name = f"sample-{idx}.tsv"
        path = (
            tmp_path
            / "TCGA-TEST"
            / "Transcriptome_Profiling"
            / "Gene_Expression_Quantification"
            / file_id
        )
        path.mkdir(parents=True)
        (path / file_name).write_text(STAR_CONTENT, encoding="utf-8")
        rows.append(
            {
                "project": "TCGA-TEST",
                "data_category": "Transcriptome Profiling",
                "data_type": "Gene Expression Quantification",
                "analysis_workflow_type": "STAR - Counts",
                "file_id": file_id,
                "file_name": file_name,
                "cases": case,
            }
        )

    query = GDCQuery(
        results=pd.DataFrame(rows),
        project=["TCGA-TEST"],
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
    )

    prepared = GDCprepare(query, directory=tmp_path)

    assert set(assays(prepared)) == {
        "unstranded",
        "stranded_first",
        "stranded_second",
        "tpm_unstrand",
        "fpkm_unstrand",
        "fpkm_uq_unstrand",
    }
    assert assay(prepared).loc["ENSG00000000001.1", "TCGA-AA-0001-01A-11R-0000-07"] == 10
    assert assay(prepared, "fpkm_uq_unstrand").loc[
        "ENSG00000000002.1",
        "TCGA-AA-0002-01A-11R-0000-07",
    ] == 6.5
    assert prepared.col_data.loc["TCGA-AA-0001-01A-11R-0000-07", "project_id"] == "TCGA-TEST"
    assert prepared.metadata["data_release"] == "Data Release Test"
