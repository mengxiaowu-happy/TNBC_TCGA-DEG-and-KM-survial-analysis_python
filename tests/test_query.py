from __future__ import annotations

from ptcgabiolinks import GDCquery, getManifest, getResults


def test_gdcquery_flattens_and_filters_hits(monkeypatch):
    hit = {
        "file_id": "file-1",
        "file_name": "sample.tsv",
        "md5sum": "abc",
        "file_size": 123,
        "state": "released",
        "data_category": "Transcriptome Profiling",
        "data_type": "Gene Expression Quantification",
        "experimental_strategy": "RNA-Seq",
        "data_format": "TSV",
        "access": "open",
        "analysis": {"analysis_id": "analysis-1", "workflow_type": "STAR - Counts"},
        "cases": [
            {
                "submitter_id": "TCGA-XX-0001",
                "samples": [
                    {
                        "submitter_id": "TCGA-XX-0001-01A",
                        "sample_type": "Primary Tumor",
                        "is_ffpe": False,
                        "portions": [
                            {
                                "analytes": [
                                    {
                                        "aliquots": [
                                            {
                                                "submitter_id": "TCGA-XX-0001-01A-11R-0000-07"
                                            }
                                        ]
                                    }
                                ]
                            }
                        ],
                    }
                ],
            }
        ],
    }

    monkeypatch.setattr("ptcgabiolinks.query.is_serve_ok", lambda: True)
    monkeypatch.setattr(
        "ptcgabiolinks.query._fetch_files_for_project",
        lambda **kwargs: [hit],
    )

    query = GDCquery(
        project="TCGA-TEST",
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
        barcode="TCGA-XX-0001",
    )

    results = getResults(query)
    assert len(results) == 1
    assert results.loc[0, "cases"] == "TCGA-XX-0001-01A-11R-0000-07"
    assert results.loc[0, "sample_type"] == "Primary Tumor"

    manifest = getManifest(query)
    assert manifest.to_dict("records")[0] == {
        "id": "file-1",
        "filename": "sample.tsv",
        "md5": "abc",
        "size": 123,
        "state": "released",
    }
