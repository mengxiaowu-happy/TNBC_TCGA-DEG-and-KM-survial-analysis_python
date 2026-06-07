from __future__ import annotations

import hashlib
import io
import tarfile

import pandas as pd

from ptcgabiolinks import GDCQuery, GDCdownload


class FakeResponse:
    def __init__(self, content: bytes):
        self.content = content

    def iter_content(self, chunk_size: int):
        for start in range(0, len(self.content), chunk_size):
            yield self.content[start : start + chunk_size]


def test_gdcdownload_api_extracts_tar_layout(tmp_path, monkeypatch):
    files = {
        "file-1/sample-1.tsv": b"gene_id\tvalue\nENSG1\t1\n",
        "file-2/sample-2.tsv": b"gene_id\tvalue\nENSG1\t2\n",
    }
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w:gz") as tar:
        for name, content in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            tar.addfile(info, io.BytesIO(content))

    rows = []
    for idx, (name, content) in enumerate(files.items(), start=1):
        file_id, file_name = name.split("/")
        rows.append(
            {
                "project": "TCGA-TEST",
                "data_category": "Transcriptome Profiling",
                "data_type": "Gene Expression Quantification",
                "file_id": file_id,
                "file_name": file_name,
                "md5sum": hashlib.md5(content).hexdigest(),
                "file_size": len(content),
                "state": "released",
            }
        )

    query = GDCQuery(
        results=pd.DataFrame(rows),
        project=["TCGA-TEST"],
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
    )
    monkeypatch.setattr("ptcgabiolinks.download.is_serve_ok", lambda: True)
    monkeypatch.setattr(
        "ptcgabiolinks.download.post_data",
        lambda ids: FakeResponse(archive.getvalue()),
    )

    GDCdownload(query, directory=tmp_path, files_per_chunk=2)

    root = (
        tmp_path
        / "TCGA-TEST"
        / "Transcriptome_Profiling"
        / "Gene_Expression_Quantification"
    )
    assert (root / "file-1" / "sample-1.tsv").read_bytes() == files["file-1/sample-1.tsv"]
    assert (root / "file-2" / "sample-2.tsv").read_bytes() == files["file-2/sample-2.tsv"]
