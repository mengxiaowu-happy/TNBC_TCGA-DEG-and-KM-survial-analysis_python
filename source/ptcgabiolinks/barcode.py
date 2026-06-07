"""Helpers for sample metadata derived from barcodes."""

from __future__ import annotations

import pandas as pd

TCGA_SAMPLE_TYPES = {
    "01": ("TP", "Primary Tumor"),
    "02": ("TR", "Recurrent Tumor"),
    "03": ("TB", "Primary Blood Derived Cancer - Peripheral Blood"),
    "04": ("TRBM", "Recurrent Blood Derived Cancer - Bone Marrow"),
    "05": ("TAP", "Additional - New Primary"),
    "06": ("TM", "Metastatic"),
    "07": ("TAM", "Additional Metastatic"),
    "08": ("THOC", "Human Tumor Original Cells"),
    "09": ("TBM", "Primary Blood Derived Cancer - Bone Marrow"),
    "10": ("NB", "Blood Derived Normal"),
    "11": ("NT", "Solid Tissue Normal"),
    "12": ("NBC", "Buccal Cell Normal"),
    "13": ("NEBV", "EBV Immortalized Normal"),
    "14": ("NBM", "Bone Marrow Normal"),
    "20": ("CELLC", "Control Analyte"),
    "40": ("TRB", "Recurrent Blood Derived Cancer - Peripheral Blood"),
    "50": ("CELL", "Cell Lines"),
    "60": ("XP", "Primary Xenograft Tissue"),
    "61": ("XCL", "Cell Line Derived Xenograft Tissue"),
}


def col_data_prepare(samples: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample in samples:
        parts = sample.split("-")
        row: dict[str, object] = {
            "barcode": sample,
            "sample_submitter_id": sample,
        }
        if len(parts) >= 4 and parts[0] in {"TCGA", "TARGET"}:
            row["project_id"] = parts[0] if parts[0] == "TARGET" else f"{parts[0]}-{parts[1]}"
            row["patient"] = "-".join(parts[:3])
            row["sample"] = "-".join(parts[:4])
            row["tss"] = parts[1]
            row["participant"] = parts[2]
            tissue_code = parts[3][:2]
            row["tissue_code"] = tissue_code
            short_code, definition = TCGA_SAMPLE_TYPES.get(tissue_code, (None, None))
            row["short_letter_code"] = short_code
            row["sample_type"] = definition
            if len(parts) >= 7:
                row["portion"] = parts[4]
                row["plate"] = parts[5]
                row["center"] = parts[6]
        else:
            row["project_id"] = None
            row["sample"] = sample
        rows.append(row)
    data = pd.DataFrame(rows)
    data.index = samples
    return data
