"""Download and prepare the full TCGA-BRCA STAR Counts expression matrix.

Run from the repository root:

    uv run python examples/download_brca_star_counts.py

If GDC is blocked by your network, run:

    uv run python examples/download_brca_star_counts.py --proxy
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from ptcgabiolinks import GDCdownload, GDCprepare, GDCquery, assay, getManifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download all TCGA-BRCA STAR Counts files and save the unstranded expression matrix."
    )
    parser.add_argument(
        "--directory",
        default="GDCdata",
        help="Directory for downloaded GDC files. Default: GDCdata",
    )
    parser.add_argument(
        "--files-per-chunk",
        type=int,
        default=10,
        help="Number of files per GDC API download chunk. Default: 10",
    )
    parser.add_argument(
        "--output",
        default="BRCA_expre_matrix.pkl",
        help="Pickle file for assay(data), equivalent to the R expre_matrix object.",
    )
    parser.add_argument(
        "--proxy",
        action="store_true",
        help="Use local proxy http://127.0.0.1:7890 and socks5://127.0.0.1:7890.",
    )
    parser.add_argument(
        "--skip-prepare",
        action="store_true",
        help="Only download files; skip GDCprepare and matrix serialization.",
    )
    return parser.parse_args()


def configure_proxy(enabled: bool) -> None:
    if not enabled:
        return
    os.environ.setdefault("https_proxy", "http://127.0.0.1:7890")
    os.environ.setdefault("http_proxy", "http://127.0.0.1:7890")
    os.environ.setdefault("all_proxy", "socks5://127.0.0.1:7890")


def main() -> None:
    args = parse_args()
    configure_proxy(args.proxy)

    print("Creating GDC query for TCGA-BRCA STAR - Counts...")
    query = GDCquery(
        project="TCGA-BRCA",
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
    )

    manifest = getManifest(query)
    total_bytes = int(manifest["size"].astype(int).sum())
    print(f"Files to download: {len(manifest)}")
    print(f"Total download size: {total_bytes / 1e9:.2f} GB")
    print(f"Download directory: {Path(args.directory).resolve()}")

    GDCdownload(
        query,
        files_per_chunk=args.files_per_chunk,
        directory=args.directory,
    )

    if args.skip_prepare:
        print("Download finished. Skipping GDCprepare because --skip-prepare was set.")
        return

    print("Preparing SummarizedExperiment-like object...")
    data = GDCprepare(query, directory=args.directory)
    expre_matrix = assay(data)

    output = Path(args.output)
    expre_matrix.to_pickle(output)
    print(f"Expression matrix shape: {expre_matrix.shape[0]} genes x {expre_matrix.shape[1]} samples")
    print(f"Saved expression matrix: {output.resolve()}")


if __name__ == "__main__":
    main()
