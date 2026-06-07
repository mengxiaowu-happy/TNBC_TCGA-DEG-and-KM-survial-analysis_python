"""Download GDC files into the TCGAbiolinks directory layout."""

from __future__ import annotations

import hashlib
import shutil
import tarfile
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from .client import is_serve_ok, post_data
from .models import GDCQuery


def GDCdownload(
    query: GDCQuery,
    token_file: str | None = None,
    method: str = "api",
    directory: str | Path = "GDCdata",
    files_per_chunk: int | None = None,
    **kwargs: Any,
) -> None:
    """Download files returned by :func:`GDCquery`.

    The API method is implemented for open data. ``method="client"`` delegates
    to an installed ``gdc-client`` binary when present.
    """

    files_per_chunk = _pop_alias(kwargs, files_per_chunk, "files.per.chunk")
    token_file = _pop_alias(kwargs, token_file, "token.file")
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise TypeError(f"Unknown GDCdownload argument(s): {unknown}")
    if query is None:
        raise ValueError("Please set query argument")
    if method not in {"api", "client"}:
        raise ValueError("method arguments possible values are: 'api' or 'client'")
    if token_file and method == "api":
        raise ValueError("token_file is supported only with method='client'")
    data_types = query.results["data_type"].dropna().unique()
    if len(data_types) > 1:
        raise ValueError(
            "We can only download one data type. Please use data_type in GDCquery."
        )

    is_serve_ok()
    base_directory = Path(directory)
    base_directory.mkdir(parents=True, exist_ok=True)

    if method == "client":
        _download_with_client(query, token_file, base_directory)
        return

    for project in query.results["project"].dropna().unique():
        project_results = query.results.loc[query.results["project"] == project]
        manifest = getManifest(query.copy_with_results(project_results))
        path = _download_root(base_directory, project_results)
        manifest = _filter_already_downloaded(path, manifest)
        if manifest.empty:
            print("All samples have been already downloaded")
            continue

        chunks = _chunks(manifest, files_per_chunk or len(manifest))
        total_chunks = len(chunks)
        for idx, chunk in enumerate(chunks, start=1):
            size = _human_readable_bytes(int(chunk["size"].astype(int).sum()))
            print(
                f"Downloading chunk {idx} of {total_chunks} "
                f"({len(chunk)} files, size = {size})"
            )
            _download_manifest_chunk(chunk, path)


def getManifest(query: GDCQuery, save: bool = False) -> pd.DataFrame:
    manifest = query.results.loc[:, ["file_id", "file_name", "md5sum", "file_size", "state"]].copy()
    manifest.columns = ["id", "filename", "md5", "size", "state"]
    if save:
        manifest.to_csv("gdc_manifest.txt", sep="\t", index=False)
    return manifest


def _download_manifest_chunk(manifest: pd.DataFrame, path: Path) -> None:
    ids = manifest["id"].astype(str).tolist()
    response = post_data(ids)
    path.mkdir(parents=True, exist_ok=True)
    if len(ids) == 1:
        row = manifest.iloc[0]
        destination = path / str(row["id"]) / str(row["filename"])
        destination.parent.mkdir(parents=True, exist_ok=True)
        _write_response(response, destination)
        _verify_md5(destination, str(row["md5"]))
        return

    with tempfile.TemporaryDirectory(prefix="ptcgabiolinks_") as tmp:
        tmpdir = Path(tmp)
        archive = tmpdir / "download.tar.gz"
        _write_response(response, archive)
        with tarfile.open(archive, mode="r:*") as tar:
            _safe_extract(tar, tmpdir)
        for _, row in manifest.iterrows():
            source = tmpdir / str(row["id"]) / str(row["filename"])
            destination = path / str(row["id"]) / str(row["filename"])
            if not source.exists():
                raise RuntimeError(f"GDC archive did not contain {source.name}")
            _verify_md5(source, str(row["md5"]))
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), destination)


def _download_with_client(
    query: GDCQuery,
    token_file: str | None,
    directory: Path,
) -> None:
    gdc_client = shutil.which("gdc-client")
    if not gdc_client:
        raise RuntimeError("method='client' requires an installed gdc-client binary")
    manifest = getManifest(query, save=True)
    project_results = query.results
    path = _download_root(directory, project_results)
    path.mkdir(parents=True, exist_ok=True)
    config = Path("gdc_client_configuration.dtt")
    config.write_text(f"[download]\nretry_amount = 6\ndir = {path}\n", encoding="utf-8")
    command = [gdc_client, "download", "-m", "gdc_manifest.txt", "--config", str(config)]
    if token_file:
        command.extend(["-t", token_file])
    import subprocess

    subprocess.run(command, check=True)
    _filter_already_downloaded(path, manifest)


def _download_root(directory: Path, results: pd.DataFrame) -> Path:
    first = results.iloc[0]
    return (
        directory
        / str(first["project"])
        / _space_to_underscore(str(first["data_category"]))
        / _space_to_underscore(str(first["data_type"]))
    )


def _filter_already_downloaded(path: Path, manifest: pd.DataFrame) -> pd.DataFrame:
    missing = []
    for _, row in manifest.iterrows():
        filename = str(row["filename"])
        file_id = str(row["id"])
        missing.append(not ((path / file_id / filename).exists() or (path / filename).exists()))
    return manifest.loc[missing].reset_index(drop=True)


def _write_response(response: Any, destination: Path) -> None:
    with destination.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)


def _verify_md5(path: Path, expected: str) -> None:
    if not expected or expected == "nan":
        return
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    observed = digest.hexdigest()
    if observed.lower() != expected.lower():
        path.unlink(missing_ok=True)
        raise RuntimeError(f"File corrupted: {path}")


def _safe_extract(tar: tarfile.TarFile, path: Path) -> None:
    root = path.resolve()
    for member in tar.getmembers():
        target = (path / member.name).resolve()
        if not str(target).startswith(str(root)):
            raise RuntimeError(f"Refusing to extract unsafe path {member.name!r}")
    try:
        tar.extractall(path, filter="data")
    except TypeError:
        tar.extractall(path)


def _chunks(frame: pd.DataFrame, size: int) -> list[pd.DataFrame]:
    if size <= 0:
        raise ValueError("files_per_chunk must be positive")
    return [frame.iloc[start : start + size].reset_index(drop=True) for start in range(0, len(frame), size)]


def _space_to_underscore(value: str) -> str:
    return value.replace(" ", "_")


def _human_readable_bytes(bytes_count: int) -> str:
    unit = 1000
    if bytes_count < unit:
        return f"{bytes_count} B"
    exp = 0
    value = float(bytes_count)
    prefixes = ["", "KB", "MB", "GB", "TB", "PB"]
    while value >= unit and exp < len(prefixes) - 1:
        value /= unit
        exp += 1
    return f"{value:.1f} {prefixes[exp]}"


def _pop_alias(kwargs: dict[str, Any], value: Any, *aliases: str) -> Any:
    for alias in aliases:
        if alias in kwargs:
            if value is not None:
                raise TypeError(f"Received both {alias!r} and Python alias")
            value = kwargs.pop(alias)
    return value
