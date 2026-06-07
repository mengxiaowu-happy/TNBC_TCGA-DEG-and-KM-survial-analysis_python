"""Lightweight data objects used by the Python port."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import pandas as pd


@dataclass(slots=True)
class GDCQuery:
    """Result returned by :func:`GDCquery`.

    The shape intentionally follows the R object: query parameters are retained
    alongside a results table that can be passed to download and prepare.
    """

    results: pd.DataFrame
    project: list[str]
    data_category: str
    data_type: str | None = None
    workflow_type: str | None = None
    access: list[str] | None = None
    experimental_strategy: list[str] | None = None
    platform: list[str] | None = None
    sample_type: list[str] | None = None
    barcode: list[str] | None = None
    data_format: str | None = None

    def __getitem__(self, key: str) -> Any:
        aliases = {
            "data.category": "data_category",
            "data.type": "data_type",
            "workflow.type": "workflow_type",
            "experimental.strategy": "experimental_strategy",
            "sample.type": "sample_type",
            "data.format": "data_format",
        }
        if key == "results":
            return [self.results]
        return getattr(self, aliases.get(key, key))

    def copy_with_results(self, results: pd.DataFrame) -> "GDCQuery":
        return GDCQuery(
            results=results,
            project=self.project.copy(),
            data_category=self.data_category,
            data_type=self.data_type,
            workflow_type=self.workflow_type,
            access=self.access.copy() if self.access else None,
            experimental_strategy=(
                self.experimental_strategy.copy()
                if self.experimental_strategy
                else None
            ),
            platform=self.platform.copy() if self.platform else None,
            sample_type=self.sample_type.copy() if self.sample_type else None,
            barcode=self.barcode.copy() if self.barcode else None,
            data_format=self.data_format,
        )


@dataclass(slots=True)
class SummarizedExperiment:
    """Small Python analogue of Bioconductor's SummarizedExperiment."""

    assays: dict[str, pd.DataFrame]
    row_data: pd.DataFrame = field(default_factory=pd.DataFrame)
    col_data: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int]:
        first = next(iter(self.assays.values()))
        return first.shape

    @property
    def index(self) -> pd.Index:
        return next(iter(self.assays.values())).index

    @property
    def columns(self) -> pd.Index:
        return next(iter(self.assays.values())).columns

    def assay(self, name: str | None = None) -> pd.DataFrame:
        if name is None:
            name = "unstranded" if "unstranded" in self.assays else next(iter(self.assays))
        try:
            return self.assays[name]
        except KeyError as exc:
            available = ", ".join(self.assays)
            raise KeyError(f"Unknown assay {name!r}. Available assays: {available}") from exc

    def __getitem__(self, item: str | Iterable[str]) -> pd.Series | pd.DataFrame:
        return self.col_data[item]
