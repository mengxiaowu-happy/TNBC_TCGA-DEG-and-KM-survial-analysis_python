# PTCGAbiolinks

Python port of the core `TCGAbiolinks` GDC workflow.

```python
from ptcgabiolinks import GDCquery, GDCdownload, GDCprepare, assay

query = GDCquery(
    project="TCGA-BRCA",
    data_category="Transcriptome Profiling",
    data_type="Gene Expression Quantification",
    workflow_type="STAR - Counts",
)

GDCdownload(query, files_per_chunk=10, directory="GDCdata")
data = GDCprepare(query, directory="GDCdata")
expre_matrix = assay(data)
```

The initial port focuses on open GDC API data and the STAR Counts expression
workflow used by TCGA transcriptome profiling.
