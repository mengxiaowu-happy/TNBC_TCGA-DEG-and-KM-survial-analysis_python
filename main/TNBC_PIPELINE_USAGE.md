# TCGA-BRCA 三阴性乳腺癌分析脚本使用说明

本文档面向刚接触 Python 和生信数据分析的开发者，说明如何运行脚本、脚本会下载什么数据、会生成哪些结果。

## 1. 脚本位置

主脚本在：

```bash
experiments/tnbc_survival_workflow.py
```

虽然文件名里仍然保留了 `survival_workflow`，但现在脚本包含两大分析：

1. TNBC tumor vs Normal tissue 差异表达分析。
2. TNBC 患者内部的 KM 生存曲线分析。

## 2. 运行方式

在项目根目录执行：

```bash
uv run python experiments/tnbc_survival_workflow.py
```

脚本不需要任何命令行参数，也不需要手动传入 CSV 文件。

## 3. 脚本默认会做什么

脚本固定分析 TCGA-BRCA 数据。

它会自动完成：

1. 从 GDC 下载完整 TCGA-BRCA 临床病例数据。
2. 从 GDC 下载 TCGA-BRCA STAR Counts 表达量数据。
3. 解析临床数据中的 ER、PR、HER2 状态。
4. 筛选三阴性乳腺癌患者。
5. 根据 TCGA barcode 识别正常组织样本。
6. 构建差异表达分析 metadata。
7. 执行 TNBC tumor vs Normal tissue 差异表达分析。
8. 绘制火山图。
9. 匹配 TNBC 临床样本和表达矩阵，并对匹配后的基因做 TNBC 患者内部 DFS Kaplan-Meier 生存分析。

## 4. 下载目录

表达量数据会下载到：

```bash
GDCdata/
```

典型文件结构如下：

```text
GDCdata/
  TCGA-BRCA/
    Transcriptome_Profiling/
      Gene_Expression_Quantification/
        <file_id>/
          *.rna_seq.augmented_star_gene_counts.tsv
```

如果文件已经存在，下载函数会自动跳过，不会重复下载同一个文件。

## 5. 输出目录

所有分析结果默认写入：

```bash
TNBC_output/
```

主要输出文件如下：

| 文件 | 含义 |
|---|---|
| `gdc_clinical_cases_raw.json` | GDC 原始临床病例 JSON |
| `clinical_gdc.csv` | 整理后的一行一个患者的临床表 |
| `TNBC_clin.pkl` | 筛选出的 TNBC 临床样本 |
| `BRCA_raw_counts.pkl` | BRCA 原始 unstranded counts 表达矩阵 |
| `deg_metadata_tnbc_vs_normal.csv` | DEG 分析使用的样本分组表 |
| `star_gene_name_map.csv` | STAR Counts 中 Ensembl ID 到 gene symbol 的映射 |
| `TNBC_vs_Normal_deg_results.csv` | TNBC vs Normal 差异表达结果 |
| `TNBC_vs_Normal_volcano.png` | 差异表达火山图 |
| `TNBC_upregulated_genes.csv` | TNBC 中上调的显著基因 |
| `TNBC_downregulated_genes.csv` | TNBC 中下调的显著基因 |
| `TNBC_gene_exp_mat.pkl` | TNBC 患者生存分析使用的 log2 表达矩阵 |
| `TNBC_matched_gene_list.csv` | TNBC 临床样本和表达矩阵匹配后得到的基因列表 |
| `merged_survival_expression.csv` | TNBC 生存字段与匹配后基因表达合并后的数据 |
| `significant_genes.csv` | 每个匹配基因的 KM/log-rank 结果 |
| `survival_curve_<GENE>.png` | 每个匹配基因的 KM 生存曲线 |

## 6. 差异表达分析的分组规则

差异表达分析只比较两组：

```text
TNBC tumor samples vs Normal tissue samples
```

其中：

### TNBC tumor samples

必须同时满足：

```text
barcode sample type code == 01
ER_STATUS_BY_IHC == Negative
PR_STATUS_BY_IHC == Negative
IHC_HER2 == Negative
```

### Normal tissue samples

必须满足：

```text
barcode sample type code == 11
```

TCGA barcode 中：

- `01` 表示 Primary Solid Tumor。
- `11` 表示 Solid Tissue Normal。

非 TNBC 肿瘤样本不会进入 DEG 分析。

## 7. log2FoldChange 怎么解释

脚本固定比较方向为：

```text
TNBC vs Normal
```

所以：

```text
log2FoldChange = log2(TNBC / Normal)
```

解释规则：

- `log2FoldChange > 0`：基因在 TNBC 中表达更高。
- `log2FoldChange < 0`：基因在 Normal 中表达更高，或者说在 TNBC 中下调。

## 8. 火山图怎么看

火山图文件：

```bash
TNBC_output/TNBC_vs_Normal_volcano.png
```

图中：

- 横轴：`log2FoldChange = log2(TNBC / Normal)`
- 纵轴：`-log10(adjusted p-value)`
- 红色：TNBC 中上调的显著基因
- 蓝色：TNBC 中下调的显著基因
- 灰色：未达到显著阈值的基因

默认显著阈值：

```text
padj < 0.05
abs(log2FoldChange) >= 1
```

## 9. KM 生存分析做什么

KM 生存分析只在 TNBC 患者内部进行。

脚本不会把 Normal tissue 样本放进 KM 生存分析。

脚本不再使用预设基因列表。它会先把 TNBC 临床表中的 `patientId` 和表达矩阵列名匹配起来，再把匹配后的表达矩阵行名作为 KM 分析的基因列表。

每个匹配基因会按表达量中位数分成：

```text
High expression
Low expression
```

然后计算 DFS 生存曲线和 log-rank p 值。

## 10. 常见问题

### 为什么临床样本有一千多个，但 TNBC 只有一百多个？

GDC 下载的是完整 TCGA-BRCA 临床病例。TNBC 是其中 ER、PR、HER2 全阴性的子集，所以数量会明显减少。

### 为什么表达文件数和临床病例数不同？

临床数据是患者级别，表达数据是样本/文件级别。一个患者可能有多个样本，也可能没有对应表达文件。

### 这个脚本是不是完整复刻 DESeq2？

不是。脚本使用 raw counts，并实现了轻量级 median-of-ratios 样本归一化、log2 fold change、近似 Wald 检验和 BH 校正。它保持了 DESeq2 常用输出字段和比较方向，但不是完整 PyDESeq2/DESeq2 模型。
