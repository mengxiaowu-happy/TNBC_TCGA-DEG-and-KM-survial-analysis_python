"""TCGA-BRCA 三阴性乳腺癌差异表达和生存分析流程。


1. 下载并整理 TCGA-BRCA 临床数据。
2. 下载并读取 TCGA-BRCA STAR Counts 表达量数据。
3. 按 ER/PR/HER2 全阴性规则筛选 TNBC 肿瘤样本。
4. 按 TCGA barcode sample type code == 11 识别正常组织样本。
5. 对 TNBC tumor vs Normal tissue 做差异表达分析。
6. 绘制火山图。
7. 在 TNBC 患者内部，对临床样本和表达矩阵成功匹配到的基因做 KM 生存曲线分析。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ptcgabiolinks import GDCdownload, GDCprepare, GDCquery, assay, getManifest
from tools.tools import (
    analyze_gene_survival,
    build_query_from_gdc_directory,
    build_survival_expression_data,
    convert_ensembl_to_hgnc,
    drop_duplicate_sample_columns,
    fetch_gdc_clinical,
    filter_tnbc_samples,
    gene_info_from_star_row_data,
    log2_transform,
    match_expression_to_clinical,
    simplify_tcga_sample_ids,
    strip_ensembl_versions,
)


# ==============================
# 固定配置：脚本不使用命令行参数
# ==============================

PROJECT = "TCGA-BRCA"
GDC_DIRECTORY = Path("GDCdata")
OUTPUT_DIR = Path("TNBC_output")
FILES_PER_CHUNK = 10

CASE_GROUP = "TNBC"
CONTROL_GROUP = "Normal"
VOLCANO_PADJ_THRESHOLD = 0.05
VOLCANO_LOG2FC_THRESHOLD = 1.0


@dataclass(frozen=True)
class DegConfig:
    """差异表达分析的比较方向配置。"""

    case_group: str
    control_group: str
    contrast_name: str


def parse_tcga_sample_type(barcode: str) -> Literal["tumor", "normal", "other"]:
    """根据 TCGA barcode 的 sample type code 判断样本类型。

    TCGA barcode 第 4 段的前两位是 sample type code：
    - 01 表示 Primary Solid Tumor，这里归为 tumor。
    - 11 表示 Solid Tissue Normal，这里归为 normal。
    - 其他编码不参与本脚本的 DEG 分析，归为 other。
    """

    parts = str(barcode).split("-")
    if len(parts) < 4:
        return "other"
    sample_code = parts[3][:2]
    if sample_code == "01":
        return "tumor"
    if sample_code == "11":
        return "normal"
    return "other"


def standardize_clinical_columns(clinical: pd.DataFrame) -> pd.DataFrame:
    """把临床表整理成测试和分组函数使用的标准字段名。"""

    return pd.DataFrame(
        {
            "patient_id": clinical["patientId"].astype(str),
            "er_status": clinical["ER_STATUS_BY_IHC"],
            "pr_status": clinical["PR_STATUS_BY_IHC"],
            "her2_status": clinical["IHC_HER2"],
        }
    )


def select_tnbc_patients(clinical: pd.DataFrame) -> set[str]:
    """筛选 ER、PR、HER2 同时为 Negative 的 TNBC 患者。"""

    required = {"patient_id", "er_status", "pr_status", "her2_status"}
    missing = required - set(clinical.columns)
    if missing:
        raise KeyError(f"临床表缺少字段: {sorted(missing)}")

    mask = (
        clinical["er_status"].eq("Negative")
        & clinical["pr_status"].eq("Negative")
        & clinical["her2_status"].eq("Negative")
        & clinical["patient_id"].notna()
    )
    return set(clinical.loc[mask, "patient_id"].astype(str))


def build_deg_metadata(expression: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """构建 DEG 分析 metadata，只保留 TNBC tumor 和 Normal tissue。

    
    - TNBC 组：barcode 为肿瘤样本，并且该患者满足三阴性临床规则。
    - Normal 组：barcode sample type code 为 11 的 Solid Tissue Normal 样本。
    - 非 TNBC 肿瘤样本必须排除，不能进入 DEG 的 case 组。
    """

    tnbc_patients = select_tnbc_patients(clinical)
    rows: list[dict[str, str]] = []

    for sample_id in expression.columns.astype(str):
        patient_id = sample_id[:12]
        sample_type = parse_tcga_sample_type(sample_id)

        if sample_type == "tumor" and patient_id in tnbc_patients:
            condition = CASE_GROUP
        elif sample_type == "normal":
            condition = CONTROL_GROUP
        else:
            continue

        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": patient_id,
                "sample_type": sample_type,
                "condition": condition,
            }
        )

    metadata = pd.DataFrame(rows)
    if metadata.empty:
        raise RuntimeError("没有找到可用于 DEG 的 TNBC/Normal 样本。")
    return metadata


def build_deg_config(
    case_group: str = CASE_GROUP,
    control_group: str = CONTROL_GROUP,
) -> DegConfig:
    """显式声明 DEG 比较方向。"""

    return DegConfig(
        case_group=case_group,
        control_group=control_group,
        contrast_name=f"{case_group}_vs_{control_group}",
    )


def normalize_counts_median_ratio(counts: pd.DataFrame) -> pd.DataFrame:
    """用简化版 median-of-ratios 方法做样本层面归一化。

    这里的输入是 raw counts。为了降低测序深度差异的影响，我们先估计每个样本的
    size factor，再用 raw counts 除以 size factor。这个思想来自 DESeq2，
    但这里实现的是轻量级版本，便于单文件脚本直接运行。
    """

    numeric_counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    numeric_counts = numeric_counts.clip(lower=0)

    positive = numeric_counts.replace(0, np.nan)
    log_positive = np.log(positive)
    geometric_mean = np.exp(log_positive.mean(axis=1, skipna=True))
    valid_genes = geometric_mean.replace([np.inf, -np.inf], np.nan).dropna()
    valid_genes = valid_genes[valid_genes > 0]

    if valid_genes.empty:
        library_sizes = numeric_counts.sum(axis=0)
        median_library_size = library_sizes.median()
        size_factors = library_sizes / median_library_size
    else:
        ratios = numeric_counts.loc[valid_genes.index].div(valid_genes, axis=0)
        ratios = ratios.replace(0, np.nan)
        size_factors = ratios.median(axis=0, skipna=True)

    size_factors = size_factors.replace([0, np.inf, -np.inf], np.nan)
    size_factors = size_factors.fillna(size_factors.dropna().median())
    size_factors = size_factors.fillna(1.0)
    return numeric_counts.div(size_factors, axis=1)


def run_differential_expression(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_info: pd.DataFrame,
    config: DegConfig,
) -> pd.DataFrame:
    """执行 TNBC vs Normal 的轻量级差异表达分析。

    输出中的 log2FoldChange 固定解释为：

        log2FoldChange = log2(TNBC / Normal)

    因此：
    - log2FoldChange > 0 表示在 TNBC 中上调。
    - log2FoldChange < 0 表示在 Normal 中更高，或在 TNBC 中下调。
    """

    required_metadata = {"sample_id", "condition"}
    missing_metadata = required_metadata - set(metadata.columns)
    if missing_metadata:
        raise KeyError(f"metadata 缺少字段: {sorted(missing_metadata)}")

    metadata = metadata.loc[metadata["sample_id"].isin(counts.columns)].copy()
    if set(metadata["condition"].unique()) != {config.case_group, config.control_group}:
        raise ValueError("DEG metadata 中 condition 必须且只能包含 TNBC 和 Normal。")

    selected_counts = counts.loc[:, metadata["sample_id"]]
    selected_counts = selected_counts.loc[~selected_counts.index.duplicated(keep="first")]
    selected_counts = selected_counts.loc[selected_counts.sum(axis=1) > 0]

    normalized_counts = normalize_counts_median_ratio(selected_counts)
    log_counts = np.log2(normalized_counts + 1)

    case_samples = metadata.loc[metadata["condition"].eq(config.case_group), "sample_id"]
    control_samples = metadata.loc[metadata["condition"].eq(config.control_group), "sample_id"]

    case_values = log_counts.loc[:, case_samples]
    control_values = log_counts.loc[:, control_samples]

    base_mean = normalized_counts.mean(axis=1)
    case_mean = case_values.mean(axis=1)
    control_mean = control_values.mean(axis=1)
    log2_fold_change = case_mean - control_mean

    case_var = case_values.var(axis=1, ddof=1).fillna(0)
    control_var = control_values.var(axis=1, ddof=1).fillna(0)
    lfc_se = np.sqrt(case_var / max(len(case_samples), 1) + control_var / max(len(control_samples), 1))
    lfc_se = lfc_se.replace(0, np.nan)

    stat = log2_fold_change / lfc_se
    pvalue = stat.abs().map(lambda value: math.erfc(value / math.sqrt(2)) if pd.notna(value) else 1.0)
    padj = benjamini_hochberg(pvalue)

    gene_map = gene_info.copy()
    gene_map["ensembl_gene_id"] = gene_map["ensembl_gene_id"].astype(str)
    gene_map = gene_map.drop_duplicates("ensembl_gene_id", keep="first")
    symbol_by_gene = gene_map.set_index("ensembl_gene_id")["hgnc_symbol"]

    result = pd.DataFrame(
        {
            "gene_id": log_counts.index.astype(str),
            "gene_symbol": [symbol_by_gene.get(gene_id, gene_id) for gene_id in log_counts.index.astype(str)],
            "baseMean": base_mean,
            "log2FoldChange": log2_fold_change,
            "lfcSE": lfc_se,
            "stat": stat,
            "pvalue": pvalue,
            "padj": padj,
            "contrast": config.contrast_name,
        }
    )
    return result.sort_values("padj", na_position="last").reset_index(drop=True)


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg 多重检验校正。"""

    p = pd.Series(pvalues, dtype=float).fillna(1.0)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    running_min = 1.0

    for index in range(n - 1, -1, -1):
        rank = index + 1
        value = ranked[index] * n / rank
        running_min = min(running_min, value)
        adjusted[index] = running_min

    output = np.empty(n, dtype=float)
    output[order] = np.clip(adjusted, 0, 1)
    return pd.Series(output, index=p.index)


def filter_upregulated_genes(
    deg_result: pd.DataFrame,
    padj_threshold: float = VOLCANO_PADJ_THRESHOLD,
    log2fc_threshold: float = VOLCANO_LOG2FC_THRESHOLD,
) -> pd.DataFrame:
    """筛选 TNBC 中上调的基因。"""

    return deg_result.loc[
        (deg_result["padj"] < padj_threshold)
        & (deg_result["log2FoldChange"] >= log2fc_threshold)
    ].copy()


def filter_downregulated_genes(
    deg_result: pd.DataFrame,
    padj_threshold: float = VOLCANO_PADJ_THRESHOLD,
    log2fc_threshold: float = VOLCANO_LOG2FC_THRESHOLD,
) -> pd.DataFrame:
    """筛选 TNBC 中下调的基因。"""

    return deg_result.loc[
        (deg_result["padj"] < padj_threshold)
        & (deg_result["log2FoldChange"] <= -log2fc_threshold)
    ].copy()


def plot_volcano(
    deg_result: pd.DataFrame,
    output_file: Path,
    padj_threshold: float = VOLCANO_PADJ_THRESHOLD,
    log2fc_threshold: float = VOLCANO_LOG2FC_THRESHOLD,
) -> None:
    """绘制 TNBC vs Normal 差异表达火山图。"""

    plot_data = deg_result.copy()
    plot_data["padj_for_plot"] = plot_data["padj"].clip(lower=1e-300).fillna(1.0)
    plot_data["minus_log10_padj"] = -np.log10(plot_data["padj_for_plot"])
    plot_data["category"] = "Not significant"
    plot_data.loc[
        (plot_data["padj"] < padj_threshold)
        & (plot_data["log2FoldChange"] >= log2fc_threshold),
        "category",
    ] = "Up in TNBC"
    plot_data.loc[
        (plot_data["padj"] < padj_threshold)
        & (plot_data["log2FoldChange"] <= -log2fc_threshold),
        "category",
    ] = "Down in TNBC"

    colors = {
        "Up in TNBC": "#d73027",
        "Down in TNBC": "#4575b4",
        "Not significant": "#bdbdbd",
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    for category, frame in plot_data.groupby("category"):
        ax.scatter(
            frame["log2FoldChange"],
            frame["minus_log10_padj"],
            s=10,
            alpha=0.75,
            color=colors[category],
            label=category,
        )

    ax.axvline(log2fc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axvline(-log2fc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axhline(-math.log10(padj_threshold), color="black", linestyle="--", linewidth=1)
    ax.set_title("TNBC tumor vs Normal tissue")
    ax.set_xlabel("log2FoldChange = log2(TNBC / Normal)")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.legend()
    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def prepare_gene_symbol_counts(raw_counts: pd.DataFrame, gene_info: pd.DataFrame) -> pd.DataFrame:
    """把 Ensembl gene id 行名转换为 gene symbol，供生存分析使用。"""

    stripped_counts = strip_ensembl_versions(raw_counts)
    return convert_ensembl_to_hgnc(stripped_counts, gene_info)


def match_tnbc_expression_and_gene_list(
    log_expression: pd.DataFrame,
    tnbc_clinical: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
  # 匹配 TNBC 临床样本和表达矩阵，并从匹配后的矩阵生成基因列表。
    

    matched_expression = match_expression_to_clinical(log_expression, tnbc_clinical)
    matched_expression = matched_expression.apply(pd.to_numeric, errors="coerce")
    matched_expression = matched_expression.loc[matched_expression.notna().any(axis=1)].copy()
    gene_list = matched_expression.index.astype(str).tolist()
    if not gene_list:
        raise RuntimeError("TNBC 临床样本和表达矩阵匹配后没有可分析的基因。")
    return matched_expression, gene_list


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. 从 GDC 下载完整 TCGA-BRCA 临床病例，并整理出 ER/PR/HER2、OS、DFS 等字段。
    clinical = fetch_gdc_clinical(
        project=PROJECT,
        raw_output_file=OUTPUT_DIR / "gdc_clinical_cases_raw.json",
    )
    clinical.to_csv(OUTPUT_DIR / "clinical_gdc.csv", index=False)
    print(f"GDC 临床病例数: {len(clinical)}")

    # 2. 按 ER、PR、HER2 均为 Negative 的规则筛选 TNBC 患者。
    tnbc_clinical = filter_tnbc_samples(clinical)
    tnbc_clinical.to_pickle(OUTPUT_DIR / "TNBC_clin.pkl")
    print(f"三阴性临床样本数: {len(tnbc_clinical)}")

    # 3. 查询并下载 TCGA-BRCA STAR Counts 表达量文件；已存在文件会自动跳过。
    query = GDCquery(
        project=PROJECT,
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
    )
    manifest = getManifest(query)
    total_gb = int(manifest["size"].astype(int).sum()) / 1e9
    print(f"表达量文件数: {len(manifest)}")
    print(f"表达量文件总大小: {total_gb:.2f} GB")
    GDCdownload(query, files_per_chunk=FILES_PER_CHUNK, directory=GDC_DIRECTORY)

    # 4. 从本地 GDC 下载目录重新构建 query，确保列名使用真实 sample barcode。
    local_query = build_query_from_gdc_directory(GDC_DIRECTORY, project=PROJECT)
    print(f"本地发现表达量文件数: {len(local_query.results)}")

    # 5. 读取 STAR Counts 文件，得到 raw unstranded count matrix。
    prepared = GDCprepare(local_query, directory=GDC_DIRECTORY)
    raw_counts = assay(prepared)
    raw_counts.to_pickle(OUTPUT_DIR / "BRCA_raw_counts.pkl")
    print(f"BRCA raw counts 维度: {raw_counts.shape[0]} genes x {raw_counts.shape[1]} samples")

    # 6. 构建 DEG metadata。这里严格只保留 TNBC tumor 和 barcode code 为 11 的 Normal tissue。
    clinical_for_deg = standardize_clinical_columns(clinical)
    deg_metadata = build_deg_metadata(raw_counts, clinical_for_deg)
    deg_metadata.to_csv(OUTPUT_DIR / "deg_metadata_tnbc_vs_normal.csv", index=False)
    print("DEG 分组数量:")
    print(deg_metadata["condition"].value_counts())

    # 7. 构建基因 ID 到 gene symbol 的映射表。STAR Counts 自带 gene_name，可避免额外访问 BioMart。
    gene_info = gene_info_from_star_row_data(prepared.row_data)
    gene_info.to_csv(OUTPUT_DIR / "star_gene_name_map.csv", index=False)

    # 8. 执行 TNBC vs Normal 差异表达分析，并明确 log2FoldChange = log2(TNBC / Normal)。
    deg_config = build_deg_config()
    deg_counts = strip_ensembl_versions(raw_counts)
    deg_result = run_differential_expression(deg_counts, deg_metadata, gene_info, deg_config)
    deg_result.to_csv(OUTPUT_DIR / "TNBC_vs_Normal_deg_results.csv", index=False)
    print("DEG 结果前 10 行:")
    print(deg_result.head(10).loc[:, ["gene_id", "gene_symbol", "log2FoldChange", "pvalue", "padj"]])

    # 9. 绘制火山图。红色表示 TNBC 上调，蓝色表示 TNBC 下调。
    plot_volcano(deg_result, OUTPUT_DIR / "TNBC_vs_Normal_volcano.png")
    filter_upregulated_genes(deg_result).to_csv(OUTPUT_DIR / "TNBC_upregulated_genes.csv", index=False)
    filter_downregulated_genes(deg_result).to_csv(OUTPUT_DIR / "TNBC_downregulated_genes.csv", index=False)

    # 10. KM 生存分析只在 TNBC 肿瘤患者内部进行，不把 Normal tissue 放入生存分析。
    tnbc_tumor_sample_ids = deg_metadata.loc[deg_metadata["condition"].eq(CASE_GROUP), "sample_id"]
    tnbc_tumor_counts = raw_counts.loc[:, tnbc_tumor_sample_ids]
    tnbc_patient_counts = simplify_tcga_sample_ids(tnbc_tumor_counts)
    tnbc_patient_counts = drop_duplicate_sample_columns(tnbc_patient_counts)
    tnbc_gene_symbol_counts = prepare_gene_symbol_counts(tnbc_patient_counts, gene_info)
    tnbc_log_expression = log2_transform(tnbc_gene_symbol_counts)
    matched_tnbc_log_expression, survival_gene_list = match_tnbc_expression_and_gene_list(
        tnbc_log_expression,
        tnbc_clinical,
    )
    matched_tnbc_log_expression.to_pickle(OUTPUT_DIR / "TNBC_gene_exp_mat.pkl")
    pd.DataFrame({"gene": survival_gene_list}).to_csv(OUTPUT_DIR / "TNBC_matched_gene_list.csv", index=False)
    print(f"TNBC 临床和表达矩阵匹配后的样本数: {matched_tnbc_log_expression.shape[1]}")
    print(f"TNBC 临床和表达矩阵匹配后的基因数: {len(survival_gene_list)}")

    survival_data = build_survival_expression_data(
        tnbc_clinical,
        matched_tnbc_log_expression,
        survival_gene_list,
        status_col="DFS_STATUS",
        time_col="DFS_MONTHS",
    )
    survival_data.to_csv(OUTPUT_DIR / "merged_survival_expression.csv", index=False)
    print(f"可用于 TNBC DFS 生存分析的样本数: {len(survival_data)}")

    survival_results = analyze_gene_survival(
        survival_data,
        survival_gene_list,
        output_dir=OUTPUT_DIR,
        status_col="DFS_STATUS",
        time_col="DFS_MONTHS",
        event_value="1:Recurred/Progressed",
        plot=True,
        risk_table=True,
    )
    survival_results.to_csv(OUTPUT_DIR / "significant_genes.csv", index=False)
    print(survival_results.loc[:, ["gene", "status", "p_value"]])


if __name__ == "__main__":
    main()
