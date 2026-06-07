# TCGA-BRCA TNBC 分析脚本代码设计说明

本文档解释 `experiments/tnbc_survival_workflow.py` 中每一部分代码为什么这样设计。目标读者是新手 Python 开发者。

## 1. 为什么做成无参数脚本

这份脚本服务于一个固定分析任务：

```text
TCGA-BRCA 的 TNBC tumor samples vs Normal tissue samples
```

所以脚本不接收命令行参数，不提供跳过步骤，也不让用户手动选择输入输出路径。

这样做的好处：

- 初学者只需要运行一条命令。
- 每次运行的分析规则一致。
- 不容易因为参数传错导致分组错误。

## 2. 顶部常量区

脚本顶部定义了固定配置：

```python
PROJECT = "TCGA-BRCA"
GDC_DIRECTORY = Path("GDCdata")
OUTPUT_DIR = Path("TNBC_output")
FILES_PER_CHUNK = 10
CASE_GROUP = "TNBC"
CONTROL_GROUP = "Normal"
```

这些常量控制项目名、下载目录、输出目录和 DEG 分组名称。

新手可以先把它们理解为“脚本的默认设置”。

脚本没有写死固定基因列表。KM 生存分析使用的基因列表，会在 TNBC 临床样本和表达矩阵匹配之后动态生成。

## 3. barcode 解析函数

函数：

```python
parse_tcga_sample_type(barcode)
```

作用：

从 TCGA barcode 中识别样本类型。

核心规则：

```text
01 -> tumor
11 -> normal
其他 -> other
```

为什么需要它：

差异表达分析中的 Normal tissue 不能靠临床 receptor 状态判断，而要靠 TCGA barcode sample type code 判断。

## 4. TNBC 患者筛选函数

函数：

```python
select_tnbc_patients(clinical)
```

作用：

从临床表中筛选 ER、PR、HER2 全部为 Negative 的患者。

核心规则：

```text
er_status == Negative
pr_status == Negative
her2_status == Negative
```

为什么单独写函数：

TNBC 识别是整个分析的核心规则，单独写函数后可以很容易做单元测试，避免误把 ER 阳性或 HER2 Unknown 患者纳入 TNBC。

## 5. DEG metadata 构建函数

函数：

```python
build_deg_metadata(expression, clinical)
```

作用：

构建差异表达分析需要的样本分组表。

输出字段：

```text
sample_id
patient_id
sample_type
condition
```

设计规则：

- 如果样本是 tumor，并且患者是 TNBC，则 `condition = TNBC`。
- 如果样本是 normal，则 `condition = Normal`。
- 其他样本全部排除。

为什么这样设计：

用户明确要求 DEG contrast 是：

```text
TNBC tumor samples vs Normal tissue samples
```

而不是：

- 所有 BRCA tumor vs Normal
- TNBC vs 非 TNBC tumor
- TNBC high-risk vs low-risk

所以 metadata 中只允许出现：

```text
TNBC
Normal
```

## 6. DEG 配置对象

类：

```python
DegConfig
```

作用：

显式记录 case/control 和 contrast 名称。

默认：

```text
case_group = TNBC
control_group = Normal
contrast_name = TNBC_vs_Normal
```

为什么需要它：

差异表达最容易出错的地方就是比较方向。如果方向写反，`log2FoldChange` 的解释也会完全相反。

## 7. counts 归一化函数

函数：

```python
normalize_counts_median_ratio(counts)
```

作用：

对 raw counts 做样本层面的 size factor 校正。

设计思路：

- 用每个基因的几何均值作为参考。
- 每个样本计算 count / 几何均值的比例。
- 每个样本取比例的中位数作为 size factor。
- raw counts 除以 size factor 得到 normalized counts。

为什么这样做：

不同样本测序深度可能不同，直接比较 raw counts 会有偏差。这个方法借鉴了 DESeq2 的 median-of-ratios 思路。

注意：

脚本没有完整复刻 DESeq2 的负二项模型，只实现了轻量级、可直接运行的近似分析。

## 8. 差异表达函数

函数：

```python
run_differential_expression(counts, metadata, gene_info, config)
```

作用：

计算 TNBC vs Normal 的差异表达结果。

主要步骤：

1. 检查 metadata 中是否只有 TNBC 和 Normal。
2. 只保留 metadata 中出现的样本。
3. 删除全 0 基因。
4. 对 raw counts 做 size factor 归一化。
5. 计算 `log2(normalized_count + 1)`。
6. 计算 TNBC 组均值和 Normal 组均值。
7. 计算：

```text
log2FoldChange = mean(log2 TNBC) - mean(log2 Normal)
```

8. 计算近似 Wald statistic。
9. 计算 p value。
10. 使用 Benjamini-Hochberg 方法计算 `padj`。

输出字段：

```text
gene_id
gene_symbol
baseMean
log2FoldChange
lfcSE
stat
pvalue
padj
contrast
```

## 9. log2FoldChange 方向

脚本固定：

```text
log2FoldChange = log2(TNBC / Normal)
```

这意味着：

- 正值：TNBC 中更高。
- 负值：Normal 中更高，或 TNBC 中更低。

相关测试：

```python
test_log2fc_direction_tnbc_over_normal
test_run_differential_expression_outputs_tnbc_over_normal_direction
```

## 10. 火山图函数

函数：

```python
plot_volcano(deg_result, output_file)
```

作用：

把 DEG 结果画成火山图。

颜色规则：

- 红色：TNBC 上调
- 蓝色：TNBC 下调
- 灰色：不显著

默认显著标准：

```text
padj < 0.05
abs(log2FoldChange) >= 1
```

## 11. KM 生存分析设计

KM 分析和 DEG 分析不同。

DEG 分析比较：

```text
TNBC tumor vs Normal tissue
```

KM 分析只在：

```text
TNBC 患者内部
```

进行。

原因：

生存分析需要患者级临床随访数据，Normal tissue 不是肿瘤患者预后分组，不应该作为 KM 分析组。

KM 分析步骤：

1. 只取 DEG metadata 中 `condition == TNBC` 的 tumor 样本。
2. barcode 截取前 12 位变成 patientId。
3. 去重，保证一个患者只出现一次。
4. 转换 gene symbol。
5. 取 `log2(x + 1)`。
6. 用 TNBC 临床表中的 `patientId` 匹配表达矩阵列名。
7. 从匹配后的表达矩阵行名生成 KM 分析基因列表。
8. 合并 `DFS_STATUS` 和 `DFS_MONTHS`。
9. 每个匹配基因按中位数分为高低表达组。
10. 画 KM 曲线并计算 log-rank p 值。

## 12. 单元测试设计

测试文件：

```bash
tests/test_tnbc_deg_pipeline.py
```

测试覆盖：

- metadata 只能包含 TNBC 和 Normal。
- 非 TNBC 肿瘤不能进入 DEG。
- barcode `01` 识别为 tumor。
- barcode `11` 识别为 normal。
- ER/PR/HER2 全阴性才是 TNBC。
- ER Positive 患者被排除。
- HER2 Unknown 患者被排除。
- contrast 固定为 `TNBC_vs_Normal`。
- `log2FoldChange` 方向固定为 `log2(TNBC / Normal)`。
- KM 基因列表来自 TNBC 临床样本和表达矩阵匹配后的结果，不来自预设列表。

这些测试的目的，是防止后续维护时把核心分组逻辑改错。
