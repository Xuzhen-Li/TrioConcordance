# TrioConcordance

[English](README.md)

Version: 1.0.0

Author: Xuzhen Li

## 概述

TrioConcordance 是一个 Python 工具，用于在三个独立的基因型数据集（如全基因组测序、目标区域捕获和参考面板）之间进行逐位点和逐样本的一致性分析。

该工具包含两个主要分析模块。10-Way State Classification 根据状态同源（IBS）逻辑将等位基因状态分为 10 个互斥类别。Error Transition Profiling 提取特定的不匹配模式（例如目标捕获错误，即 `WGS == Reference != Capture`）并分析转换频率。

## 依赖项

需要 Python 3.6+、pandas 和 matplotlib。可以使用 pip 安装所需包：

```bash
pip install pandas matplotlib
```

## 数据准备

软件需要三个制表符分隔的基因型矩阵。可以使用 `bcftools` 直接从 VCF 文件中提取这些矩阵。

注意：三个输入矩阵必须具有相同的维度和结构（即完全相同顺序的位点和样本）。

```bash
# Standard command to extract genotype matrix from VCF
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input_wgs.vcf > wgs_matrix.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input_cap.vcf > cap_matrix.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input_ref.vcf > ref_matrix.txt
```

## 用法

```bash
python TrioConcordance.py \
    -w wgs_matrix.txt \
    -c cap_matrix.txt \
    -r ref_matrix.txt \
    -o ./Results \
    -p MyAnalysis
```

### 参数

| 参数 | 描述 | 是否必需 | 默认值 |
| --- | --- | --- | --- |
| `-w`, `--wgs` | WGS 基因型矩阵的路径 | 是 | - |
| `-c`, `--cap` | 捕获数据基因型矩阵的路径 | 是 | - |
| `-r`, `--crr` | 参考基因型矩阵的路径 | 是 | - |
| `-o`, `--outdir` | 输出目录路径 | 否 | `Results` |
| `-p`, `--prefix` | 输出文件前缀 | 否 | `Trio` |

## 输出文件

流程将在指定目录下生成四个标准输出文件。

- `[prefix]_Global_Stats.tsv`: 所有样本和位点 10 分类的全局摘要。
- `[prefix]_Sample_Stats.tsv`: 按单个样本划分的详细 10 分类统计。
- `[prefix]_Capture_Error_Transitions.tsv`: 特定基因型错误转换的频率。
- `[prefix]_Concordance_Distribution.pdf`: 可视化每个样本在各一致性类别中位点比例的堆叠柱状图。
