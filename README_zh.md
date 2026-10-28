# aquifer-defensome

> **The defensome of prokaryotes in aquifers** 论文的代码仓库
> （Li, P. *et al.*, *Nature Communications* **16**, 6482 (2025)）

本仓库收录了从原始宏基因组测序数据到论文全部图件的所有代码：上游生物信息学流程（从质控、
组装、分箱到 27,578 个去冗余 MAG，以及移动遗传元件、防御/抗防御系统注释），以及下游 R
统计分析与绘图脚本。

## 目录结构

```
aquifer-defensome/
├── pipelines/          # 上游流程：Slurm 作业脚本 + conda 环境（shell / python）
│   ├── 01_reads_to_MAG_catalog/       质控、组装、分箱、dRep 去冗余、GTDB-Tk、CheckM2
│   ├── 02_gene_prediction/            Prodigal 预测 ORF、Prokka 注释
│   ├── 03_abundance_and_function/     DIAMOND 功能筛选、coverm 丰度
│   ├── 04_mobile_genetic_elements/    geNomad、IntegronFinder、ICEfinder2、minced、BLAST
│   ├── 05_defense_systems/            DefenseFinder 2.0.0、DefensePredictor
│   ├── 06_phylogeny_and_traits/       抗 CRISPR 系统发育树、OGT 预测
│   └── 07_table_postprocessing/       防御岛与 integron 中间表（python / awk）
├── analyses/           # 下游分析：R 统计与绘图
│   ├── 01_main_figures/               Fig01 – Fig06
│   └── 02_supplementary_figures/      FigS02_B – FigS18
└── docs/               # 流程、目录、文件对照表、数据与软件说明
```

每个 R 脚本以它对应的图件命名，脚本头部同样写明所属图件；脚本与新目录的完整对照见
[docs/03_file_mapping.md](docs/03_file_mapping.md)。

## 如何复现

**1. 上游流程**：`pipelines/` 下是当时实际投递的 SLURM 作业脚本（保留 `#SBATCH` 头）。
与运行环境相关的绝对路径已改为变量：导出 `LUSTRE_HOME`、`HOME_ROOT`、`DB_ROOT`（必要时
还有 `SHARED_HOME`），改好队列名，按 [docs/01](docs/01_pipeline_overview.md) 的顺序提交。

**2. R 分析**：需要 R ≥ 4.2 及 [docs/05](docs/05_software_and_environment.md) 中列出的 R 包。
输入数据通过 `GTR_DATA` 定位，默认取当前工作目录：

```bash
cd /path/to/data && Rscript /path/to/aquifer-defensome/analyses/01_main_figures/Fig01.R
```

另有 4 个脚本需要 `GTR_REV`（审稿修订分析的工作目录）。

### 两个脚本必须先运行

其余脚本相互独立，但下面两个脚本同时承担"数据准备"角色，它们创建的对象被别的脚本使用：

| 先运行 | 之后才能运行 | 因为它定义了 |
| --- | --- | --- |
| `analyses/02_supplementary_figures/FigS10.R` | `Fig02.R`、`Fig03.R`、`Fig05_A-D.R`、`FigS02_B.R`、`FigS07.R` | `ds_density_tax`（以及为 Fig02 提供 `md_p_val1`） |
| `analyses/02_supplementary_figures/FigS17-18.R` | `Fig06.R` | `prophage_anti_coloca_4` 及基因图数据 |

详见 [docs/01](docs/01_pipeline_overview.md) 末尾的依赖表。

## 数据

仓库只放代码。MAG 目录与注释结果随论文以 Supplementary Data 提供；原始测序数据在 NCBI SRA
的登录号见论文。脚本所需的数据目录结构见
[docs/04_data_availability.md](docs/04_data_availability.md)。

## 关于旧版仓库

旧版本是两个平铺的目录（`R codes/`、`Shell scripts/`）。本次整理为上面的分层结构：脚本
一个都没有丢（23 个 R 脚本、27 个流程脚本全部保留），路径改为可移植，并修复了 1 处导致
`Fig.S17-18.R` 无法解析的全角冒号笔误。仓库中已不含任何绝对路径、个人主目录或账号名：
全部由 `GTR_DATA`、`GTR_REV`、`LUSTRE_HOME`、`HOME_ROOT`、`DB_ROOT`、`SHARED_HOME`、
`GTR_DI_DIR`、`GTR_INTEGRON_DIR` 变量控制。

## 引用

Li, P. *et al.* The defensome of prokaryotes in aquifers. *Nat. Commun.* **16**, 6482
(2025).
