# Software and environment

## Cluster / job scheduler

All `pipelines/` scripts were submitted with **Slurm** (`#SBATCH` headers). The production
runs used two queues, `cpu` and `fat`. Job names, memory requests and thread counts are
kept in the scripts as a record of the actual runs.

The site-specific absolute paths are variables, so each script reads them from its
environment (or from the `# Site-specific settings` block just below the `#SBATCH`
header):

| Variable | Meaning | Default in the script |
| --- | --- | --- |
| `LUSTRE_HOME` | home directory on the large-storage filesystem | `/path/to/lustre/home` |
| `HOME_ROOT` | home directory on the login-node filesystem | `/path/to/home` |
| `DB_ROOT` | reference databases (geNomad, CARD, ...) | `${HOME_ROOT}/database` |
| `SHARED_HOME` | shared faculty storage account holding the raw reads | `/path/to/shared/home` |
| `GTR_DI_DIR` | directory with the per-batch DefenseFinder island TSVs | `/path/to/.../DIs` |
| `GTR_INTEGRON_DIR` | directory with the integron defense-gene tables | `.` (current directory) |

To re-run a script elsewhere:

```bash
export LUSTRE_HOME=/your/lustre/home
export HOME_ROOT=/your/home
export DB_ROOT=/your/databases
sbatch pipelines/01_reads_to_MAG_catalog/05_gtdbtk_taxonomy.sh
```

Also adapt:

* `#SBATCH --partition=...` (queue name),
* `source activate <env>` — replace with `conda activate <env>` on modern conda,
* `# #SBATCH --reservation=...` — commented out in the repository (the reservation name
  of the production account is site-specific); uncomment and set your own if needed.

Python is `python3` with `pandas` (and `csv`/`re` from the standard library); the
table post-processing scripts run without a scheduler.

## Conda environments and modules (as used)

| Environment / module | Used by | Software |
| --- | --- | --- |
| `module load metaWRAP/1.2.2` | `01_quality_control.sh`, `02_assembly.sh`, `03_binning.sh` | metaWRAP (read_qc, assembly, binning) |
| `drep` + `module load CheckM/1.1.2` | `04_dereplicate_bins.sh` | dRep |
| `gtdbtk2.3` | `05_gtdbtk_taxonomy.sh` | GTDB-Tk (r207 release) |
| `checkm2` | `06_checkm2_quality.sh` | CheckM2 |
| `module load prodigal/2.6.3-gcc-4.8.5` | `07_prodigal_orf.sh` | Prodigal 2.6.3 |
| `prokka` | `22_prokka_annotate.sh` | Prokka |
| `diamond` | `08_diamond_screen.sh` | DIAMOND |
| `coverm` + `module load bowtie2`, `samtools` | `09_coverm_abundance.sh` | bowtie2, samtools, coverm |
| `genomad` | `10_genomad.sh` | geNomad |
| `integron_finder` | `11_integron_finder.sh` | IntegronFinder |
| `icefinder` | `12_icefinder.sh` | ICEfinder2 |
| *(none)* | `30_minced_crispr.sh` | minced |
| `module load BLAST+/2.11.0-gompi-2021a` | `31_blast_spacers.sh` | BLAST+ 2.11 |
| `defensefinder2.0.0` | `20_defensefinder.sh` | DefenseFinder 2.0.0 |
| `defense_predictor` + `module load gcc/8.3.0` | `21_defense_predictor.sh` | DefensePredictor |
| *(seqkit, muscle, trimAl, FastTree)* | `40_acr_phylogeny.sh` | anti-CRISPR tree |
| `OGT_prediction` | `50_ogt_prediction.sh` | OGT_prediction |

Database paths (geNomad db, CARD, MRG, mibPOP, VFG, PaCRISPR) appear as variables or
inline paths in the scripts and must be downloaded separately.

## R

The analyses are plain R scripts (no `renv`, no project file). They were developed under
R 4.x; all 23 scripts parse cleanly under R 4.6.0. Packages used, by frequency:

* data handling: `dplyr`, `tidyr`, `tidyverse`, `stringr`, `reshape2`, `reshape`, `readxl`
* plotting: `ggplot2`, `cowplot`, `Cairo`, `ggpubr`, `ggsignif`, `ggsci`, `ggsankeyfier`,
  `gggenes`, `RColorBrewer`, `scales`, `treemap`, `treemapify`, `pheatmap`,
  `ComplexHeatmap`, `circlize`, `UpSetR`, `export`
* statistics / networks: `rstatix`, `igraph`

`Cairo` needs the system `libcairo`; `UpSetR` and `export` are CRAN packages.

A minimal install:

```r
install.packages(c("ggplot2", "dplyr", "cowplot", "ggpubr", "tidyr", "Cairo",
                   "rstatix", "stringr", "ggsignif", "pheatmap", "scales",
                   "reshape2", "RColorBrewer", "readxl", "UpSetR"))
install.packages(c("ComplexHeatmap", "circlize", "treemapify", "ggsankeyfier",
                   "gggenes", "igraph", "treemap", "export"))
```

## Reproducibility notes

* Scripts are self-contained: no shared helper file, no `source()` between analyses
  (except the two data-prep scripts described in
  [01_pipeline_overview.md](01_pipeline_overview.md)).
* Figures are written as PDF with explicit `width`/`height` in inches, matching the
  submission requirements; several scripts also export the underlying numbers as CSV
  (source-data files).
* `GTR_DATA` is resolved with `normalizePath()` at load time, so the dozens of `setwd()`
  calls inside the scripts cannot invalidate it.
