# Directory layout

The repository has two code layers plus documentation. Both layers use the same
conventions: a two-digit directory prefix gives the execution / narrative order, and a
file prefix gives the order inside a directory.

```
aquifer-defensome/
├── README.md                  English overview
├── README_zh.md               中文说明
├── docs/                      workflow, layout, mapping, data, software
├── pipelines/                 layer 1: bioinformatics, Slurm + conda
└── analyses/                  layer 2: statistics and figures, R
```

## `pipelines/` — upstream, ordered by stage

| Directory | Scope |
| --- | --- |
| `01_reads_to_MAG_catalog/` | QC, assembly, binning, dRep dereplication, GTDB-Tk taxonomy, CheckM2 quality |
| `02_gene_prediction/` | Prodigal ORF prediction, Prokka annotation of the high-quality MAGs |
| `03_abundance_and_function/` | DIAMOND screening against the ARG/MRG/mibPOP/VFG databases, coverm abundance |
| `04_mobile_genetic_elements/` | geNomad, IntegronFinder, ICEfinder2, minced CRISPR, spacer BLAST |
| `05_defense_systems/` | DefenseFinder 2.0.0, DefensePredictor |
| `06_phylogeny_and_traits/` | anti-CRISPR phylogeny (muscle/trimAl/FastTree), OGT prediction |
| `07_table_postprocessing/` | defense-island and integron tables (python scripts + awk/seqkit shell steps) |

The scripts are the Slurm job scripts as they were submitted, with their `#SBATCH`
headers preserved, so that the computational environment is documented. The site-specific
paths are variables — see the `# Site-specific settings` block below each `#SBATCH`
header, and [05_software_and_environment.md](05_software_and_environment.md).

File names keep their original pipeline number (`01_` … `07_`, then `08_`, `09_`, `10_`,
`11_`, `12_`, `20_`, `21_`, `22_`, `30_`, `31_`, `40_`, `50_`), so the order inside a
stage directory still reflects the global pipeline order.

## `analyses/` — downstream, ordered by the manuscript

| Directory | Content |
| --- | --- |
| `01_main_figures/` | `Fig01.R` – `Fig06.R` (Fig. 5 is split into `Fig05_A-D.R` and `Fig05_E-I.R`) |
| `02_supplementary_figures/` | `FigS02_B.R` – `FigS18.R` |

The file name is the figure name: `Fig05_E-I.R` is Fig. 5 panels e–i, `FigS17-18.R`
produces both Fig. S17 and Fig. S18. Two scripts are not only figure scripts but also
data-prep steps (`FigS10.R`, `FigS17-18.R`); the dependency table is at the end of
[01_pipeline_overview.md](01_pipeline_overview.md).

## Script header

Every R script starts with the same block, which is also the fastest way to navigate the
repository:

```r
# aquifer-defensome -- main figures / Defensome overview
# Prevalence of defense systems in the 27,578 MAGs, ...
#
# Manuscript      : Fig. 1
# Main outputs    : 1_DSs_distribution.pdf, DSs_percentage.pdf, ...
# Original script : R codes/Fig1.R
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
```

`GTR_DATA` is the root of the data directory; `GTR_REV` (declared only by the four scripts
that need it) is the working copy of the peer-review revision analyses.

## `docs/`

| File | Content |
| --- | --- |
| `01_pipeline_overview.md` | end-to-end workflow, stage by stage, with dependencies |
| `02_directory_layout.md` | this file |
| `03_file_mapping.md` | old path -> new path for every script, plus the one code correction |
| `04_data_availability.md` | working directory and input tables of every R script |
| `05_software_and_environment.md` | conda environments, R packages, cluster notes |
