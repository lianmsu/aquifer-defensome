# File mapping

Every script of this repository, its location in the previous (flat) layout and
its purpose. The scripts were **moved, renamed and made path-portable**; the
statistical and bioinformatic logic is unchanged except for the typo listed
under [Corrections](#corrections).


Code edits applied to every file:

1. the site-specific absolute paths were replaced by variables -- `GTR_DATA` and
   `GTR_REV` in the R scripts, `LUSTRE_HOME` / `HOME_ROOT` / `DB_ROOT` / `SHARED_HOME`
   in the pipeline scripts (see [below](#site-specific-paths-and-how-they-are-hidden));
2. a header was added to every R script (purpose, manuscript figure, outputs, original
   path);
3. `#SBATCH --reservation=<account>...` directives were commented out;
4. line endings were normalised to LF (`.gitattributes`) and the R/python/shell
   artefacts were added to `.gitignore`.

## Upstream pipelines

| New path | Previous path | Purpose |
| --- | --- | --- |
| `pipelines/01_reads_to_MAG_catalog/01_quality_control.sh` | `Shell scripts/01Quality_control.sh` | Read QC with metaWRAP (unpigz decompression + `metawrap read_qc`, host decontamination skipped). |
| `pipelines/01_reads_to_MAG_catalog/02_assembly.sh` | `Shell scripts/02Assembly_mh.sh` | Per-sample metagenomic assembly with metaWRAP/megahit; keeps a 1 kb and a 0.5 kb contig set. |
| `pipelines/01_reads_to_MAG_catalog/03_binning.sh` | `Shell scripts/03binning.sh` | metaWRAP binning with metabat2, maxbin2 and concoct, followed by bin refinement. |
| `pipelines/01_reads_to_MAG_catalog/04_dereplicate_bins.sh` | `Shell scripts/04drep.sh` | dRep dereplication of the refined bins from all samples (-sa 0.98 -nc 0.30, >=50% completeness, <=5% contamination). |
| `pipelines/01_reads_to_MAG_catalog/05_gtdbtk_taxonomy.sh` | `Shell scripts/05GTDB.sh` | GTDB-Tk (r207) taxonomic classification of the dereplicated MAGs. |
| `pipelines/01_reads_to_MAG_catalog/06_checkm2_quality.sh` | `Shell scripts/06checkm.sh` | CheckM2 completeness/contamination estimates of the dereplicated MAGs. |
| `pipelines/02_gene_prediction/07_prodigal_orf.sh` | `Shell scripts/07Prodigal-ORF.sh` | Prodigal gene prediction per assembled sample (meta mode), producing the ORF protein sets. |
| `pipelines/02_gene_prediction/22_prokka_annotate.sh` | `Shell scripts/22prokka.sh` | Prokka annotation of the 1,626 high-quality MAGs (needed by DefensePredictor). |
| `pipelines/03_abundance_and_function/08_diamond_screen.sh` | `Shell scripts/08diamond.sh` | DIAMOND blastp of the ORFs against the ARG / MRG / mibPOP / VFG databases, run in MAG batches. |
| `pipelines/03_abundance_and_function/09_coverm_abundance.sh` | `Shell scripts/09coverm.sh` | bowtie2 index of the dereplicated MAGs and coverm read recruitment, giving per-sample MAG abundance (RPKM/TPM). |
| `pipelines/04_mobile_genetic_elements/10_genomad.sh` | `Shell scripts/10genomad.sh` | geNomad detection of plasmids and viruses/proviruses in the MAG contigs. |
| `pipelines/04_mobile_genetic_elements/11_integron_finder.sh` | `Shell scripts/11integronfinder.sh` | IntegronFinder detection of integrons in the 2,000-MAG batches. |
| `pipelines/04_mobile_genetic_elements/12_icefinder.sh` | `Shell scripts/12ICEfinder.sh` | ICEfinder2 detection of integrative conjugative elements, 2,000 MAGs per job. |
| `pipelines/04_mobile_genetic_elements/30_minced_crispr.sh` | `Shell scripts/30minced.sh` | minced detection of CRISPR arrays and spacers on the Cas-containing contigs. |
| `pipelines/04_mobile_genetic_elements/31_blast_spacers.sh` | `Shell scripts/31blast.sh` | Spacer-to-prophage BLAST search: builds the prophage database and maps CRISPR spacers onto it (historical phage-host links). |
| `pipelines/05_defense_systems/20_defensefinder.sh` | `Shell scripts/20defensefinder.sh` | DefenseFinder 2.0.0 run on the MAG protein sets, 1,000 MAGs per job (`--preserve-raw`, needed for defense islands). |
| `pipelines/05_defense_systems/21_defense_predictor.sh` | `Shell scripts/21defensepredictor.sh` | DefensePredictor on the Prokka GFF/FFN/FAA of the MAGs; used to cross-check DefenseFinder. |
| `pipelines/06_phylogeny_and_traits/40_acr_phylogeny.sh` | `Shell scripts/40phylogenetic.sh` | Collect the anti-CRISPR proteins of prophages, add the PaCRISPR reference set, then muscle + trimAl (-gt 0.5) + FastTree. |
| `pipelines/06_phylogeny_and_traits/50_ogt_prediction.sh` | `Shell scripts/50OGT.sh` | OGT_prediction pipeline: optimal growth temperature for each MAG, using tax_info.txt as the taxonomy input. |
| `pipelines/07_table_postprocessing/01_collect_DI_hits.py` | `Shell scripts/others/01defenseislands.py` | Collect the per-batch DefenseFinder island TSV files, keep Defense rows and standardise the contig IDs. |
| `pipelines/07_table_postprocessing/02_DI_summary.py` | `Shell scripts/others/02DI_res.py` | Parse `defense_islands.txt` into a tidy per-island summary (gene counts, defense types, boundaries). |
| `pipelines/07_table_postprocessing/03_DI_nondG_orf_ids.py` | `Shell scripts/others/03DI_nonDG.sh` (lines 1-33 of that file) | Expand every defense-island record into the ORF IDs of its contig (Python block extracted from the original mixed scratch file). |
| `pipelines/07_table_postprocessing/04_DI_nondG_filter.sh` | `Shell scripts/others/03DI_nonDG.sh` (lines 35-71 of that file) | Subtract the defense-gene IDs from the DI ORF list and extract the remaining ORFs with seqkit (shell block of the same scratch file). |
| `pipelines/07_table_postprocessing/10_integron_DG_extract.sh` | `Shell scripts/others/10integron_DGextract.sh` | Extract the defense systems carried by the integrons from the per-batch DefenseFinder output. |
| `pipelines/07_table_postprocessing/11_DS_filter.py` | `Shell scripts/others/11DS_filter.py` | Clean the integron defense-gene table: extract the genome name from each `hit_id`. |
| `pipelines/07_table_postprocessing/12_DS_bind.py` | `Shell scripts/others/12DS_abind.py` | Bind the cleaned integron tables and count defense genes per genome and type. |
| `pipelines/07_table_postprocessing/20_DP_filter.sh` | `Shell scripts/others/20DP_filter.sh` | Merge the per-MAG DefensePredictor output, keep predictions with log-odds >= 7.2 and drop the sequence columns. |

## R analyses

| New path | Previous path | Manuscript | Title |
| --- | --- | --- | --- |
| `analyses/01_main_figures/Fig01.R` | `R codes/Fig1.R` | Fig. 1 | Defensome overview |
| `analyses/01_main_figures/Fig02.R` | `R codes/Fig2.R` | Fig. 2 | Defensome of CPR, Pseudomonadota and other phyla |
| `analyses/01_main_figures/Fig03.R` | `R codes/Fig3.R` | Fig. 3 | VHR/DHR and the defensome |
| `analyses/01_main_figures/Fig04.R` | `R codes/Fig4.R` | Fig. 4 | Antibiotic resistance genes and the defensome |
| `analyses/01_main_figures/Fig05_A-D.R` | `R codes/Fig5_A-D.R` | Fig. 5a-d | Mobile genetic elements and the defensome |
| `analyses/01_main_figures/Fig05_E-I.R` | `R codes/Fig5_E-I.R` | Fig. 5e-i | Defense islands |
| `analyses/01_main_figures/Fig06.R` | `R codes/Fig6.R` | Fig. 6 | Arms-race figure assembly |
| `analyses/02_supplementary_figures/FigS02_B.R` | `R codes/Fig.S2_B.R` | Fig. S2b | MAG abundance by niche breadth |
| `analyses/02_supplementary_figures/FigS02_C.R` | `R codes/Fig.S2_C.R` | Fig. S2c | DS density across geological settings |
| `analyses/02_supplementary_figures/FigS02_D.R` | `R codes/Fig.S2_D.R` | Fig. S2d | Defense families across zones |
| `analyses/02_supplementary_figures/FigS03.R` | `R codes/Fig.S3.R` | Fig. S3 | Defense genes vs genome size |
| `analyses/02_supplementary_figures/FigS04.R` | `R codes/Fig.S4.R` | Fig. S4 | Defensome of archaea |
| `analyses/02_supplementary_figures/FigS06.R` | `R codes/Fig.S6.R` | Fig. S6 | DefenseFinder vs DefensePredictor |
| `analyses/02_supplementary_figures/FigS07.R` | `R codes/Fig.S7.R` | Fig. S7 | VHR and genome size |
| `analyses/02_supplementary_figures/FigS08.R` | `R codes/Fig.S8.R` | Fig. S8 | DS density across geological settings (with source data) |
| `analyses/02_supplementary_figures/FigS09.R` | `R codes/Fig.S9.R` | Fig. S9 | Plaque assay experiment |
| `analyses/02_supplementary_figures/FigS10.R` | `R codes/Fig.S10.R` | Fig. S10 | DS distribution across phyla and genera (data-prep script) |
| `analyses/02_supplementary_figures/FigS12.R` | `R codes/Fig.S12.R` | Fig. S12 | DS distribution in selected phyla |
| `analyses/02_supplementary_figures/FigS13.R` | `R codes/Fig.S13.R` | Fig. S13 | Mobile resistance genes and the defensome |
| `analyses/02_supplementary_figures/FigS14.R` | `R codes/Fig.S14.R` | Fig. S14 | Virulence factor genes and the defensome |
| `analyses/02_supplementary_figures/FigS15.R` | `R codes/Fig.S15.R` | Fig. S15 | Prophages and plasmids across phyla |
| `analyses/02_supplementary_figures/FigS16.R` | `R codes/Fig.S16.R` | Fig. S16 | Function of non-defense genes in defense islands |
| `analyses/02_supplementary_figures/FigS17-18.R` | `R codes/Fig.S17-18.R` | Fig. S17, Fig. S18 | Gene-level maps of representative elements (data-prep script) |

## Corrections

| File | Change | Why |
| --- | --- | --- |
| `analyses/02_supplementary_figures/FigS17-18.R` | `c(first_index ：last_index)` -> `c(first_index:last_index)` | the original file contained a full-width colon (typed with a Chinese IME) and therefore did not parse in R; the surrounding code clearly builds the row index range of each gene block |

Apart from this one-character fix, no line of analysis code was altered: only the path
statements (and the headers) differ from the original files.

## Site-specific paths and how they are hidden

The absolute paths of the production runs -- user home directories, a shared faculty
storage account, local Windows analysis folders -- have been removed from every file.
They are replaced by variables, so nothing about the private set-up is published and the
same script can be re-pointed without touching the logic:

| Variable | Replaces | Used in |
| --- | --- | --- |
| `GTR_DATA` | the root of the data directory used for the study | every script in `analyses/` |
| `GTR_REV` | the working copy of the peer-review revision analyses | `FigS06.R`, `FigS09.R`, `FigS15.R`, `Fig05_A-D.R` |
| `LUSTRE_HOME` | home directory on the large-storage filesystem | `pipelines/*` |
| `HOME_ROOT` | home directory on the login-node filesystem | `pipelines/*` |
| `DB_ROOT` | reference databases (`${HOME_ROOT}/database` by default) | `pipelines/04_mobile_genetic_elements/10_genomad.sh`, `pipelines/07_table_postprocessing/20_DP_filter.sh` |
| `SHARED_HOME` | the shared faculty storage account that holds the raw reads | `pipelines/01_reads_to_MAG_catalog/01_quality_control.sh` |
| `GTR_DI_DIR` | directory of the per-batch DefenseFinder island tables | `pipelines/07_table_postprocessing/01_collect_DI_hits.py` |
| `GTR_INTEGRON_DIR` | directory of the integron defense-gene tables | `pipelines/07_table_postprocessing/11_DS_filter.py`, `12_DS_bind.py` |

Each pipeline script declares only the variables it uses, in a `# Site-specific settings`
block right below its `#SBATCH` header, with `/path/to/...` placeholders as defaults.
The R scripts read the variable at the top of the file:

```r
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
```

`normalizePath()` makes the path absolute at load time, so it stays valid after the
numerous `setwd()` calls inside the scripts.

