#!/bin/bash
#SBATCH --job-name=p1626
#SBATCH --partition=cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=28

#SBATCH --output=job.lpw-prokka_bug4.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

source activate prokka
cd ${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs
file=${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/1626mags.id
for i in $(cat $file)
do
prokka ${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/MAGs/${i} \
--prefix ${i} \
--outdir ${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/prokka_1626/${i} \
--centre X --compliant --locustag ${i} \
--kingdom Bacteria --cpus 28;
done
