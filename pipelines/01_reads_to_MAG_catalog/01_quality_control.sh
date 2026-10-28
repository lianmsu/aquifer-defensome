#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J QC
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
SHARED_HOME="${SHARED_HOME:-/path/to/shared/home}"  # shared faculty storage (raw reads)
# ---------------------------------------------------------------------------

module load metaWRAP/1.2.2

pt=${SHARED_HOME}/83LPW/02RDN/01QC

for i in A804
do
echo Running Sample ${i}
unpigz  -k ${pt}/${i}/*gz

mv ${pt}/${i}/*1.fq   ${i}.raw_1.fq
mv ${pt}/${i}/*2.fq   ${i}.raw_2.fq

metawrap read_qc -1 ${i}.raw_1.fq -2 ${i}.raw_2.fq  -o ${i} --skip-bmtagger -t $SLURM_NTASKS

rm ${i}.raw_*.fq

done

