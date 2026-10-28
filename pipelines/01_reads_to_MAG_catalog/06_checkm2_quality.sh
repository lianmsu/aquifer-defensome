#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J checkm_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --mail-type=end

source  activate checkm2

mkdir checkm2
checkm2 predict --threads $SLURM_NTASKS --input ./dereplicated_genomes --output-directory ./checkm2 -x fa
