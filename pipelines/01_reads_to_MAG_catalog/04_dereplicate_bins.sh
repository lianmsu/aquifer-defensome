#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J bin_drep_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mail-type=end

source  activate drep 
module load CheckM/1.1.2

pt=./10binning

#合并到一起
mkdir tmp
for i in ${samples}
do
ln -s  ${pt}/${i}/BIN_REASSEMBLY/reassembled_bins/*.fa tmp
done

#去冗余
dRep dereplicate res_merged -g tmp/*.fa -sa 0.98 -nc 0.30 -p $SLURM_NTASKS  -comp 50 -con 5

