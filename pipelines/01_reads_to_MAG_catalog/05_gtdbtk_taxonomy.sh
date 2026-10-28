#!/bin/bash
#SBATCH -o job.%j.%N.out 
#SBATCH --partition=fat
#SBATCH -J GTDB_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --mail-type=end
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
SHARED_HOME="${SHARED_HOME:-/path/to/shared/home}"  # shared faculty storage (raw reads)
# ---------------------------------------------------------------------------

source activate gtdbtk2.3

binning_path=./dereplicated_genomes
#out_path=${SHARED_HOME}/11basic_lake/02DB/07GTDB



start_time=$(date +%s)
 
mkdir res2
gtdbtk classify_wf --cpus $SLURM_NTASKS --pplacer_cpus  $SLURM_NTASKS  -x fa \
    --genome_dir ${binning_path}/  \
    --mash_db ./mash  \
    --out_dir res2

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



#mkdir ${out_path}/res

#gtdbtk classify_wf --cpus 14 \
#--genome_dir ${binning_path}/res_merged/dereplicated_genomes \
#--out_dir ${out_path}/res \
#-x fa \
#--pplacer_cpus 1

