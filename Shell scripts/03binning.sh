#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J binning_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mail-type=end

qc=./01QC
zz=./02Assembly

module load metaWRAP/1.2.2

#for i in A814 B814 C814 L1-814 L2-814 L3-814 L4-814 L5-814
#for i in A403 A725 B403 B725 C403 C725
for i in C725
do

start_time=$(date +%s)

echo Running Sample ${i}
mkdir ${i}

unpigz -k ../01QC/${i}/*_1.fastq.gz  -p $SLURM_NTASKS
unpigz -k ../01QC/${i}/*_2.fastq.gz  -p $SLURM_NTASKS


metawrap binning -o ${i}/initial_binning -t $SLURM_NTASKS  -m 2000  \
                 -a ${zz}/${i}/contig_1kb_*.fa  \
                 --metabat2 --maxbin2 --concoct \
                 ${qc}/${i}/final_pure_reads_*.fastq --run-checkm


metawrap bin_refinement -o ${i}/refinement -c 50 -x 10  -t $SLURM_NTASKS -m 2000 \
                        -A ${i}/initial_binning/metabat2_bins/ \
                        -B ${i}/initial_binning/maxbin2_bins/ \
                        -C ${i}/initial_binning/concoct_bins/ 

metawrap reassemble_bins  -t $SLURM_NTASKS  -m 2000 -c 50 -x 5 \
                          -o ${i}/BIN_REASSEMBLY \
                          -1 ${qc}/${i}/final_pure_reads_1.fastq \
                          -2 ${qc}/${i}/final_pure_reads_2.fastq \
                          -b ${i}/refinement/metawrap_50_10_bins

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "${i} mapping time is $(($cost_time/60))min $(($cost_time%60))s" >> run_time_binning_${SLURM_JOBID}.log

rm ../01QC/${i}/*_*.fastq

done
