#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J assembly_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
module load metaWRAP/1.2.2

#for i in A814 B814 C814 L1-814 L2-814 L3-814 L4-814 L5-814
for i in A403 A725 B403 B725 C403 C725
do
echo Running Sample ${i}
start_time=$(date +%s)

metawrap assembly -1 ../01QC/${i}/*_1.fastq.gz -2 ../01QC/${i}/*_2.fastq.gz  -o ${i}  -m 2100 -t $SLURM_NTASKS

mv ${i}/final_assembly.fasta ${i}/contig_1kb_${i}.fa

awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${i}/megahit/final.contigs.fa | awk '{if(length($5)>=500) print}'| perl -p -e 's/\t/\n/g' | fold -w 60 > ./${i}/contig_500bp_${i}.fa

rm ${i}/megahit  -rf

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"


done
