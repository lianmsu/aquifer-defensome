#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J prodigal_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7


module load prodigal/2.6.3-gcc-4.8.5
#A814 B814 C814 L1-814 L2-814 L3-814 L4-814 L5-814
#for i in A814 B814 C814 L1-814 L2-814 L3-814 L4-814 L5-814
for i in A403 A725 B403 B725 C403 C725
do
echo Running Sample ${i}

mkdir ${i}
prodigal -p meta -i ../02Assembly/${i}/contig_500bp_${i}.fa  \
				 -f gff -o ${i}/genes.gff -m -d ${i}/nucleotide_seq.fa -a ${i}/protein_seq.fa

done
