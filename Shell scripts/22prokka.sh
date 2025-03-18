#!/bin/bash
#SBATCH --job-name=p1626
#SBATCH --partition=cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=28

#SBATCH --output=job.lpw-prokka_bug4.%j.%N.out

source activate prokka
cd /lustre/home/lipengwei2024phd/prokaryote_revision/HQ_MAGs
file=/lustre/home/lipengwei2024phd/prokaryote_revision/HQ_MAGs/1626mags.id
for i in $(cat $file)
do
prokka /lustre/home/lipengwei2024phd/prokaryote_revision/HQ_MAGs/MAGs/${i} \
--prefix ${i} \
--outdir /lustre/home/lipengwei2024phd/prokaryote_revision/HQ_MAGs/prokka_1626/${i} \
--centre X --compliant --locustag ${i} \
--kingdom Bacteria --cpus 28;
done
