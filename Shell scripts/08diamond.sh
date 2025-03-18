#!/bin/sh
#SBATCH -o jobdia1060.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J diamond_10_60
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=80GB
#SBATCH --mail-type=end

source activate diamond
i=8
orf="/home/lipengwei2024phd/03GW/03mag_orfs/merged_mags_${i}_${i}"
out=/home/lipengwei2024phd/03GW/04mag_orfs_arg

#各种类型功能基因的数据库
ARGdb=/home/lipengwei2024phd/03GW/04mag_orfs_arg/database/SARG/SARG_v3.2_diamond.dmnd
MRGdb=/home/lipengwei2024phd/03GW/04mag_orfs_arg/database/MRG/MRG.dmnd
MGEdb=/home/lipengwei2024phd/03GW/04mag_orfs_arg/database/MGE/MGE.dmnd
AromaDegdb=/home/lipengwei2024phd/03GW/04mag_orfs_arg/database/AromaDeg/AromaDeg.fasta.dmnd
VFGdb=/home/lipengwei2024phd/03GW/04mag_orfs_arg/database/VFG_new/VFDB.dmnd


folder_name=$(basename ${orf})
echo "runing ${folder_name}"
mkdir -p ${out}/${folder_name}
query=${orf}/protein_seq.fa
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${ARGdb} -q ${query} -o ${out}/${folder_name}/${folder_name}_ARG_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${MGEdb} -q ${query} -o ${out}/${folder_name}/${folder_name}_MGE_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${MRGdb} -q ${query} -o ${out}/${folder_name}/${folder_name}_MRG_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${AromaDegdb} -q ${query} -o ${out}/${folder_name}/${folder_name}_AromaDeg_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${VFGdb} -q ${query} -o ${out}/${folder_name}/${folder_name}_VFG_out --id 60 --subject-cover 90
#--id 80 比较好
#The host system is detected to have 404 GB of RAM. It is recommended to increase the block size for better performance using these parameters : -b8 -c1
#--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#

#--id 80, 
echo "successfully"

