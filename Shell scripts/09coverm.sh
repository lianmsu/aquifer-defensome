#!/bin/sh
#SBATCH -o job.%j.%N.build_index.out
#SBATCH --partition=fat
#SBATCH -J index_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --mail-type=end

start_time=$(date +%s)

module load bowtie2
bins=./dereplicated_genomes

cat ${bins}/*.fa > ./merge.fa

####Step1:对contig建立index索引文件
echo "index build start"

mkdir index
bowtie2-build merge.fa index/all

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Total index build execution time is $(($cost_time/60))min " 

#!/bin/sh
#SBATCH -o job.%j.%N.coverm.out
#SBATCH --partition=fat
#SBATCH -J coverm_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --mail-type=end

start_time=$(date +%s)
module load samtools
module load bowtie2


reads=./01QC
bins=./dereplicated_genomes

mkdir sorted_bam #存放排序后的bam
#
for i in ${samples}
do
echo "${i} mapping"
start_time=$(date +%s)

mkdir ${i}
#bowtie2进行mapping
bowtie2 -x index/all -1 ${reads}/${i}/final_pure_reads_1.fastq.gz -2 ${reads}/${i}/final_pure_reads_2.fastq.gz -S ${i}/${i}.sam --threads $SLURM_NTASKS
#sam转化为bam并排序
samtools view -@ $SLURM_NTASKS -b -S ${i}/${i}.sam > ${i}/${i}.bam
samtools sort -@ $SLURM_NTASKS ${i}/${i}.bam > sorted_bam/${i}.sorted.bam

rm -rf ${i}
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "${i} $SLURM_JOB_NODELIST   $SLURM_NTASKS    execution time is $(($cost_time/60))min " >> Running_time_${SLURM_JOBID}.log
done

#计算rpkm
echo "abundance calculate start"
coverm_start_time=$(date +%s)

mkdir genomes_fna
cd genomes_fna
ln -s ${bins}/* ./
rename fa fna *fa
cd ..

source activate coverm
coverm genome \
-m rpkm \
--bam-files sorted_bam/*.bam \
--genome-fasta-directory genomes_fna  \
--min-read-percent-identity 95 \
--min-read-aligned-percent 90 \
--output-file bins_abundance_tpm.tsv \
--contig-end-exclusion 0 \
--no-zeros \
-t $SLURM_NTASKS

coverm_end_time=$(date +%s)
coverm_cost_time=$[ $coverm_end_time-$coverm_start_time ]
echo "abundance calculate time is $(($coverm_cost_time/60))min " >> Running_time_${SLURM_JOBID}.log






