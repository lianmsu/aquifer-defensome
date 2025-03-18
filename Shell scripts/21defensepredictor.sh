#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J DP_try
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100GB

# 开始计时
start_time=$(date +%s)

module load gcc/8.3.0
source activate defense_predictor
#########################################################################

#########################################################################

prodigal_gff=/home/lipengwei2024phd/03GW/try/B6_bin.5_prokka/B6_bin.5.fasta.gff
prodigal_ffn=/home/lipengwei2024phd/03GW/try/B6_bin.5_prokka/B6_bin.5.fasta.ffn
prodigal_faa=/home/lipengwei2024phd/03GW/try/B6_bin.5_prokka/B6_bin.5.fasta.faa

defense_predictor --prokka_gff ${prodigal_gff} --prokka_ffn ${prodigal_ffn} --prokka_faa ${prodigal_faa} --device cpu --output /home/lipengwei2024phd/03GW_prokaryotes/try/B6_bin5_defense_predictor_output.csv

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
