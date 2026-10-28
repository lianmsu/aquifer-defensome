#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J bac2.0.0_1
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500GB
# #SBATCH --reservation=<your-reservation>   # uncomment if your cluster requires one
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

# 开始计时
start_time=$(date +%s)


source activate defensefinder2.0.0
#########################################################################
# File Name: GO_DefenseFinder.sh
# Date: 17/4/2024
#########################################################################
cd ${HOME_ROOT}/03GW/02ANI99_MAGs_defense2.0.0/seperate
i="1"
mkdir -p merged_mags_${i}_${i}

fnas=${HOME_ROOT}/03GW/02ANI99MAGs_defense/merged_mags_${i}_${i}.fasta
magdefense=${HOME_ROOT}/03GW/02ANI99_MAGs_defense2.0.0/seperate/merged_mags_${i}_${i}
echo "当前目录为：$(pwd)"
echo "Running DefenseFinder"

defense-finder run "$fnas" --preserve-raw -o "$magdefense"

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
