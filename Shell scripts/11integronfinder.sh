#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J integron_98
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=250GB

# 开始计时
start_time=$(date +%s)



#########################################################################
cd /home/lipengwei2024phd/03GW/06mag_integron
source activate integron_finder
#########################################################################

for i in 9 8;
do
mags="/home/lipengwei2024phd/03GW/02ANI99MAGs_defense/merged_mags_${i}_${i}.fasta"
integron_finder --local-max --func-annot ${mags} --cpu $SLURM_NTASKS --outdir /home/lipengwei2024phd/03GW/06mag_integron --gbk --pdf
done

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
