#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J integron_98
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=250GB
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

# 开始计时
start_time=$(date +%s)



#########################################################################
cd ${HOME_ROOT}/03GW/06mag_integron
source activate integron_finder
#########################################################################

for i in 9 8;
do
mags="${HOME_ROOT}/03GW/02ANI99MAGs_defense/merged_mags_${i}_${i}.fasta"
integron_finder --local-max --func-annot ${mags} --cpu $SLURM_NTASKS --outdir ${HOME_ROOT}/03GW/06mag_integron --gbk --pdf
done

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
