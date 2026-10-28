#!/bin/bash
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------
###########################2##########################################

#1. 生成所有DIs所在的contig的所有orf id
python3 ${HOME_ROOT}/03GW/02ANI99MAGs_defense/DIs_unknowngene/01all_DI_orf.py
#2. 提取所有DI对应的DGs的id
cat ${HOME_ROOT}/03GW/02ANI99MAGs_defense/DIs_unknowngene/defense_islands.txt | cut -f6 | grep -v ":" | grep -v '^ *$' | tr ',' '\n' | sort | uniq > 02DI_contig_DGs.txt
#注意：defensefinder在识别防御系统时会把一个DGs重复算到不同的DSs

##################################3##################################

#1. 生成所有DIs所在的contig的所有orf id
python3 ${HOME_ROOT}/03GW/02ANI99MAGs_defense/DIs_unknowngene/01all_DI_orf.py

#2. 提取所有DI对应的DGs的id
cat ${HOME_ROOT}/03GW/02ANI99MAGs_defense/DIs_unknowngene/defense_islands.txt | cut -f6 | grep -v ":" | grep -v '^ *$' | tr ',' '\n' | sort | uniq > 02DI_contig_DGs.txt
#注意：defensefinder在识别防御系统时会把一个DGs重复算到不同的DSs

#3. 从01all_DIorf.txt中删除与02DI_contig_DGs.txt共有的行
# 首先，提取并排序两个文件的每一行
sort 01all_DIorf.txt | uniq > 03all_DIorf_sorted.txt
sort 02DI_contig_DGs.txt | uniq > 03DI_contig_DGs_sorted.txt

comm -23 03all_DIorf_sorted.txt 02DI_contig_DGs.txt > 03filtered_01all_DIorf.txt

# 找出两个排序文件共有的行
#grep -Ff 03DI_contig_DGs_sorted.txt 03all_DIorf_sorted.txt > 03common_lines.txt
# 从原始文件01all_DIorf.txt中删除这些共有行
#grep -v -Ff 03common_lines.txt 01all_DIorf.txt > 03all_DIorf_filtered.txt

#4. 根据03all_DIorf_filtered.txt提取orf序列
seqkit grep -f 03filtered_01all_DIorf.txt ${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa > 04non_DG_orf.faa

#5. 在线eggnog-mapper注释

##################################

