# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------


#一、准备自己的Acr蛋白序列

#自己的Acr蛋白################
#首先是defensefinder识别的所有prophage中的164个anti_CRISPR蛋白。
cd ${HOME_ROOT}/03GW/30Acr_tree/my_acr
grep "Anti_CRISPR" ${HOME_ROOT}/03GW/05mag_genomad_anti/all_right_virus_anti.txt > prophage_164Anti_CRISPR.txt

#手动添加标题行

#根据prophage_164Anti_CRISPR.txt的hit_id列从${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa中提取蛋白序列
awk '{print $2}' prophage_164Anti_CRISPR.txt > prophage_164Anti_CRISPR_hit_ids.txt
seqkit grep -f prophage_164Anti_CRISPR_hit_ids.txt ${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa > prophage_164Anti_CRISPR.faa

##################################

#第二是被细菌spacer靶向的prophage中找到的27个anti_CRISPR蛋白。
sed 's/|.*//' ${HOME_ROOT}/03GW/12minced/2virus_seq/targeted.prophage.id > updated_targeted.601prophage.id

while IFS=$'\t' read -r col1 _; do
    # 在第二个文件中查找包含第一个文件第一列字符串的行
    grep -F "$col1" "${HOME_ROOT}/03GW/30Acr_tree/my_acr/prophage_164Anti_CRISPR.txt" >> "spacer_targeted_prophagewithAnti_CRISPR.txt"
done < "${HOME_ROOT}/03GW/30Acr_tree/my_acr/updated_targeted.601prophage.id"

#根据spacer_targeted_prophagewithAnti_CRISPR.txt的hit_id列从${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa中提取蛋白序列
awk '{print $2}' spacer_targeted_prophagewithAnti_CRISPR.txt > prophage_27Anti_CRISPR_hit_ids.txt
seqkit grep -f prophage_27Anti_CRISPR_hit_ids.txt ${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa > prophage_27Anti_CRISPR.faa

###############################################

#二、准备建树

#0、把自己的Acr蛋白序列与PaCRISPR database中的Acr序列合并建树，用数据库作为参考。
cd ${HOME_ROOT}/03GW/30Acr_tree
cat ${HOME_ROOT}/03GW/30Acr_tree/my_acr/prophage_164Anti_CRISPR.faa ${HOME_ROOT}/03GW/30Acr_tree/PaCRISPR/Acr_sequences.fasta > ${HOME_ROOT}/03GW/30Acr_tree/all_Acr.faa


#1、先用muscle把protein序列对齐形成多序列比对（MSA）
muscle -align ${HOME_ROOT}/03GW/30Acr_tree/all_Acr.faa -output 503acr_aln.afa

#2、然后用trimal修剪
trimal -in 503acr_aln.afa -gt 0.5 > 503acr_aln.trim.afa
#-gt 0.5代表删去有50%以上gap的位点。-gt 0.5这个是对应这篇文章方法中的参数 https://www.nature.com/articles/s41467-023-38400-0

#3、再用fasttree对修建后的MSA进行建树
FastTree 503acr_aln.trim.afa > aln.trim.tree



######################上面的建树过程混入了aca、vcrx90-93
#自己的Acr蛋白################
#首先是defensefinder识别的所有prophage中的164个anti_CRISPR蛋白。
cd ${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise
grep "Anti_CRISPR" ${HOME_ROOT}/03GW/05mag_genomad_anti/all_right_virus_anti.txt > prophage_164Anti_CRISPR.txt

#提取prophage_164Anti_CRISPR.txt中只为acr的部分

grep "acr" prophage_164Anti_CRISPR.txt | grep -v "aca" > prophage_41acr.txt

#手动添加标题行

#根据prophage_41acr.txt的hit_id列从${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa中提取蛋白序列
awk '{print $2}' prophage_41acr.txt > prophage_41acr_hit_ids.txt
seqkit grep -f prophage_41acr_hit_ids.txt ${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa > prophage_41acr.faa

##################################

#第二是被细菌spacer靶向的prophage中找到的5个anti_CRISPR蛋白。
sed 's/|.*//' ${HOME_ROOT}/03GW/12minced/2virus_seq/targeted.prophage.id > updated_targeted.601prophage.id

while IFS=$'\t' read -r col1 _; do
    # 在第二个文件中查找包含第一个文件第一列字符串的行
    grep -F "$col1" "${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise/prophage_41acr.txt" >> "spacer_targeted_prophagewithAnti_CRISPR.txt"
done < "${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise/updated_targeted.601prophage.id"

#根据spacer_targeted_prophagewithAnti_CRISPR.txt的hit_id列从${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa中提取蛋白序列
awk '{print $2}' spacer_targeted_prophagewithAnti_CRISPR.txt > prophage_5Anti_CRISPR_hit_ids.txt
seqkit grep -f prophage_5Anti_CRISPR_hit_ids.txt ${HOME_ROOT}/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa > prophage_5Anti_CRISPR.faa


#二、准备建树

#0、把自己的Acr蛋白序列与PaCRISPR database中的Acr序列合并建树，用数据库作为参考。
cd ${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise
cat ${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise/prophage_41acr.faa ${HOME_ROOT}/03GW/30Acr_tree/PaCRISPR/Acr_sequences.fasta > ${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise/all_Acr.faa


#1、先用muscle把protein序列对齐形成多序列比对（MSA）
muscle -align ${HOME_ROOT}/03GW/30Acr_tree/my_acr_revise/all_Acr.faa -output 380acr_aln.afa

#2、然后用trimal修剪
trimal -in 380acr_aln.afa -gt 0.5 > 380acr_aln.trim.afa
#-gt 0.5代表删去有50%以上gap的位点。-gt 0.5这个是对应这篇文章方法中的参数 https://www.nature.com/articles/s41467-023-38400-0

#3、再用fasttree对修建后的MSA进行建树
FastTree 380acr_aln.trim.afa > aln.trim.tree
