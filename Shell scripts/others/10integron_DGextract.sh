#以下是一个 Shell 脚本，用于根据 [merged_mags_2_2_virus_summary.tsv]文件的第一列，从 [merged_mags_2_2_defense_finder_systems.tsv]中提取对应的行。匹配规则是 [merged_mags_2_2_virus_summary.tsv]的第一列的字符串被 [merged_mags_2_2_defense_finder_systems.tsv]第一列的字符串所包含：

cd /home/lipengwei2024phd/03GW_prokaryotes/10integron

for i in 1 2 3 4 5 6 7 8 9 10;
do
# 文件路径
file1="/home/lipengwei2024phd/03GW_prokaryotes/10integron/complete_integron_contig.id"
file2="/home/lipengwei2024phd/03GW/02ANI99MAGs_defense/merged_mags_${i}_${i}/merged_mags_${i}_${i}_defense_finder_systems.tsv"
output="integron_dss_${i}.txt"

# 清空输出文件
#> "$output"

#head -n 1 "$file2" > "$output"

# 读取第一个文件的第一列
while IFS=$'\t' read -r col1 _; do
    # 在第二个文件中查找包含第一个文件第一列字符串的行
    grep -F "$col1" "$file2" >> "$output"
done < "$file1"

echo "Data extraction and processing completed. Check the output file: $output"
done

cat integron_dss_* > all_integron_dss.txt
