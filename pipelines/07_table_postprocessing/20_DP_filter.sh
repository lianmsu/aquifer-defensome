# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------
cd ${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/defense_predictor_result_analyses

#列名为：
#protein_context_id,product_accession,name,symbol,genomic_accession,start,end,strand,dna_seq,protein_seq,defense_probability,defense_log_odds

#合并所有结果
find ${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/defense_predictor/out/ -name "*.csv" | xargs cat | awk 'NR==1 || FNR!=1' > 1626mag_defensepredictor.csv

#处理文件1626mag_defensepredictor.csv，删除该文件中第十二列数值小于7.2的行，之后再删除第9列和第10列。
#!/bin/bash

# 定义输入文件和输出文件
input_file="1626mag_defensepredictor.csv"
output_file="1626mag_defensepredictor_processed.csv"

# 处理文件
awk -F',' 'BEGIN {OFS=","} {
    # 跳过表头
    if (NR == 1) {
        # 删除第 9 列和第 10 列的表头
        $9 = ""; $10 = "";
        print $1, $2, $3, $4, $5, $6, $7, $8, $11, $12;
        next;
    }
    # 检查第十二列是否大于等于 7.2
    if ($12 >= 7.2) {
        # 删除第 9 列和第 10 列
        $9 = ""; $10 = "";
        print $1, $2, $3, $4, $5, $6, $7, $8, $11, $12;
    }
}' "$input_file" > "$output_file"

echo "处理完成，结果已保存到 $output_file"

#删除文件${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/defense_predictor_result_analyses/1626mag_defensepredictor_processed.csv中除了第一行之外的protein_context_id,product_accession,name,symbol,genomic_accession,start,end,strand,defense_probability,defense_log_odds

# 定义输入文件和输出文件
input_file="${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/defense_predictor_result_analyses/1626mag_defensepredictor_processed.csv"
output_file="${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/defense_predictor_result_analyses/1626mag_defensepredictor_processed_filtered.csv"

# 使用 awk 删除除了第一行之外的指定行
awk -F, 'NR == 1 || !($1 == "protein_context_id" && $2 == "product_accession" && $3 == "name" && $4 == "symbol" && $5 == "genomic_accession" && $6 == "start" && $7 == "end" && $8 == "strand" && $9 == "defense_probability" && $10 == "defense_log_odds")' "$input_file" > "$output_file"

#统计路径${defense_predictor_out_path}下所有以.csv结尾的文件中，第十二列的数值大于7.2的行的数量并输出到文件defensepredictor_gene.csv。输出时以.csv文件名为第一列，大于7.2的数量为第二列。
cut -d "," -f12 out/B10_bin.23.fasta.csv | sort -nr | head
#!/bin/bash
defense_predictor_out_path=${LUSTRE_HOME}/prokaryote_revision/HQ_MAGs/defense_predictor/out/
# 定义输出文件
output_file="defensepredictor_genecount.csv"

# 清空输出文件（如果已存在）
> "$output_file"

# 遍历路径下所有以 .csv 结尾的文件
for file in "${defense_predictor_out_path}"/*.csv; do
    # 获取文件名（不包括路径）
    filename=$(basename "$file" .csv)
    
    # 统计第十二列大于 7.2 的行数
    count=$(awk -F',' '$12 > 7.2 {count++} END {print count+0}' "$file")
    
    # 将结果写入输出文件
    echo "$filename,$count" >> "$output_file"
done

echo "统计完成，结果已保存到 $output_file"


#对应defensefinder的gene结果到defensepredictor。根据${HOME_ROOT}/03GW/60high_quality_mags/02defensefinder/defense/1626mag_defense_finder_genes.tsv的基因ID找到基因对应的位置，用位置在DP结果中匹配start,end。
