# 假设txt_files是一个包含所有txt文件路径的列表
#txt_files = ['file1.txt', 'file2.txt', '...']

# 读取所有文件并将它们合并为一个DataFrame
#df_list = [pd.read_csv(file, sep='\t') for file in txt_files]
#combined_df = pd.concat(df_list, ignore_index=True)

# 接下来使用上面的代码对combined_df进行处理
# ...

import pandas as pd

# 读取txt文件（这里假设只有一个文件，如果有多个，可以循环读取）
txt_files = ['/home/lipengwei2024phd/03GW_prokaryotes/10integron/all_integron_dss2.txt']
df_list = [pd.read_csv(file, sep='\t') for file in txt_files]
df = pd.concat(df_list, ignore_index=True)


# 确保数据列的数据类型正确
df['genome'] = df['genome'].astype(str)
df['type'] = df['type'].astype(str)

# 使用groupby对genome和gene_name进行分组，并计数每个组合的数量
gene_counts = df.groupby(['genome', 'type']).size().reset_index(name='counts')

# 找出所有唯一的genome和gene_name
unique_genomes = df['genome'].unique()
unique_gene_names = df['type'].unique()

# 创建一个以genome为行索引，gene_name为列索引的矩阵
matrix = pd.pivot_table(gene_counts, index='genome', columns='type', values='counts', fill_value=0)

# 如果需要，可以将矩阵保存为CSV文件
matrix.to_csv('all_integron_dss_counts_matrix.csv')

# 打印结果
#print(gene_counts)

# 如果需要，可以将结果保存到新的txt文件中
#gene_counts.to_csv('gene_counts.txt', sep='\t', index=False)
