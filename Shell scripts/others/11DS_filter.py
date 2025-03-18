
#要实现这个要求，我们需要对`hit_id`列的每个条目进行字符串处理，提取第二个`_`之前的部分，并将其作为新的`genome`列写入到输出文件中。以下是修改后的Python脚本：

import csv
import re
# 输入和输出文件的路径
#i = 6
input_tsv_path = f'/home/lipengwei2024phd/03GW_prokaryotes/10integron/all_integron_dss.txt'  # TSV文件的路径
output_txt_path = f'/home/lipengwei2024phd/03GW_prokaryotes/10integron/all_integron_dss2.txt'  # 输出文本文件的路径

# 要提取的列的名称
#columns_to_extract = ['hit_id', 'gene_name']
# 正则表达式模式，用于匹配直到.fasta的部分
#pattern = re.compile(r'(\S+\.fasta)')
# 正则表达式模式，用于匹配直到第一个.fasta的部分
pattern = re.compile(r'(\S+?\.fasta)')
try:
    # 打开TSV文件
    with open(input_tsv_path, 'r') as tsvfile:
        # 创建CSV阅读器
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        # 打开输出文本文件
        with open(output_txt_path, 'w') as outfile:
            outfile.write(f'sample\tgenome\tsys_id\ttype\tsubtype\tgenes_count\n')
            # 遍历TSV文件中的每一行
            for row in reader:
                # 提取hit_id列，并处理字符串以获取genome列
                sys_id = row['sys_id']

                sample = sys_id.split('_', 1)[0]
                # 使用正则表达式匹配第一个.fasta及其之前的内容
                match = pattern.search(sys_id)
                if match:
                    genome = match.group(1)  # 获取匹配的第一部分
                else:
                    genome = ""  # 如果没有匹配，设置为空字符串
                
                
                # 提取type列 subtype列 genes_count列
                type = row['type']
                subtype = row['subtype']
                genes_count = row['genes_count']
                
                # 将提取的数据以制表符分隔的方式写入新的文本文件
                # 包括新创建的genome列
                outfile.write(f'{sample}\t{genome}\t{sys_id}\t{type}\t{subtype}\t{genes_count}\n')

    print(f"Data extraction and processing completed. Check the output file: {output_txt_path}")

except FileNotFoundError:
    print(f"Error: The file {input_tsv_path} does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")


#在这个脚本中，我们使用`split('_', 1)`方法来限制分割操作只在第一个`_`处进行分割，这样`genome`变量就会包含第二个`_`之前的字符串。然后，我们将`genome`列连同`hit_id`和`gene_name`列一起写入到输出文件中，列之间用制表符分隔。

#请确保`input_tsv_path`变量的值指向正确的TSV文件路径，并且该文件具有`hit_id`和`gene_name`这两列。输出文件将包含三个列：`genome`、`hit_id`和`gene_name`。
