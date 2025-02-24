import csv,pathlib,re

# 读取PRJNA1048153_runinfo_Ac.csv的第1列和第12列，作为样本的旧ID和新ID
with open('PRJNA1048153_runinfo_Ac.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    id_map = [[row[0],row[11]] for row in reader]
# 将id_map转换为字典，键为旧ID，值为新ID
id_dict = dict(id_map)
old_ids = id_dict.keys()
print(id_dict)
dir_fastq = pathlib.Path('./fastq_data/')
list_fastq = dir_fastq.glob('*.fastq.gz')

# os修改文件名，将旧ID替换为新ID
for fastq_file in list_fastq:
    # 提取文件名中的ID部分
    id = re.search(r'(.+)_.+.fastq.gz', fastq_file.name).group(1)
    if id in old_ids:
        new_name = fastq_file.with_name(str(fastq_file.name).replace(id, id_dict[id]))
        print(f'Renaming {fastq_file.name} to {new_name.name}')
        fastq_file.rename(new_name)
        
