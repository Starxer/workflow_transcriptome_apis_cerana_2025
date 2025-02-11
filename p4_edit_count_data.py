import pathlib,re
import pandas as pd


# snakemake相关
from snakemake.script import snakemake
input_file = snakemake.input[0]
output_counts = snakemake.output.counts
output_design_matrix = snakemake.output.design_matrix
nuber_of_samples = len(snakemake.params.samples)

# # 本地测试相关
# input_file = 'result/counts/all_samples.counts'
# output_counts = 'result/counts/counts.csv'
# output_design_matrix = 'result/counts/design_matrix.txt'
# nuber_of_samples = 6

count_file = pathlib.Path(input_file)
if count_file.exists():
    counts = pd.read_csv(count_file, sep="\t",header=1)
else:
    print(f"File {count_file} does not exist")
    exit(1)
# 提取第一列和后六列
count_data = counts.iloc[:, [0] + list(range(-nuber_of_samples, 0))]
print(count_data.columns)

# 创建design_matrix
new_sample_name = [re.search(r'(.+)\.sorted\.bam',pathlib.Path(i).name).groups()[0] for i in count_data.columns[1:]]
print(new_sample_name)
samples_condition_timepoint = [[sample, sample.split('_')[1], sample.split('-')[0]] for sample in new_sample_name]
design_matrix = pd.DataFrame(samples_condition_timepoint, columns=['SampleID', 'Condition', 'TimePoint'])
design_matrix = design_matrix.sort_values(by='SampleID')
print(design_matrix)
design_matrix.to_csv(output_design_matrix, sep=",", header=False,index=False)

new_sample_name_sorted = list(design_matrix['SampleID'])
# 修改列名
count_data.columns = [count_data.columns[0]] + new_sample_name
# print(count_data.iloc[:3,:5])
# 修改列顺序
count_data = count_data[[count_data.columns[0]] + new_sample_name_sorted]
# print(count_data.iloc[:3,:5])

count_data.to_csv(output_counts,sep=",",header=True,index=False)
