import pathlib,re
import pandas as pd


# snakemake相关
from snakemake.script import snakemake
input_counts = snakemake.input[0]
design_matrix = snakemake.input.design_matrix
output_counts = snakemake.output.counts
output_dir = snakemake.params.output_dir
nuber_of_samples = len(snakemake.params.samples)


# 本地测试相关
# input_counts = 'result/counts/all_samples.counts'
# design_matrix = 'design_matrix.csv'
# output_counts = 'result/counts/counts.csv'
# output_dir = 'result/counts/'
# nuber_of_samples = 6

# 读取输入文件
count_file = pathlib.Path(input_counts)
if count_file.exists():
    counts = pd.read_csv(count_file, sep="\t", header=1)
else:
    print(f"File {count_file} does not exist")
    exit(1)

# 读取design_matrix文件
design_file = pathlib.Path(design_matrix)
if design_file.exists():
    design = pd.read_csv(design_file, sep=",")
else:
    print(f"File {design_file} does not exist")
    exit(1)

# 提取第一列gene_id和后面的样本列
count_data = counts.iloc[:, [0] + list(range(-nuber_of_samples, 0))]
print(count_data.columns)  # 打印提取后的列名


# 修改样本列名
sample_mapping = {old_name: new_name for old_name in count_data.columns[1:] for new_name in design['SampleID'] if new_name in old_name}
new_sample_name = [sample_mapping[old_name] for old_name in count_data.columns[1:]]
# 检查不匹配的样本列名，如果存在则抛出异常
unmatched_samples = set(count_data.columns[1:]) - set(sample_mapping.keys())
if unmatched_samples:
    raise ValueError(f"The following samples do not match any sample ID in the design matrix: {', '.join(unmatched_samples)}")

# 修改列名
count_data = count_data.rename(columns=sample_mapping)

# 修改样本列的顺序
count_data = count_data[[count_data.columns[0]] + design['SampleID'].to_list()]

# 保存修改后的count_data到文件
# count_data.to_csv(output_counts, sep=",", header=True, index=False)

# 按design_matrix中的Experiment分组，分别建立每个实验的counts.csv文件
for experiment in design['Experiment'].unique():
    experiment_counts = count_data[count_data.columns[0]].to_frame().join(count_data.loc[:, design[design['Experiment'] == experiment]['SampleID']])
    experiment_counts.to_csv(f"{output_dir}/{experiment}.counts.csv", sep=",", header=True, index=False)

# 保存对应的design_matrix到文件
for experiment in design['Experiment'].unique():
    experiment_design = design[design['Experiment'] == experiment]
    experiment_design.to_csv(f"{output_dir}/{experiment}.design.csv", sep=",", header=True, index=False)

