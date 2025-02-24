import csv
import pandas as pd
############
# 在运行snakemake之前，创建design_matrix
############

# 读取PRJNA1048153_runinfo_Ac.csv的第1列和第12列，作为样本的旧ID和新ID
with open('PRJNA1048153_runinfo_Ac.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    samples = [row[11] for row in reader]
# 使用列表推导式将samples中的字符串按下划线分割，取第二个元素作为第二列的值
samples_condition_timepoint = [[sample, sample.split('_')[1], 
                                sample.split('-')[0]] for sample in samples]
# 将结果写入csv文件
# with open('design_matrix.csv', 'w') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerows(samples_condition_timepoint)

# print(samples)

design_matrix = pd.DataFrame(samples_condition_timepoint, columns=['SampleID', 'Condition', 'Experiment'])
design_matrix = design_matrix.sort_values(by='SampleID')  # 按样本ID排序
print(design_matrix)  # 打印design_matrix
# 保存design_matrix到文件
design_matrix.to_csv('design_matrix.csv', header=True, index=False)
