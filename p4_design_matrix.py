import csv
############
# p4_edit_count_data.py已经实现这部分的功能，此部分可以删除
############

# 读取PRJNA1048153_runinfo_Ac.csv的第1列和第12列，作为样本的旧ID和新ID
with open('PRJNA1048153_runinfo_Ac.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    samples = [row[11] for row in reader]
# 使用列表推导式将samples中的字符串按下划线分割，取第二个元素作为第二列的值
design_matrix = [[sample, sample.split('_')[1]] for sample in samples]
for i in design_matrix:
    print(i)
# 将结果写入csv文件
with open('design_matrix.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(design_matrix)

print(samples)
