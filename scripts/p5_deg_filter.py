import pandas as pd

# 从snakemake获取输入和输出文件路径
from snakemake.script import snakemake
input_files = snakemake.input
output_upregulated = snakemake.output[0]
output_downregulated = snakemake.output[1]

# 本地测试相关
# input_files = ['result/diff_genes/DEG_results_exp1.csv', 'result/diff_genes/DEG_results_exp2.csv']
# output_upregulated = 'result/diff_genes/experiments_upregulated.csv'
# output_downregulated = 'result/diff_genes/experiments_downregulated.csv'

# 初始化DataFrame列表
upregulated_dfs = []
downregulated_dfs = []

logFC_threshold = 1.0

# 遍历每个deg分析结果文件
for file in input_files:
    # 读取文件
    df = pd.read_csv(file, header=0)
    df.columns = ['gene_id'] + list(df.columns[1:])
    
    # 筛选上调基因 (logFC > logFC_threshold)
    upregulated = df[df['logFC'] > logFC_threshold]
    upregulated_dfs.append(upregulated[['gene_id']])
    
    # 筛选下调基因 (logFC < -logFC_threshold)
    downregulated = df[df['logFC'] < -logFC_threshold]
    downregulated_dfs.append(downregulated[['gene_id']])

# 找出所有实验中共有的上调基因
common_upregulated = pd.concat(upregulated_dfs).groupby('gene_id').filter(lambda x: len(x) == len(input_files))
common_upregulated = common_upregulated.drop_duplicates()

# 找出所有实验中共有的下调基因
common_downregulated = pd.concat(downregulated_dfs).groupby('gene_id').filter(lambda x: len(x) == len(input_files))
common_downregulated = common_downregulated.drop_duplicates()

# 保存结果到指定的输出文件中
common_upregulated.to_csv(output_upregulated, index=False)
common_downregulated.to_csv(output_downregulated, index=False)