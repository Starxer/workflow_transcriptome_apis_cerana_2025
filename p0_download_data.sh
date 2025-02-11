#!/bin/bash
# 查询BioProject PRJNA1048153下的所有SRA编号
esearch -db sra -query "PRJNA1048153" | efetch -format runinfo > PRJNA1048153_runinfo.csv

# 导出Apis cerana的条目
#(head -n 1 PRJNA1048153_runinfo.csv && grep -F 'A.c.' PRJNA1048153_runinfo.csv)> PRJNA1048153_runinfo_Ac.csv
grep -F 'A.c.' PRJNA1048153_runinfo.csv | grep -v 'Am' > PRJNA1048153_runinfo_Ac.csv

# 提取SRA编号列表（例如SRR/ERR/DRR开头的编号）
#cut -d ',' -f 1 PRJNA1048153_runinfo_Ac.csv | grep -v "Run" > PRJNA1048153_sra_list.txt
cut -d ',' -f 1 PRJNA1048153_runinfo_Ac.csv > PRJNA1048153_sra_list.txt


# 查看前5个SRA编号确认
head -n 5 PRJNA1048153_sra_list.txt

# 创建下载目录
mkdir -p sra_data_Ac

# 使用prefetch批量下载
time prefetch -O ./sra_data_Ac --option-file PRJNA1048153_sra_list.txt
