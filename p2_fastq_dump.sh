#!/bin/bash
# 创建FASTQ输出目录
mkdir -p fastq_data

# 使用find获取输入文件信息，并通过管道用bgzip压缩
find ./sra_data_Ac -name "*.sra" | while read -r sra; do
  fasterq-dump --progress --mem 1000MB "$sra" -O "./fastq_data/"
done

# 使用find获取压缩后的FASTQ文件，并通过管道用bgzip压缩
find ./fastq_data/ -name "*.fastq" | while read -r fastq; do
  bgzip -@ 4 "$fastq"
done