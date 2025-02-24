## 简介
从NCBI上下载项目PRJNA1048153数据，进行中蜂感染中蜂囊状幼虫病病毒（Sacbrood Virus）前后的转录组分析

## 环境软件准备
创建conda环境
```sh
conda create -n rna-seq python pandas hisat2 samtools entrez-direct sra-tools fastp snakemake subread R -c conda-forge
```

## 参考基因组准备
- 下载连接https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029169275.1/
- 下载好后解压在工作目录的genome目录中。如果没有放到对应目录，就修改p3_hisat2_index.sh中genome和gtf的路径为实际的路径

## 运行分析流程
- 首先把p开头的脚本按0-4的顺序运行
- 运行snakemake -s snakefile.smk工作流