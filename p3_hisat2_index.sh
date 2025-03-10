genome=./genome/GCF_029169275.1/GCF_029169275.1_AcerK_1.0_genomic.fna
gtf=./genome/GCF_029169275.1/genomic.gtf
index=./genome/index_hisat2/genome_Ac

mkdir -p ./genome/index_hisat2

# 提取剪接位点（splice sites）
time hisat2_extract_splice_sites.py $gtf > splicesites.txt

# 提取外显子信息（exons）
time hisat2_extract_exons.py $gtf > exons.txt

time hisat2-build -p 8 \
    --ss splicesites.txt \
    --exon exons.txt \
    $genome \
    $index