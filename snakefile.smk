# 配置文件路径 config.yaml
configfile: "config.yaml"

# 预处理
# 自动读取 fastq 文件夹中的样本
samples = glob_wildcards(config["dir_fastq"] + "/{sample}_1.fastq.gz").sample
# 将{sample}按“-”拆分取第一个部分作为时间点
experiments = list({s.split("-")[0] for s in samples})

# 规则定义
rule all:
    input:
        #expand(config["result_dir"] + "/counts/{sample}.counts", sample=samples),
        # config["result_dir"] + "/counts/all_samples.counts",
        expand(config["result_dir"] + "/diff_genes/DEG_results_{experiment}.csv", experiment=experiments),
        expand(config["result_dir"] + "/diff_genes/volcano_plot_{experiment}.svg", experiment=experiments),
        # config["result_dir"] + "/counts/counts.csv",
        config["result_dir"] + "/diff_genes/all_upregulated.csv",
        config["result_dir"] + "/diff_genes/all_downregulated.csv",

# 1. fastp质量控制
rule fastp_qc:
    input:
        r1 = config["dir_fastq"] + "/{sample}_1.fastq.gz",
        r2 = config["dir_fastq"] + "/{sample}_2.fastq.gz"
    output:
        r1 = config["result_dir"] + "/clean/{sample}_1.clean.fastq.gz",
        r2 = config["result_dir"] + "/clean/{sample}_2.clean.fastq.gz",
        html = config["result_dir"] + "/reports/fastp/{sample}_fastp.html",
        json = config["result_dir"] + "/reports/fastp/{sample}_fastp.json"
    threads: 4
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "-j {output.json} -h {output.html} "
        "--thread {threads}"

# 2. HISAT2比对 + SAM转BAM + 排序（管道一步完成）
rule hisat2_align:
    input:
        r1 = rules.fastp_qc.output.r1,
        r2 = rules.fastp_qc.output.r2
    output:
        config["result_dir"] + "/align/{sample}.sorted.bam"  # 直接输出排序后的BAM
    params:
        index = config["hisat2_index"]
    threads: 8  # 总线程分配给HISAT2和samtools
    shell:
        "(hisat2 -x {params.index} "          # HISAT2输出到stdout
        "-1 {input.r1} -2 {input.r2} "
        "--threads {threads} | "              # 管道传递给samtools
        "samtools view -@ 2 -Sb - | "         # 使用2个线程转BAM
        "samtools sort -@ 6 -o {output})"     # 使用6个线程排序
    # 注意：总线程数= hisat2线程 + samtools view线程 + samtools sort线程 ≤ 可用核心数

# 3. featureCounts定量
rule featureCounts:
    input:
        bam = expand(config["result_dir"] + "/align/{sample}.sorted.bam", sample=samples),
        anno = config["gtf"]
    output:
        config["result_dir"] + "/counts/all_samples.counts"
    params:
        feature_type = "exon",
        gene_id = "gene_id"
    threads: 8
    shell:
        "featureCounts -p --countReadPairs "
        "-T {threads} -a {input.anno} "
        "-t {params.feature_type} -g {params.gene_id} "
        "-o {output} {input.bam}"


# 4. 编辑featureCounts输出文件all_samples.counts
rule edit_counts_file:
    input:
        config["result_dir"] + "/counts/all_samples.counts",
        design_matrix = config["design_matrix"]
    output:
        #counts = config["result_dir"] + "/counts/counts.csv", # 编辑all_samples.counts
        counts = config["result_dir"] + "/counts/{experiments}.counts.csv",  # 分实验输出
        new_design = config["result_dir"] + "/counts/{experiments}.design.csv"
    params:
        samples = samples,
        output_dir = config["result_dir"] + "/counts",
    script:
        "scripts/p4_edit_count_data.py"

# 5. DESeq2差异分析 + 火山图绘制
rule deg_analysis:
    input:
        counts = rules.edit_counts_file.output.counts,
        design = rules.edit_counts_file.output.new_design
    output:
        deg = config["result_dir"] + "/diff_genes/DEG_results_{experiments}.csv",
        plot = config["result_dir"] + "/diff_genes/volcano_plot_{experiments}.svg",
    script:
        "scripts/edgeR_analysis.R"

# 6. 差异基因筛选，从所有deg分析结果中提取同样是上调或下调的基因
rule deg_filter:
    input:
        expand(config["result_dir"] + "/diff_genes/DEG_results_{experiment}.csv", experiment=experiments),
    output:
        config["result_dir"] + "/diff_genes/all_upregulated.csv",
        config["result_dir"] + "/diff_genes/all_downregulated.csv"
    script:
        "scripts/p5_deg_filter.py"

import os
from snakemake.io import glob_wildcards
import logging

# 配置日志记录
logging.basicConfig(level=logging.INFO)

def validate_path(path):
    # 确保路径是绝对路径并且不存在路径遍历风险
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    return os.path.normpath(path)

try:
    fastq_dir = config.get("dir_fastq", "")
    if not fastq_dir:
        raise ValueError("配置中缺少 'dir_fastq' 路径")

    validated_dir = validate_path(fastq_dir)
    
    # 检查目录是否存在
    if not os.path.isdir(validated_dir):
        raise FileNotFoundError(f"目录 {validated_dir} 不存在")

    # 获取样本列表
    samples = glob_wildcards(os.path.join(validated_dir, "{sample}_1.fastq.gz")).sample
    
    if not samples:
        logging.warning("未找到任何匹配的样本文件")
        samples = []

except (KeyError, ValueError) as e:
    logging.error(f"配置错误: {e}")
    samples = []
except FileNotFoundError as e:
    logging.error(f"文件路径错误: {e}")
    samples = []
except Exception as e:
    logging.error(f"未知错误: {e}")
    samples = []

# 输出最终的样本列表
print(samples)