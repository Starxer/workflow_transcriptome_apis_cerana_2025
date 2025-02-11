import os
from snakemake.io import glob_wildcards
from snakemake.common import configfile
import logging


# 配置日志记录
logging.basicConfig(level=logging.INFO)

def validate_path(path):
    # 确保路径是绝对路径并且不存在路径遍历风险
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    return os.path.normpath(path)

try:
    fastq_dir = configfile.load_configfile("config.yaml")["dir_fastq"]
    print(fastq_dir)
    if not fastq_dir:
        raise ValueError("配置中缺少 'dir_fastq' 路径")

    validated_dir = validate_path(fastq_dir)
    
    # 检查目录是否存在
    if not os.path.isdir(validated_dir):
        raise FileNotFoundError(f"目录 {validated_dir} 不存在")

    # 获取样本列表
    samples = glob_wildcards(os.path.join(validated_dir, r"{sample}.fastq.gz")).sample
    # 使用集合推导式提取唯一前缀和后缀，并确保分割安全性
    prefixes = {s.split('_', 1)[0] for s in samples}
    suffixes = {s.split('_', 1)[1] for s in samples}

    # 转换为排序列表
    sample_names = sorted(prefixes)
    sample_paired_suffix = sorted(suffixes)  # 直接获取排序后的后缀列表

    # 打印结果（保持原代码的输出顺序）
    print(sample_names)
    print(sample_paired_suffix)  # 已排序，与原代码最终结果一致

    # 重新组合样本列表（保持后缀顺序）
    samples = [[f"{name}_{suffix}.fastq.gz" for suffix in sample_paired_suffix] for name in sample_names]
    
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