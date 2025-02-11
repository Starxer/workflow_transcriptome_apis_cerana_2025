# 判断DESeq2包是否存在，如果不存在则安装
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("DESeq2")) BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

# 读取counts数据
count_data <- read.table(snakemake@input$counts,
  header= TRUE,row.names = 1,comment.char="#",sep=",",check.names = FALSE)
design <- read.table(snakemake@input$design,header = FALSE, 
  row.names = 1, sep=",",check.names = FALSE)
colnames(design) <- c("condition", "timepoint")

# 输出count_data和design的行数和列数
print(colnames(count_data))
print(dim(design))

print(count_data[1:5, ])
print(design)

# 输出ncol(count_data)和nrow(design)
print(ncol(count_data))
print(nrow(design))

# 定义时间段
time_points <- unique(design$timepoint)

# 循环遍历每个时间段
for (time_point in time_points) {
  # 过滤出当前时间段的数据
  current_design <- subset(design, timepoint == time_point)
  current_count_data <- count_data[, rownames(current_design)]
  
  # 创建DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = current_count_data,
                                colData = current_design,
                                design = ~ condition)
  
  # 差异分析
  dds <- DESeq(dds)
  
  # 获取结果
  res <- results(dds)
  
  # 保存结果
  output_file <- paste0(snakemake@output$deg, "_", time_point, ".csv")
  write.csv(as.data.frame(res), output_file)
  
  # 绘制火山图
  res_df <- as.data.frame(res)
  # 处理NA值并创建显著性列
  res_df$significance <- ifelse(res_df$pvalue < 0.05 & !is.na(res_df$pvalue), "sig", "ns")
  
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significance), alpha = 0.6) +  # 直接使用创建好的分类列
    scale_color_manual(
      name = "Significance",  # 设置图例标题
      values = c("sig" = "red", "ns" = "gray"),  # 显式指定颜色映射
      labels = c("sig" = "Significant (pvalue < 0.05)", "ns" = "Not significant")
    ) +
    theme_bw() +
    labs(x = "log2(Fold Change)", y = "-log10(p-value)")
  
  output_plot_file <- paste0(snakemake@output$plot, "_", time_point, ".png")
  ggsave(output_plot_file, volcano_plot, width=8, height=6)
}
