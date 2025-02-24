# 判断edgeR包是否存在，如果不存在则安装
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
if (!requireNamespace("edgeR")) BiocManager::install("edgeR")

library(edgeR)
library(ggplot2)

# 读取counts数据
count_data <- read.table(snakemake@input$counts,
  header = TRUE, row.names = 1, comment.char = "#", sep = ",", check.names = FALSE
)
design <- read.table(snakemake@input$design,
  header = TRUE,
  row.names = 1, sep = ",", check.names = FALSE
)
colnames(design) <- c("condition", "experiment")

# # 输出count_data和design的行数和列数
# print(colnames(count_data))
# print(dim(design))

# print(count_data[1:5, ])
# print(design)

# # 输出ncol(count_data)和nrow(design)
# print(ncol(count_data))
# print(nrow(design))


# 过滤低计数基因
keep <- rowSums(count_data) >= 5
count_data <- count_data[keep, ]
print(dim(count_data))

# 过滤存在表达为0的基因，即去除存在表达为0的行
keep <- apply(count_data, 1, function(x) any(x > 1))
count_data <- count_data[keep, ]
print("filtered: ")
print(dim(count_data))


# edgeR差异表达分析
y <- DGEList(counts = count_data, group = design$condition)
y <- calcNormFactors(y)

bcv <- 0.1
et <- exactTest(y, dispersion = bcv^2, pair = c("Blank", "AcSBV"))
res <- topTags(et, n = Inf, sort.by = "logFC")$table


print(dim(res))
print(dim(et))
# print(et)

gene1 <- decideTests(et, p.value = 0.05)
summary(gene1)

# 输出结果到文件,确保保存列名为第一列并添加表头
write.csv(res, snakemake@output$deg, row.names = TRUE)

# 绘制火山图
res$significance <- ifelse(res$PValue < 0.05 & !is.na(res$PValue), "sig", "ns")

volcano_plot <- ggplot(res, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(
    name = "Significance",
    values = c("sig" = "red", "ns" = "gray"),
    labels = c("sig" = "Significant (pvalue < 0.05)", "ns" = "Not significant")
  ) +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10 FDR")

output_plot_file <- paste0(snakemake@output$plot)
ggsave(output_plot_file, volcano_plot, width = 8, height = 6)
