library(ggplot2)

pca <- read.table("finalSNP_PCA.eigenvec", header=F)
pca_eigenval <- read.table("finalSNP_PCA.eigenval", 
                           header = FALSE, col.names = "eigenvalue")

# 计算每个主成分的方差解释比例
variance_percent <- pca_eigenval$eigenvalue / sum(pca_eigenval$eigenvalue) * 100
# 取前两个主成分的解释比例（PC1和PC2）
pc1_var <- variance_percent[1]
pc2_var <- variance_percent[2]


# 基础绘图
p <- ggplot(pca, aes(x = V3, y = V4)) +
  geom_point(
    size = 3,          # 增大点的大小
    alpha = 0.7,       # 设置透明度
    color = "#2E75B6", # 使用更柔和的蓝色
    shape = 16         # 实心圆点
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1_var),
    y = sprintf("PC2 (%.1f%%)", pc2_var),
    title = "PCA Analysis of SNP Variants",
    subtitle = "PC1 vs PC2 Projection"
  ) +
  theme_classic() +    # 简洁的经典主题
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)  # 添加边框
  ) 
# 显示图形
print(p)
