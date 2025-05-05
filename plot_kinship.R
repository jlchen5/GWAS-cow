setwd("/Volumes/Newsmy/XJKPS162024016_KPS202412104-450例牛-全基因组重测序/analysis_result/")
# 加载必要的包
library(pheatmap)

# 读取kinship矩阵（确保路径正确）
kinship <- read.table("tassel_kinship.txt", header = FALSE, row.names = 1, skip = 3)
colnames(kinship) <- rownames(kinship)
kinship[kinship < 0] <- 0
diag(kinship) <- NA

# --------------------------
# 1. 绘制直方图（处理NA值）
# --------------------------
# 将矩阵转为向量并移除NA
kinship_values <- as.vector(as.matrix(kinship))
kinship_values <- kinship_values[!is.na(kinship_values)]

# 绘制直方图（调整分箱和颜色）
hist(
  kinship_values,
  breaks = 50,  # 增加分箱数
  xlab = "Kinship Coefficient",
  ylab = "Frequency",
  main = "Distribution of Kinship Coefficients",
  col = "skyblue",
  border = "white",
  xlim = c(0, max(kinship_values, na.rm = TRUE))  # 限制x轴范围
)
# --------------------------
# 2. 绘制热图（优化显示）
# --------------------------
# 设置热图参数（避免标签重叠）
pheatmap(
    as.matrix(kinship),
    color = colorRampPalette(c("skyblue", "white","white", "red"))(100),  # 颜色渐变
    cluster_rows = TRUE,  # 对行聚类
    cluster_cols = TRUE,  # 对列聚类
    show_rownames = FALSE,  # 隐藏行名（样本多时推荐）
    show_colnames = FALSE,  # 隐藏列名
    fontsize_row = 6,  # 行名字体大小
    fontsize_col = 6,  # 列名字体大小
    main = "Kinship Heatmap"
  )
