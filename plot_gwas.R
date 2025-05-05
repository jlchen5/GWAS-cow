library(qqman)
library(CMplot)


# 定义要处理的文件列表
gwas_files <- c("gwas_guanwei.assoc.linear",
                "gwas_weight.assoc.linear",
                "gwas_height.assoc.linear",
                "gwas_bust.assoc.linear",
                "gwas_oxhorn.assoc.logistic")

for (file in gwas_files) {
  # 从文件名中提取性状名称
  trait <- sub("gwas_(.*)\\.assoc\\..*", "\\1", file)
  
  # 读取数据
  gwas_results <- read.table(file, header = TRUE)
  
  # 数据预处理
  gwas_results$CHR <- as.numeric(gwas_results$CHR)
  gwas_results$SNP <- paste0("chr", gwas_results$CHR, "_", gwas_results$BP)
  df <- gwas_results[, c("SNP", "CHR", "BP", "P")]
  
  # 计算Bonferroni阈值
  # threshold <- 1/nrow(df[!is.na(df$BP), ])
  threshold <- 1e-5
  
  # 筛选显著SNP
  sig_snps <- subset(df, P <= threshold)
  
  # 保存显著结果
  write.csv(sig_snps, 
            file = paste0("sig_", trait, ".csv"), 
            row.names = FALSE)
  
  # 生成曼哈顿图和QQ图
  # CMplot(df,
  #        plot.type = c("m", "q"),
  #        threshold = threshold,
  #        threshold.col = c("red", "blue"),
  #        threshold.lty = c(1, 2),
  #        amplify = FALSE,
  #        file = "png",
  #        file.name = paste0(trait, "_"),
  #        dpi = 300,
  #        verbose = TRUE)
}



##############################################
################## 合并画图 ##################
##############################################

# 定义性状列表和对应文件
traits <- c(
  "guanwei" = "gwas_guanwei.assoc.linear",
  "weight" = "gwas_weight.assoc.linear",
  "height" = "gwas_height.assoc.linear",
  "bust" = "gwas_bust.assoc.linear"
)

# 初始化数据存储
combined_data <- NULL

# 循环读取并合并数据
for (i in seq_along(traits)) {
  trait_name <- names(traits)[i]
  file <- traits[trait_name]
  
  # 读取数据
  gwas <- read.table(file, header = TRUE)
  
  # 预处理
  gwas$CHR <- as.numeric(gwas$CHR)
  gwas$SNP <- paste0("chr", gwas$CHR, "_", gwas$BP)
  
  # 创建临时数据框
  temp <- data.frame(
    SNP = gwas$SNP,
    CHR = gwas$CHR,
    BP = gwas$BP,
    P = gwas$P
  )
  colnames(temp)[4] <- trait_name  # 重命名P值列为性状名称
  
  # 合并数据
  if (is.null(combined_data)) {
    combined_data <- temp
  } else {
    combined_data <- merge(combined_data, temp, by = c("SNP", "CHR", "BP"), all = TRUE)
  }
}

combined_data$CHR <- factor(combined_data$CHR, levels = 1:30)

# 设置绘图参数
CMplot(
  combined_data,
  chr.labels = as.character(1:30),  # 染色体标签
  plot.type = c("m","q"),           # 同时绘制曼哈顿图和QQ图
  multracks = F,                # 多性状绘制模式
  threshold = c(1e-5, 1e-6),    # 设置显著性阈值线
  threshold.lty = c(1, 2),      # 阈值线类型
  threshold.lwd = c(1, 1),      # 阈值线粗细
  threshold.col = c("red", "blue"),
  amplify = FALSE,              # 不放大显著点
  signal.cex = 1.5,             # 显著点大小
  signal.pch = 19,              # 统一使用实心圆点
  signal.col = c("red", "green3", "blue", "orange"), # 性状颜色
  cex = 0.5,                    # 普通点大小
  #file = "pdf",                 # 输出格式
  #file.name = "Combined_GWAS",  # 文件名前缀
  #dpi = 300,                    # 分辨率
  #width = 14,                   # 图像宽度
  #height = 10,                  # 图像高度
  main = "Combined GWAS Results",
  verbose = TRUE
)
