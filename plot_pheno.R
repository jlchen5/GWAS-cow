# 加载必要的包
library(ggplot2)
library(gridExtra)
library(ggpubr)

# 1. 数据准备

df <- read.csv("pheno_all.csv",header = T) 

# 2. 绘制组合图表函数
plot_combined <- function(var_name) {
  # 直方图+密度曲线
  p1 <- ggplot(df, aes(x = .data[[var_name]])) +
    geom_histogram(aes(y = ..density..),
                   bins = 15,
                   fill = "skyblue",
                   color = "black") +
    geom_density(color = "red", linewidth = 1) +
    labs(title = paste(var_name, "distribution"), x = var_name)
  
  # QQ图
  p2 <- ggqqplot(df, x = var_name, color = "blue") +
    labs(title = paste(var_name, "QQ plot"))
  
  grid.arrange(p1, p2, ncol = 2)
}

# 3. 选择需要分析的数值型变量（排除分类变量）
numeric_vars <- c("age", "weight", "height", "bust", "guanwei")

# 4. 生成所有数值变量的组合图
combined_plots <- lapply(numeric_vars, plot_combined)

# 5. 正态性检验
shapiro_test <- function(var_name) {
  test_result <- shapiro.test(df[[var_name]])
  data.frame(
    Variable = var_name,
    W = round(test_result$statistic, 4),
    p.value = format.pval(test_result$p.value, eps = 0.001)
  )
}

# 执行检验并格式化结果
normality_results <- do.call(rbind, lapply(numeric_vars, shapiro_test))
print(normality_results)
