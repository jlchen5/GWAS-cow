# 加载包
library(clusterProfiler)
library(org.Bt.eg.db)
library(tidyverse)
library(enrichplot)

# 1. 读取数据并提取唯一基因符号
##### sig_height_anno.txt
##### sig_guanwei_anno.txt
##### sig_oxhorn_anno.txt
##### sig_weight_anno.txt
##### sig_bust_anno.txt

data <- read.delim("./sig_bust_anno.txt", header = FALSE)  # 替换为实际文件名
genes <- unique(data$V7)  # 第7列包含基因符号（如APP, LIPI等）
cat("需要分析的基因数量:", length(genes), "\n")

# 关键修正：清理基因符号中的空格
clean_genes <- genes %>%
  str_trim() %>%          # 移除首尾空格
  str_replace_all("\\s+", "") %>%  # 移除所有空格
  toupper() %>%           # 转换为大写（数据库通常使用大写）
  unique()
genes <-  clean_genes

# 2. 基因标识符转换（基因符号 → ENTREZID → KEGG ID）
symbol_to_entrez <- bitr(genes, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = "org.Bt.eg.db")

entrez_ids <- symbol_to_entrez$ENTREZID
cat("成功转换的基因数量:", length(entrez_ids), "\n")

# 3. KEGG通路富集分析
kk <- enrichKEGG(gene = entrez_ids,
                 organism = "bta",       # 牛: Bos taurus
                 keyType = "kegg",        # 使用KEGG内部分析
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")

# 4. 结果可视化
# 基本结果表格
results <- as.data.frame(kk)
write.csv(results, "bust_kegg_results.csv", row.names = FALSE)

# 富集分析条形图
barplot(kk, showCategory = 15, font.size = 10, 
        title = "KEGG Pathway Enrichment (Bos taurus)")

# 富集分析点图
dotplot(kk, showCategory = 15, font.size = 10,title = "Bust KEGG Pathway Enrichment")


