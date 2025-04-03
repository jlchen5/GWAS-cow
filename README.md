# 牛的GWAS分析

- 分析数据结构

```
(base)  /Volumes/Newsmy/XJKPS162024016_KPS202412104-450例牛-全基因组重测序 [06:48:25]
jialechen$ tree -L 2
.
├── analysis_result
└── result
    ├── 00.rawdata
    ├── 01.cleandata
    ├── 02.QC
    ├── 03.map
    ├── 04.vcf
    │   ├── vcf.finalINDEL.vcf
    │   ├── vcf.finalINDEL.vcf.idx
    │   ├── vcf.finalSNP.vcf
    │   └── vcf.finalSNP.vcf.idx
    ├── 05.snp_indel_stat
    └── XJKPS162024016_KPS202412104-450个牛-重测序报告.zip
```



## 1 群体结构分析



### 1.1 群体结构 PCA 图

- 格式转换

~~~
### 1. 数据质控（PLINK）
# 牛基因组需添加 --chr-set 29 no-xy
plink --vcf vcf.finalINDEL.vcf \
      --allow-extra-chr \
      --chr-set 29 no-xy \
      --const-fid \
      --make-bed \
      --out vcf.finalINDEL
      
plink --vcf vcf.finalSNP.vcf \
      --allow-extra-chr \
      --chr-set 29 no-xy \
      --const-fid \
      --make-bed \
      --out vcf.finalSNP



### 2. 群体结构分析（PLINK PCA）
plink --bfile vcf.finalSNP\
      --allow-extra-chr \
      --chr-set 29 no-xy \
      --pca 3 \
      --out finalSNP_PCA
      
plink --bfile vcf.finalINDEL\
      --allow-extra-chr \
      --chr-set 29 no-xy \
      --pca 3 \
      --out finalINDEL_PCA

### 3. R脚本绘制PCA图
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



~~~



## 2 质控前后SNP的分布

### **统计质控前的SNP染色体分布**

使用 `bcftools` 统计原始VCF文件中各染色体的SNP数量：

```
# 统计原始VCF中各染色体的SNP数量
grep -v "^#" vcf.finalSNP.vcf | cut -f 1 | sort | uniq -c
```

### **数据质控（过滤低质量SNP）**

使用 `PLINK` 进行质控（按需调整参数）：

```
# 转换为PLINK格式并质控
plink --vcf vcf.finalSNP.vcf --geno 0.1 --maf 0.01 --hwe 1e-6 --make-bed --out vcf.finalSNP_filtered --allow-extra-chr --chr-set 29

# 将质控后的数据转回VCF格式（便于后续统计）
plink --bfile vcf.finalSNP_filtered --recode vcf --out vcf.finalSNP_filtered --allow-extra-chr --chr-set 29
```

### **统计质控后的SNP染色体分布**

```
# 统计质控后VCF中各染色体的SNP数量

grep -v "^#" vcf.finalSNP_filtered.vcf | cut -f 1 | sort | uniq -c

(base) ➜  analysis_result  grep -v "^#" vcf.finalSNP_filtered.vcf | cut -f 1 | sort | uniq -c
821043 1
521943 10
472248 11
456292 12
419865 13
405667 14
491719 15
364738 16
338281 17
301232 18
275383 19
658403 2
356813 20
345685 21
292026 22
324335 23
333979 24
216208 25
252730 26
235071 27
233879 28
301007 29
588284 3
636202 4
476786 5
582112 6
580904 7
581838 8
486603 9
```



## 3 遗传多样性分析

- **计算观测杂合度（Ho）、期望杂合度（He）、核苷酸多样性（π）**

```
# 使用 vcftools 计算
vcftools --vcf vcf.finalSNP_filtered.vcf \
  --het \
  --out finalSNP_filtered_hetero

vcftools --vcf vcf.finalSNP_filtered.vcf\
  --site-pi \
  --out finalSNP_filtered_nucl_diver

# 使用 PLINK 计算杂合度
plink --bfile vcf.finalSNP_filtered \
  --het --allow-extra-chr --chr-set 29 \
  --out finalSNP_filtered_het
```

- **有效群体大小（Ne）估计**

```
### 基于连锁不平衡（LD）的方法
# 使用 PLINK 计算 LD
plink --bfile vcf.finalSNP_filtered --allow-extra-chr --chr-set 29 \
  --r2 \
  --ld-window 1000 \
  --ld-window-kb 1000 \
  --ld-window-r2 0 \
  --out finalSNP_filtered_ld

# 使用 NeEstimator (需安装并配置)
# 生成 .geno 文件后运行
NeEstimator -i input.geno -o Ne_LD_estimate.txt
```




### 基于杂合度的方法
# R脚本（需安装 hierfstat）
library(hierfstat)
data <- read.plink("plink_format")
basic_stats <- basic.stats(data$genotypes)
Ne_het <- 1 / (2 * (1 - basic_stats$Hs))
write.csv(Ne_het, "Ne_heterozygosity.csv")
```





## ROH分析

```sh
plink --bfile vcf.finalSNP_filtered --allow-extra-chr --chr-set 29\
  --homozyg   --out vcf.finalSNP_filtered_roh
```





## 样本名修改

```sh
$ bcftools  query -l vcf.finalSNP_filtered.vcf  |sed 's/_/\t/' |cut -f 1  > rename
$ bcftools reheader -s rename vcf.finalSNP_filtered.vcf -o finalSNP_filtered_rename.vcf
```





## 聚类树绘制



```
## vcf文件转phy文件：转换数据，坑！样本名一定要少于10个字符
vcf2phylip.py -i  finalSNP_filtered_rename.vcf

## 利用VCF2Dis生成距离矩阵
VCF2Dis -i finalSNP_filtered_rename.vcf -o finalSNP.mat
```

