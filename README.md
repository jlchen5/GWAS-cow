# 牛的GWAS分析

## ➡️ GWAS全基因组关联分析全流程教程：

[1.GWAS：原理与目的](https://www.jianshu.com/p/e060c55283c4)  
[2.GWAS：流程与试验设计](https://www.jianshu.com/p/2949f17715b6)  
[3.1 GWAS：表型鉴定与记录的基本原则和原始数据处理](https://www.jianshu.com/p/eb893d4fe177)  
[3.2 GWAS：最佳线性无偏估计量——BLUE值计算（多年单点有重复）](https://www.jianshu.com/p/f4f1b5b75830)  
[3.4 GWAS：遗传力计算](https://www.jianshu.com/p/d116831897a2)  
[4. 标记的开发和分型](https://www.jianshu.com/p/f0464a1afaeb)  
[4.2 基因型数据描述性统计](https://www.jianshu.com/p/5b4a9566eee7)  
[5. GWAS：群体结构——Admixture](https://www.jianshu.com/p/ef1c8dcddf96)  
[6. GWAS：主成分分析——GCTA](https://www.jianshu.com/p/933eefea3fcc)  
[7.1 GWAS：系统进化树——MEGA](https://www.jianshu.com/p/5a969bcdc1fe)  
[7.2 GWAS：系统进化树美化——ITOL](https://www.jianshu.com/p/6081dda6445f)  
[8. GWAS：亲缘关系——TASSEL&GCTA](https://www.jianshu.com/p/b24ba1d50448)  
[9.1 GWAS：关联分析](https://www.jianshu.com/p/a22c72af4dea)  
[9.2 GWAS：关联分析——TASSEL（GLM/MLM/CMLM)](https://www.jianshu.com/p/22edbe46bf7d)  
[9.3 GWAS：关联分析——EMMAX](https://www.jianshu.com/p/0f39ff5a7643)  
[9.4 GWAS：显著性阈值确定——GEC](https://www.jianshu.com/p/055daa26d5c6)  
[9.5 GWAS显著SNP筛选及曼哈顿图绘制](https://www.jianshu.com/p/e8e88c54966d)  
[10.GWAS：LD decay(LD衰减）—— PopLDdecay](https://www.jianshu.com/p/50a9c66fbd2a)  
[11.GWAS：确定候选区间](https://www.jianshu.com/p/905bf1f3a798)  
[12.GWAS：候选基因挖掘](https://www.jianshu.com/p/bc836e179347)  
[CMplot报错missing value where TRUE/FALSE needed](https://www.jianshu.com/p/19784fa6bbdd)  
 


## ➡️ 分析数据结构

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

### 2.1 **统计质控前的SNP染色体分布**

使用 `bcftools` 统计原始VCF文件中各染色体的SNP数量：

```
# 统计原始VCF中各染色体的SNP数量
grep -v "^#" vcf.finalSNP.vcf | cut -f 1 | sort | uniq -c

2781070 1
2323018 2
2003432 3
2176331 4
1932744 5
2113230 6
1853304 7
1900852 8
1796157 9
1804853 10
1791273 11
1733600 12
1445094 13
1370579 14
1611254 15
1365610 16
1243811 17
1095006 18
1022131 19
1314969 20
1210699 21
1031071 22
1147312 23
1123929 24
739724 25
948282 26
889019 27
852580 28
1020097 29
```

### 2.2 **数据质控（过滤低质量SNP）**

使用 `PLINK` 进行质控（按需调整参数）：

```
# 转换为PLINK格式并质控
plink --vcf vcf.finalSNP.vcf \
			--geno 0.1  # 过滤缺失率>10%的SNP \
			--maf 0.01  # 过滤次要等位基因频率<5% \
			--hwe 1e-6  # 过滤哈迪-温伯格平衡检验p<1e-6 \
			--make-bed  \
			--out vcf.finalSNP_filtered  \
			--allow-extra-chr --chr-set 29

# 将质控后的数据转回VCF格式（便于后续统计）
plink --bfile vcf.finalSNP_filtered --recode vcf --out vcf.finalSNP_filtered --allow-extra-chr --chr-set 29
```

### 2.3 **统计质控后的SNP染色体分布**

```
# 统计质控后VCF中各染色体的SNP数量
grep -v "^#" vcf.finalSNP_filtered.vcf | cut -f 1 | sort | uniq -c

821043 1
658403 2
588284 3
636202 4
476786 5
582112 6
580904 7
581838 8
486603 9
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



```R
### 基于杂合度的方法
# R脚本（需安装 hierfstat）
library(hierfstat)
data <- read.plink("plink_format")
basic_stats <- basic.stats(data$genotypes)
Ne_het <- 1 / (2 * (1 - basic_stats$Hs))
write.csv(Ne_het, "Ne_heterozygosity.csv")
```

### maf计算

```
plink --bfile finalSNP_gwas_filtered --freq --allow-extra-chr --chr-set 29 --out finalSNP_gwas_filtered_maf_freq 
```





## 4 ROH分析

```sh
plink --bfile vcf.finalSNP_filtered --allow-extra-chr --chr-set 29\
  --homozyg   --out vcf.finalSNP_filtered_roh
```





## 5 样本名修改

```sh
$ bcftools  query -l vcf.finalSNP_filtered.vcf  |sed 's/_/\t/' |cut -f 1  rename
$ bcftools reheader -s rename vcf.finalSNP_filtered.vcf -o finalSNP_filtered_rename.vcf
```





## 6 聚类树绘制



```
## 利用VCF2Dis生成距离矩阵
VCF2Dis -i finalSNP_filtered_rename.vcf -o finalSNP.mat

## 然后用fast2.0在线工具绘制，生成nwk文件
```



## 7 GWAS 分析

- 文件过滤：去除了重复值和缺失值的样本

  ```sh
  # 找出covar.txt中已有的样本。
  cat covar.txt |cut -f 1 |sed '1d'  filtered_sample.txt
  
  # 提取已有样本
  vcftools --vcf finalSNP_filtered_rename.vcf --recode --recode-INFO-all --stdout  --keep filtered_sample.txt finalSNP_gwas_filtered.vcf
  
  # 后续分析使用 finalSNP_gwas_filtered。
  
  ### 文件格式转换
  plink --vcf finalSNP_gwas_filtered.vcf --make-bed --out finalSNP_gwas_filtered --allow-extra-chr --chr-set 29
  ```

- weight

  ```
  plink --bfile finalSNP_gwas_filtered --pheno pheno.txt --pheno-name weight --covar covar.txt --covar-name age --linear hide-covar --ci 0.95 --allow-extra-chr --chr-set 29 --out gwas_weight --allow-no-sex
  
  ```

- height

  ```
  plink --bfile finalSNP_gwas_filtered --pheno pheno.txt --pheno-name height --covar covar.txt --covar-name age --linear hide-covar --ci 0.95 --allow-extra-chr --chr-set 29 --out gwas_height --allow-no-sex
  ```

- bust

  ```
  plink --bfile finalSNP_gwas_filtered --pheno pheno.txt --pheno-name bust --covar covar.txt --covar-name age --linear hide-covar --ci 0.95 --allow-extra-chr --chr-set 29 --out gwas_bust --allow-no-sex
  ```

- guanwei

  ```
  plink --bfile finalSNP_gwas_filtered --pheno pheno.txt --pheno-name guanwei --covar covar.txt --covar-name age --linear hide-covar --ci 0.95 --allow-extra-chr --chr-set 29 --out gwas_guanwei --allow-no-sex
  ```

- oxhorn，hide``-``covar指的是不要对我没加入的协变量进行分析，所以这里去除这个参数；由于 `oxhorn` 是二分类变量（0/1），使用 `--logistic` 代替 `--linear；表格为oxhorn.txt的前几行 (1：有角；2：无角)

  ```
  plink --bfile finalSNP_gwas_filtered --pheno oxhorn.txt --pheno-name oxhorn --logistic hide-covar --ci 0.95 --allow-extra-chr --chr-set 29 --out gwas_oxhorn --allow-no-sex 
  ```

  | FID     | IID     | oxhorn |
  | ------- | ------- | ------ |
  | R252901 | R252901 | 2      |
  | Q017506 | Q017506 | 2      |
  | Q018278 | Q018278 | 1      |
  | R253789 | R253789 | 1      |
  | Q007402 | Q007402 | 1      |

- 可视化：曼哈顿图和qq图

  ```
  # 加载包
  library(qqman)
  library(CMplot)
  
  ### 读取结果文件:
  #### gwas_guanwei.assoc.linear
  #### gwas_weight.assoc.linear
  #### gwas_height.assoc.linear
  #### gwas_bust.assoc.linear
  #### gwas_oxhorn.assoc.logistic
  
  setwd("/Volumes/Newsmy/XJKPS162024016_KPS202412104-450例牛-全基因组重测序/analysis_result/")
  gwas_results <- read.table("gwas_weight.assoc.linear", header=TRUE)
  gwas_results$CHR <- as.numeric(gwas_results$CHR)
  gwas_results$SNP <- paste0("chr", gwas_results$CHR, "_", gwas_results$BP)
  
  # QQ图
  # qq(gwas_results$P)
  
  # manhattan图
  df <- gwas_results[,c("SNP", "CHR", "BP", "P")]
  
  # 画图
  threshold <- 1/nrow(df[!is.na(df$BP),])
  
  CMplot(df, threshold = threshold ,
         amplify = F, file = "png", plot.type=c("m","q"))
  ```

  

## 8 亲缘关系计算

~~~
#### TASSEL进行亲缘关系计算
#计算IBS亲缘关系矩阵
run_pipeline.pl -Xmx32g -Xms512m -importGuess finalSNP_gwas_filtered.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export tassel_kinship.txt -exportType SqrMatrix
~~~

- 在R中可视化

  ```
  library("pheatmap")
  kinship<-read.table("tassel_kinship.txt",header = F,row.names = 1,skip = 3)
  colnames(kinship) <- row.names(kinship)
  kinship[kinship < 0] <- 0
  diag(kinship) <- NA
  pdf("Histogram of Kinship.pdf",width=10,height=8)
  hist_data <- hist(as.matrix(kinship), xlab = "Kinship", col = "grey", main = "Histogram of Kinship")
  pdf("heapmap of kinship.pdf",width=8,height=8)
  pheatmap(kinship, fontsize_row = 0.3, fontsize_col = 0.3)
  ```

  

