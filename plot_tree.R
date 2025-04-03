library(ggtree)
library(treeio)
tree<-read.newick("./p_dis_mat_fastme/p_dis_mat_fastme-tree.nwk")
ggtree(tree,layout = "roundrect")+
  geom_tiplab(size=3)


# Load ggplot2
library(ggplot2)

chromosomes <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", 
                 "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", 
                 "7A", "7B", "7D")
values <- c(103253, 142867, 61724, 133574, 138410, 73883, 81156, 142995, 44853, 
            85717, 69301, 24668, 76073, 107902, 40195, 88747, 133461, 46187, 
            114031, 99998, 56378)

# Plot
bp <- barplot(values, names.arg = chromosomes, las = 2, 
              main = "Chromosome Values", xlab = "Chromosome", ylab = "SNP numbers")

# Adding values on top of each bar
# text(x = bp, y = values, labels = values, pos = 3, srt = 90, cex = 1, xpd = TRUE)
text(x = bp, y = values, labels = values, srt = 90, adj = c(1, 0.5), cex = 0.8, xpd = TRUE)


# plot hr
het <- read.table("wheat_160k_xinjiang_25_id.het", head=TRUE)

het$HET_RATE = (het$N_SITES - het$O.HOM.)/het$N_SITES
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()

###### plot mds
dat <- read.table( "plink.mds", header = T)
pop <- dat[ ,1]
mds1 <- dat[ ,4]
mds2 <- dat[ ,5]
par(mar=c(5, 4, 4, 5) + 0.1) # 设置边距为下、左、上、右（单位：行数）
plot(mds1, mds2, col = as.factor(pop), xlab = "C1", ylab = "C2", pch = 20, main = "MDS plot")
#legend("bottomleft", bty = "n", cex = 0.8, legend = unique(pop), border = F, fill = unique(as.factor(pop))) 
#加上图注，并以第一列的FamilyID对个体标记不同的颜色。
legend("right", inset=c(-0.2,0), # 负inset值会将图例移动到绘图区域之外
       legend = unique(pop), fill = unique(as.factor(pop)), 
       bty = "n", cex = 0.8, xpd=TRUE)


