library(ggtree)
library(treeio)
tree<-read.newick("/Users/jialechen/Desktop/Analysis_work/czy/gwas-cow-450/cow-tree.nwk")
ggtree(tree,layout = "circular")+
  geom_tiplab(size=3)


# Load ggplot2
library(ggplot2)

chromosomes <- c("1",  "2", "3", "4", "5", "6", "7", "8", "9","10", "11",
                 "12", "13", "14", "15", "16", "17", "18","19","20", "21",
                 "22", "23", "24", "25", "26", "27", "28", "29" )
values <- c(2781070,2323018,2003432,2176331,1932744,2113230,1853304,1900852,
            1796157,1804853,1791273,1733600,1445094,1370579,1611254,1365610,
            1243811,1095006,1022131,1314969,1210699,1031071,1147312,1123929,
            739724,948282,889019,852580,1020097)

values_filter <- c(821043,658403,588284,636202,476786,582112,580904,581838,
                   486603,521943,472248,456292,419865,405667,491719,364738,
                   338281,301232,275383,356813,345685,292026,324335,333979,
                   216208,252730,235071,233879,301007)
# Plot
bp <- barplot(values, names.arg = chromosomes, las = 2, 
              main = "Chromosome Values", xlab = "Chromosome", ylab = "SNP numbers")
text(x = bp, y = values, labels = values, srt = 90, adj = c(1, 0.5), cex = 0.8, xpd = TRUE)

# Plot
bp2 <- barplot(values_filter, names.arg = chromosomes, las = 2, 
              main = "Chromosome Values after filtered", xlab = "Chromosome", ylab = "SNP numbers after filtered")
text(x = bp, y = values, labels = values, srt = 90, adj = c(1, 0.5), cex = 0.8, xpd = TRUE)




