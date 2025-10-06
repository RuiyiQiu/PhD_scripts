library(LEA)
library(mapplots)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggforce)

#arguments
args <- commandArgs(trailingOnly = T)
input.file <- arg[1]
#manual
input.file <- '/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/CostaRica.thin1000.vcf'

WD <- paste0(dirname(input.file),'/01_pop_structure')
#filename <- basename(input.file)
setwd(WD)


# get a genotype file ('geno' format) as input
vcf2geno(input.file, "input.geno")

# calculate
obj.at = snmf("input.geno", K = 1:10, ploidy = 2, entropy = T,
              CPU = 4, project = "new")
#obj.at = load.snmfProject("input.snmfProject")
png(filename="cross-entropy.png")
plot(obj.at, col = "blue4", cex = 1.4, pch = 19)
dev.off()
#pick up best K
cross_entropy_values <- sapply(1:10, function(k) {
  min(cross.entropy(obj.at, K = k))
})
best_K <- which.min(cross_entropy_values)
#manual:
#best_K <- 2
cat("Best K based on minimum cross-entropy is:", best_K, "\n")

# visualize the matrix of ancestry coefficients
qmatrix = Q(obj.at, K = best_K) 

# sample IDs
id_order <- paste0(dirname(input.file),'/id_order.list')
ids <- read.table(id_order, header=F)
rownames(qmatrix) <- ids$V1
qmatrix_long <- melt(qmatrix)

# barplot representation
ggplot(qmatrix_long, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", width = 1, color='black', linewidth=0.1) +  # Stacked bars
  scale_fill_manual(values = c("#00B0F6", "#F8766D")) +  # Custom colors
  theme_minimal() +
  labs(x = "Sample ID", y = "Admixture Coefficients", fill = "Ancestry") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(expand=c(0,0))
ggsave('Admixture.barplot.png', height=3, width=9)
