library(Hmisc) 

setwd("/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/04_summary")

inPATH="/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/04_summary/signals"
win <- "/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/Rciliata.h1.chr.500kb.windows"
gene <- paste(inPATH, "gene.list", sep = "/")
SNP <- paste(inPATH, "SNP.list", sep = "/")
TYBA <- paste(inPATH, "TYBA.list", sep = "/")
TRC_1 <- paste(inPATH, "TRC_1.list", sep = "/")
TE <- paste(inPATH, "TE.list", sep = "/")
Fst <- paste(inPATH, "Fst.list", sep = "/")
TajimaD <- paste(inPATH, "TajimaD.list", sep = "/")
pi_all <- paste(inPATH, "pi_all.list", sep = "/")
pi_sub1 <- paste(inPATH, "pi_n.list", sep = "/")
pi_sub2 <- paste(inPATH, "pi_s.list", sep = "/")
rho_all <- paste(inPATH, "rho_all.list", sep = "/")
rho_sub1 <- paste(inPATH, "rho_n.list", sep = "/")
rho_sub2 <- paste(inPATH, "rho_s.list", sep = "/")
c_CO <- paste(inPATH, "Rciliata_CO.list", sep = "/")  #pollen CO by Meng



windows <- read.table(win)
gene_list <- read.table(gene)
SNP_list <- read.table(SNP)
TYBA_list <- read.table(TYBA)
TRC_1_list <- read.table(TRC_1)
TE_list <- read.table(TE)
Fst_list <- read.table(Fst)
TajimaD_list <- read.table(TajimaD)
pi_all_list <- read.table(pi_all)
pi_sub1_list <- read.table(pi_sub1)
pi_sub2_list<- read.table(pi_sub2)
rho_all_list <- read.table(rho_all)
rho_sub1_list <- read.table(rho_sub1)
rho_sub2_list <- read.table(rho_sub2)
c_CO_list <- read.table(c_CO)


df <- cbind(windows, gene_list, SNP_list,  
            TYBA_list, TRC_1_list, TE_list,
            Fst_list, TajimaD_list,
            pi_all_list, pi_sub1_list, pi_sub2_list, 
            rho_all_list, rho_sub1_list, rho_sub2_list,
            c_CO_list)
colnames(df) <- c('chr', 'win_s', 'win_e',
                  'gene', 'SNP', 
                  'TYBA', 'TRC_1', 'TE',
                  'Fst', 'TajimaD',
                  'pi_all', 'pi_sub1', 'pi_sub2', 
                  'rho_all', 'rho_sub1', 'rho_sub2', 
                  'ciliata_CO')

# Keep only the statistics of interest
stats <- c('SNP', 'Fst', 'TajimaD',
           'pi_all', 'pi_sub1', 'pi_sub2', 
           'rho_all', 'rho_sub1', 'rho_sub2', 
           'ciliata_CO')
df_num <- df[ , stats]

#Remove rows with missing values
df_num <- na.omit(df_num)

#Correlation matrix with p-values
res <- rcorr(as.matrix(df_num), type = "pearson")

# rcorr returns:
#   res$r : matrix of correlation coefficients
#   res$P : matrix of p-values
#   res$n : matrix of pairwise sample sizes

# Save results
write.table(res$r, file = "correlation_coefficients.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# plot heatmap
library(pheatmap)

cols <- colorRampPalette(c("navy","white","firebrick3"))(200)
# create evenly spaced breaks from -1 to 1 so that 0 is the midpoint
brks <- seq(-1, 1, length.out = 201)   # one more than number of colors
png(filename = 'correlation.heatmap.png', width = 500, height = 500)
pheatmap(res$r,
         color          = cols,
         breaks         = brks,        # fixes scale from -1 to 1
         display_numbers = T,
         number_format   = "%.2f",
         number_color = "black",
         cluster_rows    = F,
         cluster_cols    = F)
dev.off()

png(filename = 'pi-SNP.points.png', width = 500, height = 500)
plot(df$pi_all, df$SNP, xlab="pi", ylab="SNP density")
dev.off()

png(filename = 'pi-Fst.points.png', width = 500, height = 500)
plot(df$pi_all, df$Fst, xlab="pi", ylab="Fst")
dev.off()

png(filename = 'pi-rho.points.png', width = 500, height = 500)
plot(df$pi_all, df$rho_all, xlab="pi", ylab="rho")
dev.off()

