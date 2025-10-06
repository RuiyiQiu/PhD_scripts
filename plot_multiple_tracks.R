library(png)
library(scales)

setwd("/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/04_summary")

inPATH="/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/04_summary/signals"
win <- "/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/Rciliata.h1.chr.500kb.windows"
gene <- paste(inPATH, "gene.list", sep = "/")
SNP <- paste(inPATH, "SNP.list", sep = "/")
TYBA <- paste(inPATH, "TYBA.list", sep = "/")
TRC_1 <- paste(inPATH, "TRC_1.list", sep = "/")
TE <- paste(inPATH, "TE.list", sep = "/")
Fst <- paste(inPATH, "Fst.list", sep = "/")
pi_all <- paste(inPATH, "pi_all.list", sep = "/")
pi_sub1 <- paste(inPATH, "pi_n.list", sep = "/")
pi_sub2 <- paste(inPATH, "pi_s.list", sep = "/")
rho_all <- paste(inPATH, "rho_all.list", sep = "/")
rho_sub1 <- paste(inPATH, "rho_n.list", sep = "/")
rho_sub2 <- paste(inPATH, "rho_s.list", sep = "/")
c_CO <- paste(inPATH, "Rciliata_CO.list", sep = "/") 
n_CO <- paste(inPATH, "Rnervosa_CO.list", sep = "/") #pollen CO by Meng


windows <- read.table(win)
gene_list <- read.table(gene)
SNP_list <- read.table(SNP)
TYBA_list <- read.table(TYBA)
TRC_1_list <- read.table(TRC_1)
TE_list <- read.table(TE)
Fst_list <- read.table(Fst)
pi_all_list <- read.table(pi_all)
pi_sub1_list <- read.table(pi_sub1)
pi_sub2_list<- read.table(pi_sub2)
rho_all_list <- read.table(rho_all)
rho_sub1_list <- read.table(rho_sub1)
rho_sub2_list <- read.table(rho_sub2)
c_CO_list <- read.table(c_CO)
n_CO_list <- read.table(n_CO)


df <- cbind(windows, gene_list, SNP_list,  
            TYBA_list, TRC_1_list, TE_list,
            Fst_list, 
            pi_all_list, pi_sub1_list, pi_sub2_list, 
            rho_all_list, rho_sub1_list, rho_sub2_list,
            c_CO_list, n_CO_list)
colnames(df) <- c('chr', 'win_s', 'win_e',
                  'gene', 'SNP', 
                  'TYBA', 'TRC_1', 'TE',
                  'Fst', 
                  'pi_all', 'pi_sub1', 'pi_sub2', 
                  'rho_all', 'rho_sub1', 'rho_sub2', 
                  'ciliata_CO', 'nervosa_CO')

repeat_counts <- scan("signals/TRC_1.list")
repeat_df <- cbind(windows, count = repeat_counts)
mask_df <- subset(repeat_df, count>0)
colnames(mask_df) <- c("chr", "start", "end", "count")

# ----- functions -----
scale_values <- function(x, scaling_factor){
  if(length(scaling_factor)==1) return(x/scaling_factor)
  # (x-min(x))/(max(x)-min(x))
  else if(length(scaling_factor)==2) return( (x-scaling_factor[2]) / (scaling_factor[1]-scaling_factor[2]))
  else print("Wrong number of scaling  factors!")
} # Feature scaling is used to bring all values into the range [0,1]. This is also called unity-based normalization. 


plot_multiple_tracks <- function(chr, max_length, lw, text_size, scale_factors){
  chrsize <- chrsizes[which(chrsizes$V1==chr),]$V2
  max_length <- chrsize
  this_df <- df[df$chr==chr,]
  this_mask_df <- mask_df[mask_df$chr==chr,]

  scale_factors <- list(max(this_df$gene, na.rm = T),              #1-gene
                        max(this_df$SNP, na.rm = T),               #2-SNP
                        c(max(this_df$TYBA, na.rm = T),
                          min(this_df$TYBA, na.rm = T)),           #3-TYBA
                        c(max(this_df$TRC_1, na.rm = T),
                          min(this_df$TRC_1, na.rm = T)),          #4-TRC_1
                        c(max(this_df$TE, na.rm = T),
                          min(this_df$TE, na.rm = T)),             #5-TE
                        c(max(this_df$Fst, na.rm = T),
                          min(this_df$Fst, na.rm = T)),            #6-Fst
                        c(max(this_df$pi_all, na.rm = T), 
                          min(this_df$pi_all, na.rm = T)),         #7-pi_all
                        c(max(this_df$pi_sub1, na.rm = T), 
                          min(this_df$pi_sub1, na.rm = T)),        #8-pi_sub1
                        c(max(this_df$pi_sub2, na.rm = T), 
                          min(this_df$pi_sub2, na.rm = T)),        #9-pi_sub2
                        max(this_df$rho_all, na.rm = T),           #10-rho_all
                        max(this_df$rho_sub1, na.rm = T),          #11-rho_sub1
                        max(this_df$rho_sub2, na.rm = T),          #12-rho_sub2
                        max(this_df$ciliata_CO, na.rm = T),        #13-ciliata CO
                        max(this_df$nervosa_CO), na.rm = T         #14-nervosa CO
  )
  
  ########## Subfig 1: genes and SNPs ##########
  
  ### Compute genes, SNPs denisty
  this_gene <- scale_values(this_df$gene, scale_factors[[1]])
  this_snp  <- scale_values(this_df$SNP, scale_factors[[2]])
  
  ### Plot genes, SNPs denisty in one track
  
  plot(x=this_df$win_s, y=this_gene,
       xlim=c(1, max_length), ylim = c(0,1), xlab="", ylab="",
       type="n", axes=F)
  
  # masked repeat regions
  for(i in seq_len(nrow(this_mask_df))) {
    rect(xleft  = this_mask_df$start[i],
         ybottom = 0,
         xright = this_mask_df$end[i],
         ytop   = 1,
         col    = rgb(0.5, 0.5, 0.5, 0.3),  # semi-transparent grey
         border = NA)
  }
  
  # gene density
  polygon(c(this_df$win_s,chrsize,0),
          c(this_gene, 0, 0),
          col=rgb(0.2,0.1,0.5,0.2) , border=F)
  # SNPs
  lines(x=this_df$win_s, y=this_snp, col="navy", lwd=lw)
  
  axis(1, at=c(seq(0, chrsize, 20000000), chrsize),
       labels = F, col = "black", lwd=lw)
  axis(2, at = c(0, 1), labels = c(0, 1), lwd=lw,
       cex.axis=text_size)
  ########## Subfig 2: TYBA and TRC_1 repeats ##########
  
  ### Compute TYBA and TRC_1
  this_TYBA <- scale_values(this_df$TYBA, scale_factors[[3]])
  this_TRC_1  <- scale_values(this_df$TRC_1, scale_factors[[4]])
  this_TE  <- scale_values(this_df$TE, scale_factors[[5]])
  
  ### Plot TYBA and TRC_1 in one track
  # TYBA
  plot(x=this_df$win_s, y=this_TYBA,
       xlim=c(1, max_length), ylim = c(0,1), xlab="", ylab="",
       type="n", axes=F)
  lines(x=this_df$win_s, y=this_TYBA, col="purple", lwd=lw)
  # TE
  polygon(c(this_df$win_s,chrsize,0),
          c(this_TE, 0, 0),
          col=alpha("purple", 0.2) , border=F)
  # TRC_1 
  lines(this_df$win_s, this_TRC_1, col=alpha("#1874CD", 0.8), lwd=lw) #dodgerblue3
  
  axis(1, at=c(seq(0, chrsize, 20000000), chrsize),
       labels = F, col = "black", lwd=lw)
  axis(2, at = c(0, 1), labels = c(0, 1), lwd=lw,
       cex.axis=text_size)
  
  ########## Subfig 3: Fst ##########
  
  ### normalization
  this_Fst <- scale_values(this_df$Fst, scale_factors[[6]])
  
  ### Plot frame
  plot(x=this_df$win_s, y=this_Fst,
       xlim=c(1, max_length), ylim = c(0,1), xlab="", ylab="",
       type="n", axes=F)
  
  # masked repeat regions
  for(i in seq_len(nrow(this_mask_df))) {
    rect(xleft  = this_mask_df$start[i],
         ybottom = 0,
         xright = this_mask_df$end[i],
         ytop   = 1,
         col    = rgb(0.5, 0.5, 0.5, 0.3),  # semi-transparent grey
         border = NA)
  }
  
  # Fst
  lines(x=this_df$win_s, y=this_Fst, col="firebrick3", lwd=lw)
  axis(1, at=c(seq(0, chrsize, 20000000), chrsize),
       labels = F, col = "black", lwd=lw)
  axis(2, at = c(0, 1), labels = c(0, 1), lwd=lw,
       cex.axis=text_size)
  
  ########## Subfig 4: pi ##########
  this_pi_all <- scale_values(this_df$pi_all, scale_factors[[7]])
  this_pi_sub1 <- scale_values(this_df$pi_sub1, scale_factors[[8]])
  this_pi_sub2 <- scale_values(this_df$pi_sub2, scale_factors[[9]])
  
  plot(x=this_df$win_s, y=this_pi_all,
       xlim=c(1, max_length), ylim = c(0,1), xlab="", ylab="",
       type="n", lwd=lw, axes=F)
  
  # masked repeat regions
  for(i in seq_len(nrow(this_mask_df))) {
    rect(xleft  = this_mask_df$start[i],
         ybottom = 0,
         xright = this_mask_df$end[i],
         ytop   = 1,
         col    = rgb(0.5, 0.5, 0.5, 0.3),  # semi-transparent grey
         border = NA)
  }
  
  # pi_all
  lines(this_df$win_s, this_pi_all, col="forestgreen", lwd=lw)
  # pi_sub1
  lines(this_df$win_s, this_pi_sub1, col=alpha("#6495ED", 0.6), lwd=lw-0.6)
  # pi_sub2
  lines(this_df$win_s, this_pi_sub2, col=alpha("#959c35", 0.6), lwd=lw-0.6)
  
  
  axis(1, at=c(seq(0, chrsize, 20000000), chrsize),
       labels = F, col = "black", lwd=lw)
  axis(2, at = c(0, 1), labels = c(0, 1), lwd=lw,
       cex.axis=text_size)
  
  ########## Subfig 5: rho ##########
  this_rho_all <- scale_values(this_df$rho_all, scale_factors[[10]])
  this_rho_sub1 <- scale_values(this_df$rho_sub1, scale_factors[[11]])
  this_rho_sub2 <- scale_values(this_df$rho_sub2, scale_factors[[12]])
  
  plot(x=this_df$win_s, y=this_rho_all,
       xlim=c(1, max_length), ylim = c(0,1), xlab="", ylab="",
       type="n", lwd=lw, axes=F)
  
  # masked repeat regions
  for(i in seq_len(nrow(this_mask_df))) {
    rect(xleft  = this_mask_df$start[i],
         ybottom = 0,
         xright = this_mask_df$end[i],
         ytop   = 1,
         col    = rgb(0.5, 0.5, 0.5, 0.3),  # semi-transparent grey
         border = NA)
  }
  
  # rho_all
  lines(x=this_df$win_s, y=this_rho_all, col="darkgrey", lwd=lw)
  # rho_sub1
  lines(x=this_df$win_s, y=this_rho_sub1, col=alpha("#F8766D", 0.6), lwd=lw-0.6)
  # rho_sub2
  lines(x=this_df$win_s, y=this_rho_sub2, col=alpha("#00B0F6", 0.6), lwd=lw-0.6)
  
  axis(1, at=c(seq(0, chrsize, 20000000), chrsize),
       labels = F, col = "black", lwd=lw)
  axis(2, at = c(0, 1), labels = c(0, 1), lwd=lw,
       cex.axis=text_size)
  
  ########## Subfig 5: CO landscape ##########
  this_c_CO <- scale_values(this_df$ciliata_CO, scale_factors[[13]])
  this_n_CO <- scale_values(this_df$nervosa_CO, scale_factors[[14]])
  
  plot(x=this_df$win_s, y=this_c_CO,
       xlim=c(1, max_length), ylim = c(0,1), xlab="", ylab="",
       type="n", lwd=lw, axes=F)
  # CO 
  #lines(x=this_df$win_s, y=this_n_CO, col="#E9967A", lwd=lw)
  lines(x=this_df$win_s, y=this_c_CO, col="#563429", lwd=lw)
  
  axis(1, at=c(seq(0, chrsize, 20000000), chrsize),
       labels=round(c(seq(0, chrsize, 20000000), chrsize)/1000000),
       lwd=lw, cex.axis=text_size,  mgp=c(0,2,0))
  axis(2, at = c(0, 1), labels = c(0, 1), lwd=lw,
       cex.axis=text_size)
  
  # close
  #dev.off()
  return(0)
}

# ----- end of function definition -----

#####
# Read some other related files
chrsizes <- read.table("/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/chr_size.txt", 
                       sep = "\t")

#####

# --------------------- Figure plotting ----------------------
chr='Chr5_h1'
pdf(paste(chr,".pdf", sep=""),
    family="Helvetica", height=6, width=4)
par(mai = c(0.2, 0.4, 0, 0.1)); # margin: bottom, left, top, right
m <- matrix(1:7, nrow=7, byrow=F)
layout(m, heights = c(1,1,1,1,1,1,0.2), widths = 3.6)
plot_multiple_tracks(chr=chr, lw=2, text_size=2)
dev.off()

for (chr in chrsizes$V1) {
  pdf(paste(chr,".pdf", sep=""),
      family="Helvetica", height=6, width=4)
  par(mai = c(0.2, 0.4, 0, 0.1)); # margin: bottom, left, top, right
  m <- matrix(1:7, nrow=7, byrow=F)
  layout(m, heights = c(1,1,1,1,1,1,0.2), widths = 3.6)
  plot_multiple_tracks(chr=chr, lw=2, text_size=2)
  dev.off()
}
