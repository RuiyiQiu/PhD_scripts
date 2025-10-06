library(GenomicRanges)
library(dplyr)

# Inputs
args <- commandArgs(trailingOnly=TRUE)
win_size <- as.numeric(args[1])
locs_file <- args[2]
LDhat_file <- args[3]
chr_name <- args[4]
out_file <- args[5]


# Read LDhat + locs
LD <- read.table(LDhat_file, header=TRUE) %>% filter(Loci != 0)
locs <- scan(locs_file, skip=1)
LD$start <- locs[-length(locs)]*1e3
LD$end   <- locs[-1]*1e3

# Make GRanges
gr <- GRanges(chr_name, IRanges(LD$start, LD$end), rho=LD$Mean_rho)

# Windows
chr_len <- max(LD$end)
wins <- GRanges(chr_name, 
                IRanges(seq(1, chr_len, by=win_size),
                        pmin(seq(win_size, chr_len+win_size, by=win_size), 
                             chr_len)))

# Overlaps → weighted mean rho
ov <- findOverlaps(gr, wins)
df <- data.frame(win=subjectHits(ov), rho=mcols(gr)$rho[queryHits(ov)],
                 overlap=width(pintersect(gr[queryHits(ov)], wins[subjectHits(ov)])))
res <- df %>% group_by(win) %>%
  summarise(LDhat_rho=weighted.mean(rho, overlap), .groups="drop")

# Add coords
out <- data.frame(CHROM=chr_name, BIN_START=start(wins), BIN_END=end(wins))
out$LDhat_rho <- res$LDhat_rho[match(seq_along(wins), res$win)]

write.table(out, out_file, sep="\t", quote=F, row.names=F, col.names=F)


