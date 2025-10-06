library(ggplot2)
library(ggrepel)

setwd('/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/01_pop_structure')

# Read Variance Explained from .eigenval
eigenval_file <- "/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/01_pop_structure/CostaRica.thin1000.pop.K2.eigenval"

# Read the eigenval file
eigenval <- read.table(eigenval_file, header = F, sep = "\t", comment.char = "#")

# Ensure column names are consistent
colnames(eigenval) <- c("eval", "eval_percent")  # First column = eigenvalues, Second column = variance %

# === Read PCA Eigenvec File ===
pca_file <- "/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/01_pop_structure/CostaRica.thin1000.pop.K2.eigenvec"

# Read PCA data
PCA <- read.table(pca_file, header = TRUE)
PCA$Sample <- gsub("^.*/|.bam$", "", PCA$SampleName)  # Clean sample names

# === Specify PCs to Plot ===
x_PC <- "PC2"  # Replace with the desired x-axis PC
y_PC <- "PC3"  # Replace with the desired y-axis PC

# Extract Variance Explained% for the selected PCs
x_PC_index <- as.numeric(gsub("PC", "", x_PC))  # Get PC index (e.g., PC1 -> 1)
y_PC_index <- as.numeric(gsub("PC", "", y_PC))

x_var_explained <- eigenval$eval_percent[x_PC_index]
y_var_explained <- eigenval$eval_percent[y_PC_index]

# === Create the Plot ===
# Dynamically adjust labels for axes
x_label <- paste0(x_PC, " (Variance Explained: ", round(x_var_explained, 2), "%)")
y_label <- paste0(y_PC, " (Variance Explained: ", round(y_var_explained, 2), "%)")

# refine PCA plotting
ggplot(data = PCA, mapping = aes_string(x = x_PC, y = y_PC, 
                                        colour = "Group", 
                                        fill = "Group",
                                        shape = "Cluster", 
                                        label = "Sample")) +
  geom_point(size = 2) +
  ggtitle(paste("PCA:", x_PC, "vs", y_PC))+
  labs(x = x_label, y = y_label) +
  geom_text_repel(size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_colour_manual(values = c("#F8766D", "#00B0F6")) +
  scale_fill_manual(values = c("#F8766D", "#00B0F6")) +
  theme_bw() +
  theme(text = element_text(size=15))

# Save the plot with dynamic filename
ggsave(paste0("CostaRica.thin1000.pop.K2.", x_PC, "_", y_PC, ".png"), height=6, width=7)

