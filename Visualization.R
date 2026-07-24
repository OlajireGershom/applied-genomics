# Visualization -----------------------------------------------------------
#Step 1: Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

#Step 2: Load saved objects
dds <- readRDS("data/dds_object.rds")
metadata <- read.csv("data/metadata/sample_info.tsv", stringsAsFactors = FALSE, sep = "\t")

#Set condition as a factor with Control as the reference level
metadata$Condition <- factor(metadata$Condition, levels = c("Control", "FluD4", "FluD8"))
rownames(metadata) <- metadata$SampleID

#Step 3: Visualization
##1 VSI transformation (for all visualizations)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

## 2 Principal Component Analysis Plot
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = Condition, label = name)) +
  geom_point(size = 6, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.5) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_manual(values = c("Control" = "#2E6DA4", 
                                "FluD4" = "#E67E22", 
                                "FluD8" = "#C0392B")) +
  theme_bw(base_size = 13) +
  ggtitle("PCA of VSI-normalized counts",
          subtitle = "Mouse Cerebellum - Influenza infection time course")

#save plot
ggsave("Results/Figures/PCA_plot.png", pca_plot, width = 7, height = 5, dpi = 300)
cat("PCA plot saved! \n")

##Step 3 Sample distance heatmap
samp_dists <- dist(t(vst_mat))
samp_mat <- as.matrix(samp_dists)

annot_col <- data.frame(
  Condition = metadata[colnames(vst_mat), "Condition"],
  row.names = colnames(vst_mat)
)

annot_colors <- list(
  Condition = c("Control" = "#2E6da4",
                "FluD4" = "#e67e22",
                "FluD8" = "#c0392b"
                )
)

jpeg("Results/Figures/Sample_distance_heatmap.jpg", width = 400, height = 400, 
     units = "px")

pheatmap(
  samp_mat,
  clustering_distance_rows = samp_dists,
  clustering_distance_cols = samp_dists,
  annotation_col = annot_col,
  annotation_colors = annot_colors,
  colot = colorRampPalette(brewer.pal(9, "Blues"))(100),
  main = "Sample-to-sample distances"
)

dev.off()
cat("Sample distance heatmap saved! \n")

##4. Volcano plots - one per comparison
comparisons <- list(
  list(file = "Results/dataFluD4vsControlall_genes.tsv",
       title = "FluD4 vs Control",
       name = "FluD4_vs_Control"),
  list(file = "Results/dataFluD8vsControlall_genes.tsv",
       title = "FluD8 vs Control",
       name = "FluD8_vs_Control"),
  list(file = "Results/dataFluD8vsFluD4all_genes.tsv",
       title = "FluD8 vs FluD4",
       name = "FluD8_vs_FluD4")
)

for (comp in comparisons) {
  res_df <- read.delim(comp$file)
  
  v <- EnhancedVolcano(
    res_df,
    lab = res_df$gene_id,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    title = comp$title,
    subtitle = "DESeq2 | ashr shrinkage | FDR < 0.05, |LFC| >1",
    caption = paste("Total genes tested", nrow(res_df)),
    col = c("grey70", "grey40", "#2E6DA4", "#C0392B"),
    colAlpha = 0.7,
    pointSize = 2,
    labSize = 3,
    legendPosition = "bottom",
    drawConnectors = TRUE,
    widthConnectors = 0.4
  )
  
  ggsave(paste0("Results/Figures/Volcano_plot", comp$name, ".png"),
         v, width = 8, height = 7, dpi = 300)
  cat("Volcano plot saved!", comp$name, "\n")
}

## 5. Heatmap of all significant DEGs combined
#Pool all unique significant gene IDs across comparisons
sig_d4 <- read.delim("Results/dataFluD4vsControlSig_DEGs.tsv")
sig_d8 <- read.delim("Results/dataFluD8vsControlSig_DEGs.tsv")

all_sig_genes <- unique(c(sig_d4$gene_id, sig_d8$gene_id))
cat("Total unique significant genes across comparisons:", length(all_sig_genes), "\n")

#Subset VST matrix to significant genes
heat_map <- vst_mat[all_sig_genes, ]

# Z-score scale each row so colours reflect relative expression
heat_scaled <- t(scale(t(heat_map)))
jpeg("Results/Figures/Heatmap_sig_DEGs.jpg", width = 480, height = 480, units = "px")
pheatmap(
  heat_scaled,
  annotation_col = annot_col,
  annotation_colors = annot_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Significant DEGs ~ z-score of VST counts"
  )                 
dev.off()
cat("DEG heatmap saved! \n")
