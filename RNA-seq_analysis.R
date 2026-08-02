if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ashr")

BiocManager::install("EnhancedVolcano")

BiocManager::install("clusterProfiler")

BiocManager::install("GO.db")
BiocManager::install("GO.db")

BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(DESeq2)
library(ashr)
library(tidyverse)
library(GO.db)

files <- list.files(path = "counts/",
                    pattern = "\\.tabular$",
                    full.names = TRUE
                    )
print(files)

#Read each count file
count_list <- lapply (files, function(f) {
  read.delim(f, header = TRUE)
})

#Merge all count tables by Geneid
count_matrix <- Reduce(function(x,y) merge(x, y, by = "Geneid", all = TRUE),
                       count_list)

#set gene IDs as row names
rownames(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL #remove the redundant Geneid column

#Convert to integer matrix (required by DeSEq2)
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

View(count_matrix)
dim(count_matrix)
colSums(count_matrix) #library sizes per sample

#Confirm it loaded properly
metadata <- read.csv("data/metadata/sample_info.tsv", stringsAsFactors = FALSE, sep = "\t")
print(metadata)

#Set condition as a factor with Control as the reference level
metadata$Condition <- factor(metadata$Condition, levels = c("Control", "FluD4", "FluD8"))

#Set SampleID as rownames (DESeq uses this to match colData to countData)
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- NULL

#Confirm alignment before proceeding
cat("Count matrix columns:\n"); print(colnames(count_matrix))
cat("Metadata rows:\n"); print(rownames(metadata))
cat("Do they match?", all(colnames(count_matrix) == rownames(metadata)), "\n")


#Step4...create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design =  ~ Condition #three-level factor: Control, FluD4, FluD8
)

#Set reference level (Control = baseline for comparison)
dds$Condition <- relevel(dds$Condition, ref = "Control")

#Step 5: Prefilter low-count genes
#Remove genes with fewer than 10 counts across all samples combined
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")

#Step 6: Run DESeq2
dds <- DESeq(dds)

#Step 7: Extract all results--pair wise comparisons
#FluD4 vs Control
res_D4 <- results(dds,
                  contrast = c("Condition", "FluD4", "Control"),
                  alpha = 0.05)
#FluD8 vs Control
res_D8 <- results(dds,
                  contrast = c("Condition", "FluD8", "Control"),
                  alpha = 0.05)
#FluD8 vs FluD4
res_D8vsD4 <- results(dds,
                  contrast = c("Condition", "FluD8", "FluD4"),
                  alpha = 0.05)

cat("≡ FluD4 vs Control ≡\n"); summary(res_D4)
cat("≡ FluD8 vs Control ≡\n"); summary(res_D8)
cat("≡ FluD8 vs FluD4 ≡\n"); summary(res_D8vsD4)

#Step 8: LFC shrinkage and export for each comparison
resultsNames(dds)  #run this first to see exact coefficient names
shrink_and_export <- function(dds, contrast, filename_prefix) {
  #Shrinkage using ashr for contrasts
  res_s <- lfcShrink(dds,
                     contrast = contrast,
                     type = "ashr")
  #Full results Table
  df_all <- as.data.frame(res_s) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)
  
  write_tsv(
    df_all,
    paste0("C:/Users/GERSHOM/Desktop/data/rnaseq/applied-genomics/data", filename_prefix, "all_genes.tsv")
  )
  
  #Significant DEGs
  df_sig <- df_all %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  
  write_tsv(
    df_sig,
    paste0("C:/Users/GERSHOM/Desktop/data/rnaseq/applied-genomics/data", filename_prefix, "Sig_DEGs.tsv")
  )
  
  cat(filename_prefix, "significant DEGs:", nrow(df_sig),
      "(up:", sum(df_sig$log2FoldChange > 0),
      "|down:", sum(df_sig$log2FoldChange < 0), ")\n")
  return(invisible(df_all))
}

#Run for each comparison
shrink_and_export(dds, c("Condition", "FluD4", "Control"), "FluD4vsControl")
shrink_and_export(dds, c("Condition", "FluD8", "Control"), "FluD8vsControl")
shrink_and_export(dds, c("Condition", "FluD8", "FluD4"), "FluD8vsFluD4")

saveRDS(dds, "C:/Users/GERSHOM/Desktop/data/rnaseq/applied-genomics/data/dds_object.rds")
cat("Analysis complete! \n")




# GENE ONTOLOGY -----------------------------------------------------------
library(biomaRt)
mouse_mart <- useMart(biomart = "ensembl",
                      host = "ensembl.org",
                      dataset = "mmusculus_gene_ensembl")

#Read TSV file for differentially expressed genes in FluD4vsControl file...
DEGs_FluD4 <- read.delim("./Results/dataFluD4vsControlSig_DEGs.tsv", header = TRUE, stringsAsFactors = FALSE)

# Extract ensembl IDs from its column in the FluD4vsControl sig. DEGs file
ensembl_ids4 <- DEGs_FluD4$gene_id 

#Read TSV file for differentially expressed genes in FluD8vsControl file...
DEGs_FluD8 <- read.delim("./Results/dataFluD8vsControlSig_DEGs.tsv", header = TRUE, stringsAsFactors = FALSE)

# Extract ensembl IDs from its column in the FluD8vsControl sig. DEGs file
ensembl_ids8 <- DEGs_FluD8$gene_id

#Read TSV file for differentially expressed genes in FluD4vsFluD8 file...
DEGs_FluD48 <- read.delim("./Results/dataFluD8vsFluD4Sig_DEGs.tsv", header = TRUE, stringsAsFactors = FALSE)

# Extract ensembl IDs from its column in the file
ensembl_ids48 <- DEGs_FluD48$gene_id 



#Ensemble query to get mouse gene official name for FluD4vsControl sig. DEGs file
result_FluD4 <- getBM(
  mart = mouse_mart,
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids4
)


#Ensemble query to get mouse gene official name for FluD8vsControl sig. DEGs file
result_FluD8 <- getBM(
  mart = mouse_mart,
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids8
)



#Ensemble query to get mouse gene official name for FluD4vsFluD8 sig. DEGs file
result_FluD48 <- getBM(
  mart = mouse_mart,
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids48
)


write_tsv(result_FluD8, file = "annotated_FluD8.tsv")
write_tsv(result_FluD4, file = "annotated_FluD4.tsv")
write_tsv(result_FluD48, file = "annotated_FluD48.tsv")

# GENE ONTOLOGY -----------------------------------------------------------
mouse_mart <- useMart(biomart = "ensembl",
                      host = "ensembl.org",
                      dataset = "mmusculus_gene_ensembl")
#Convert gene names to entrez ids for FluD4vs control 
annot_entrez4 <- getBM(values = result_FluD4$mgi_symbol,
                      mart = mouse_mart,
                      attributes = c("mgi_symbol",
                                     "entrezgene_id",
                                     "description"),
                      filters = "mgi_symbol")

#Convert gene names to entrez ids for FluD8 vs control 
annot_entrez8 <- getBM(values = result_FluD8$mgi_symbol,
                       mart = mouse_mart,
                       attributes = c("mgi_symbol",
                                      "entrezgene_id",
                                      "description"),
                       filters = "mgi_symbol")

#Convert gene names to entrez ids for FluD48 vs control 
annot_entrez48 <- getBM(values = result_FluD48$mgi_symbol,
                       mart = mouse_mart,
                       attributes = c("mgi_symbol",
                                      "entrezgene_id",
                                      "description"),
                       filters = "mgi_symbol")

#Convert entrez ids in the files to be character values instead of numerical values
annot_entrez4$entrezgene_id <- as.character(annot_entrez4$entrezgene_id)
annot_entrez8$entrezgene_id <- as.character(annot_entrez8$entrezgene_id)
annot_entrez48$entrezgene_id <- as.character(annot_entrez48$entrezgene_id)

library(clusterProfiler)
library(org.Mm.eg.db)
#Biological process analysis for the different analysis: FluD4; FluD8 & FluD48
BP_FluD4 <- enrichGO(gene = annot_entrez4$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #ontology is Biological Processing (BP)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
BP_FluD4_table <- as.data.frame(BP_FluD4)
view(BP_FluD4_table)

BP_FluD4_final <- clusterProfiler::simplify(BP_FluD4)
 #save file locally...
write_delim(
  x = as.data.frame(BP_FluD4@result),
  file = "FluD4_Biological process.csv",
  delim = ","
)



#FluD8 vs control
BP_FluD8 <- enrichGO(gene = annot_entrez8$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #ontology is Biological Processing (BP)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
BP_FluD8_table <- as.data.frame(BP_FluD8)
view(BP_FluD8_table)

BP_FluD8_final <- clusterProfiler::simplify(BP_FluD8)
#save file locally...
write_delim(
  x = as.data.frame(BP_FluD8@result),
  file = "FluD8_Biological process.csv",
  delim = ","
)


#FluD4 vs FluD8
BP_FluD48 <- enrichGO(gene = annot_entrez48$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #ontology is Biological Processing (BP)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
BP_FluD48_table <- as.data.frame(BP_FluD48)
view(BP_FluD48)

BP_FluD48_final <- clusterProfiler::simplify(BP_FluD48)
#save file locally...
write_delim(
  x = as.data.frame(BP_FluD48@result),
  file = "FluD48_Biological_process.csv",
  delim = ","
)
#This particular one yielded zero esults from the enrichGO() function


#dotplot
library(enrichplot)
#FluD4
png("BP_FluD4.png", width = 460, height = 550)
dotplot(BP_FluD4_final, showCategory = 10)
dev.off()

#FluD8
png("BP_FluD8.png", width = 460, height = 550)
dotplot(BP_FluD8_final, showCategory = 10)
dev.off()

#FluD48
png("BP_FluD48.png", width = 460, height = 550)
dotplot(BP_FluD48_final, showCategory = 10)
dev.off()


#Bar plot
#FluD4
png("BP_FluD4_barplot.png", width = 460, height = 550)
barplot(BP_FluD4_final, showCategory = 10)
dev.off()

#FluD8
png("BP_FluD8_barplot.png", width = 460, height = 550)
barplot(BP_FluD8_final, showCategory = 10)
dev.off()



# Cellular Component ------------------------------------------------------
CC_FluD4 <- enrichGO(gene = annot_entrez4$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC", #ontology is Cellular Component(CC)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
CC_FluD4_table <- as.data.frame(CC_FluD4)
view(CC_FluD4_table)

CC_FluD4_final <- clusterProfiler::simplify(CC_FluD4)
#save file locally...
write_delim(
  x = as.data.frame(CC_FluD4@result),
  file = "FluD4_CellularComponent.csv",
  delim = ","
)



#FluD8 vs control
CC_FluD8 <- enrichGO(gene = annot_entrez8$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC", #ontology is Cellular Component(CC)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
CC_FluD8_table <- as.data.frame(CC_FluD8)
view(CC_FluD8_table)

CC_FluD8_final <- clusterProfiler::simplify(CC_FluD8)
#save file locally...
write_delim(
  x = as.data.frame(CC_FluD8@result),
  file = "FluD8_CellularComponent.csv",
  delim = ","
)


#FluD4 vs FluD8
CC_FluD48 <- enrichGO(gene = annot_entrez48$entrezgene_id,
                      OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = "CC", #ontology is Cellular component
                      pAdjustMethod = "BH", #Benjamin H...method
                      qvalueCutoff = 0.05,
                      readable = FALSE,
                      pool = FALSE)
CC_FluD48_table <- as.data.frame(CC_FluD48)
view(CC_FluD48)

CC_FluD48_final <- clusterProfiler::simplify(CC_FluD48)
#save file locally...
write_delim(
  x = as.data.frame(CC_FluD48@result),
  file = "FluD48_CellularComponent.csv",
  delim = ","
)

#dotplot
library(enrichplot)
#FluD4
png("CC_FluD4.png", width = 460, height = 550)
dotplot(CC_FluD4_final, showCategory = 10)
dev.off()

#FluD8
png("CC_FluD8.png", width = 460, height = 550)
dotplot(CC_FluD8_final, showCategory = 10)
dev.off()

#FluD48
png("CC_FluD48.png", width = 460, height = 550)
dotplot(CC_FluD48_final, showCategory = 10)
dev.off()


#Bar plot
#FluD4
png("CC_FluD4_barplot.png", width = 460, height = 550)
barplot(CC_FluD4_final, showCategory = 10)
dev.off()

#FluD8
png("CC_FluD8_barplot.png", width = 460, height = 550)
barplot(CC_FluD8_final, showCategory = 10)
dev.off()

#FluD48
png("CC_FluD48_barplot.png", width = 460, height = 550)
barplot(CC_FluD48_final, showCategory = 10)
dev.off()



# Molecular Function ------------------------------------------------------
MF_FluD4 <- enrichGO(gene = annot_entrez4$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF", #ontology is Molecular Function(MF)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
MF_FluD4_table <- as.data.frame(MF_FluD4)
view(MF_FluD4_table)

MF_FluD4_final <- clusterProfiler::simplify(MF_FluD4)
#save file locally...
write_delim(
  x = as.data.frame(MF_FluD4@result),
  file = "FluD4_MolecularFunction.csv",
  delim = ","
)



#FluD8 vs control
MF_FluD8 <- enrichGO(gene = annot_entrez8$entrezgene_id,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF", #ontology is Molecular Function(MF)
                     pAdjustMethod = "BH", #Benjamin H...method
                     qvalueCutoff = 0.05,
                     readable = FALSE,
                     pool = FALSE)
MF_FluD8_table <- as.data.frame(MF_FluD8)
view(MF_FluD8_table)

MF_FluD8_final <- clusterProfiler::simplify(MF_FluD8)
#save file locally...
write_delim(
  x = as.data.frame(MF_FluD8@result),
  file = "FluD8_MolecularFunction.csv",
  delim = ","
)


#FluD4 vs FluD8
MF_FluD48 <- enrichGO(gene = annot_entrez48$entrezgene_id,
                      OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = "MF", #ontology is Molecular Function
                      pAdjustMethod = "BH", #Benjamin H...method
                      qvalueCutoff = 0.05,
                      readable = FALSE,
                      pool = FALSE)
MF_FluD48_table <- as.data.frame(MF_FluD48)
view(MF_FluD48)
#No data on Molecular Function betweeen FluD4 and FluD8

MF_FluD48_final <- clusterProfiler::simplify(MF_FluD48)
#save file locally...
write_delim(
  x = as.data.frame(CC_FluD48@result),
  file = "FluD48_MolecularFunction.csv",
  delim = ","
)

#dotplot
library(enrichplot)
#FluD4
png("MF_FluD4.png", width = 460, height = 550)
dotplot(MF_FluD4_final, showCategory = 10)
dev.off()

#FluD8
png("MF_FluD8.png", width = 460, height = 550)
dotplot(MF_FluD8_final, showCategory = 10)
dev.off()

#FluD48
png("MF_FluD48.png", width = 460, height = 550)
dotplot(MF_FluD48_final, showCategory = 10)
dev.off()


#Bar plot
#FluD4
png("MF_FluD4_barplot.png", width = 460, height = 550)
barplot(MF_FluD4_final, showCategory = 10)
dev.off()

#FluD8
png("MF_FluD8_barplot.png", width = 460, height = 550)
barplot(MF_FluD8_final, showCategory = 10)
dev.off()

#FluD48
png("MF_FluD48_barplot.png", width = 460, height = 550)
barplot(MF_FluD48_final, showCategory = 10)
dev.off()




# KEGG Analysis -----------------------------------------------------------
#FluD4 vs Control
KEGG_FluD4 <- enrichKEGG(gene = annot_entrez4$entrezgene_id,
                         organism = "mmu",
                         pAdjustMethod = "BH", #Benjamin H...method
                         qvalueCutoff = 0.05)
KEGG_FluD4_table <- as.data.frame(KEGG_FluD4)
view(KEGG_FluD4_table)
write.csv(KEGG_FluD4_table, file = "KEGGpathway_FluD4.csv")


#FluD8 vs control
KEGG_FluD8 <- enrichKEGG(gene = annot_entrez8$entrezgene_id,
                         organism = "mmu",
                         pAdjustMethod = "BH", #Benjamin H...method
                         qvalueCutoff = 0.05)
KEGG_FluD8_table <- as.data.frame(KEGG_FluD8)
view(KEGG_FluD8_table)
write.csv(KEGG_FluD8_table, file = "KEGGpathway_FluD8.csv")



#FluD4 vs FluD8
KEGG_FluD48 <- enrichKEGG(gene = annot_entrez48$entrezgene_id,
                         organism = "mmu",
                         pAdjustMethod = "BH", #Benjamin H...method
                         qvalueCutoff = 0.05)
KEGG_FluD48_table <- as.data.frame(KEGG_FluD48)
view(KEGG_FluD48_table)
#Yielded no results too!
write.csv(KEGG_FluD48_table, file = "KEGGpathway_FluD48.csv")



#dotplot
#library(enrichplot)
#FluD4
png("KEGG_FluD4.png", width = 460, height = 550)
dotplot(KEGG_FluD4, showCategory = 10)
dev.off()

#FluD8
png("KEGG_FluD8.png", width = 460, height = 550)
dotplot(KEGG_FluD8, showCategory = 10)
dev.off()

#FluD48
png("KEGG_FluD48.png", width = 460, height = 550)
dotplot(KEGG_FluD48, showCategory = 10)
dev.off()


#Bar plot
#FluD4
png("KEGG_FluD4_barplot.png", width = 460, height = 550)
barplot(KEGG_FluD4, showCategory = 10)
dev.off()

#FluD8
png("KEGG_FluD8_barplot.png", width = 460, height = 550)
barplot(KEGG_FluD8, showCategory = 10)
dev.off()

#FluD48
png("KEGG_FluD48_barplot.png", width = 460, height = 550)
barplot(KEGG_FluD48, showCategory = 10)
dev.off()

