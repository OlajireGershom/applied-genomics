if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ashr")

BiocManager::install("EnhancedVolcano")

BiocManager::install("clusterProfiler")

BiocManager::install("org.Mm.eg.db")

library(DESeq2)
library(ashr)
library(tidyverse)

files <- list.files(path = "counts/",
                    pattern = "\\.tabular$",
                    full.names = TRUE
                    )
print(files)

#Read eac count file
count_list <- lapply (files, function(f) {
  read.delim(f, header = TRUE)
})

#Merge all count tables by Geneid
count_matrix <- Reduce(function(x,y) merge(x, y, by = "Geneid", all = TRUE),
                       count_list)

#set gene IDs as row names
rownames(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL #remove te redundant Geneid column

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
