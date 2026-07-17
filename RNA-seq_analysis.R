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
