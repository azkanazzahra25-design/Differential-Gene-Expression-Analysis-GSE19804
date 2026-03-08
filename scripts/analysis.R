# Differential Gene Expression Analysis
# Dataset: GSE19804

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GEOquery","limma"))
BiocManager::install("hgu133plus2.db")

install.packages(c("pheatmap","ggplot2","dplyr"))

if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
library(GEOquery)
library(limma)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
library(umap)

gset <- getGEO("GSE19804", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

library(GEOquery)
library(limma)
gset <- getGEO("GSE19804", GSEMatrix = TRUE)
gset <- gset[[1]]

expr <- exprs(gset)
pheno <- pData(gset)

head(expr)
head(pheno)
#library() digunakan agar fungsi di dalam package bisa digunakan
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

colnames(pheno)
head(pheno$source_name_ch1)
unique(pheno$source_name_ch1)

unique(pheno$source_name_ch1)

group <- factor(pheno$source_name_ch1)

table(group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(
  Tumor_vs_Normal = `frozen tissue of primary tumor` - `frozen tissue of adjacent normal`,
  levels = design
)

contrast.matrix
group <- factor(ifelse(pheno$source_name_ch1 == "frozen tissue of primary tumor",
                       "Tumor","Normal"))

table(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

contrast.matrix <- makeContrasts(
  Tumor - Normal,
  levels = design
)

contrast.matrix

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="fdr", number=Inf)

head(results)

install.packages("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(
  results,
  lab = rownames(results),
  x = 'logFC',
  y = 'P.Value',
  title = 'Differential Gene Expression: Tumor vs Normal',
  pCutoff = 0.05,
  FCcutoff = 1
)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
n
library(EnhancedVolcano)
EnhancedVolcano(
  results
  lab = rownames(results)
  x = 'logFC',
  y = 'P.Value',
  title = 'Differential Gene Expression: Tumor vs Normal',
  pCutoff = 0.05,
  FCcutoff = 1
)
EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "P.Value",
  title = "Differential Gene Expression: Tumor vs Normal",
  pCutoff = 0.05,
  FCcutoff = 1
)
top_genes <- rownames(results)[1:50]
heatmap_data <- expr[top_genes,]

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE
)

write.csv(results, "DEG_results_GSE19804.csv")
