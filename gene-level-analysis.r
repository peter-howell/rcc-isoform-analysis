###########
#imports
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install("GEOquery")
BiocManager::install("biomaRt")
BiocManager::install("reshape2")
BiocManager::install("ggplot2")
BiocManager::install("pasilla")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("DEGreport")
install.packages("pheatmap")
install.packages("tidyverse")

BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

library("GEOquery")
library(ggplot2)
library("biomaRt")
library(reshape2)
library(dplyr)

library("pasilla")
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")

library("DEGreport")
library(ggrepel)

# Get data
gds <- getGEO("GSE261172", GSEMatrix = TRUE)
gse <- gds[[1]]
coldata <- pData(gse)

cts <- as.matrix(read.csv("C:\\Users\\User\\OneDrive - Worcester Polytechnic Institute (wpi.edu)\\WPI\\2024-2025\\C_D term\\BCB 590\\BCB 590 Project\\GSE261172_RSEM_counts.csv",
                            row.names = 1))

# Remove the "X" prefix (e.g., "X20" -> "20")
colnames(cts) <- sub("^X", "", colnames(cts))

# Extract numeric sample number from the title (e.g., "SarcomatoidLow29" -> 29)
coldata$sample_number <- as.character(str_extract(coldata$title, "\\d+$"))

# Optional: Confirm this matches columns in the counts file
# You'll need to make sure the column names in cts are also characters:
colnames(cts) <- as.character(colnames(cts))
# Set row names to be sample id
rownames(coldata) <- coldata$sample_number


# Check if the row names in coldata match the column names in cts
all(rownames(coldata) %in% colnames(cts))  # Should return TRUE

# Drop sample 27 from both datasets since its HighLow
coldata <- coldata[coldata$sample_number != "27", ]
cts <- cts[, colnames(cts) != "27"]


# Check data
as_tibble(cts)
coldata
rownames(coldata)
colnames(cts)

# Ensure all values in the count matrix are rounded to integers
cts <- round(cts)


# Extract Sarcomatoid condition from the title (e.g., "SarcomatoidHigh33")
coldata$condition <- ifelse(str_detect(coldata$title, "High"), "High", "Low")

# Set as factor (Low as reference)
coldata$condition <- factor(coldata$condition, levels = c("Low", "High"))

# Check if column names in cts match row names in coldata
setdiff(colnames(cts), rownames(coldata))  # Sample names in cts not in coldata
setdiff(rownames(coldata), colnames(cts))  # Sample names in coldata not in cts

# Reorder coldata to match the order of columns in cts
coldata <- coldata[match(colnames(cts), rownames(coldata)), ]

# Check if order of rows in coldata matches columns in cts
all(rownames(coldata) == colnames(cts))  # Should return TRUE

# Output rownames of coldata and colnames of cts- should be the same
rownames(coldata)
colnames(cts)

#remaining code for the DEseq2 tutorial
# Construct a DESeqDataSet:
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Factor levels
dds$condition <- factor(dds$condition, levels = c("Low","High"))

# Genes differentially expressed between conditions
dds <- DESeq(dds)

res <- results(dds)
res

res <- results(dds, name="condition_High_vs_Low")

#Log fold change shrinkage for visualization and ranking
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_High_vs_Low", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]
summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


#37 differentially expressed genes

plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Prepare volcano plot data (remove NA padj)
volcano_df <- as.data.frame(resLFC) %>%
  rownames_to_column(var = "Gene") %>%
  filter(!is.na(padj)) %>%
  mutate(Significance = case_when(
    padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
    TRUE ~ "Not Significant"
  ))

# Select all genes that are significant
labeled_genes <- volcano_df %>%
  filter(Significance == "Significant")

# Volcano plot with labels
ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = labeled_genes,
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf  # Allow labeling all
  )



mcols(res)$description


resSig <- subset(resOrdered, padj < 0.1)
resSig

vsd <- vst(dds)
head(assay(vsd), 3)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])

# Gene heatmap
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Heatmap of the sample-to-sample distances
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# PCA
plotPCA(vsd, intgroup=c("condition", "type"))

degs <- degComps(dds, combs = "condition",
                 contrast = list( c("condition", "Low", "High") ))
names(degs)
# The first result set which is -> condition_tr_vs_nt
deg(degs[[1]])
deg(degs[[1]], "raw", "tibble")

significants(degs, fc = 0, fdr = 0.05, full = TRUE)

# Gene plot
degPlot(dds = dds, res = res, n = 6, xs = "condition")

rownames_to_column(as.data.frame(res), var ="Gene") %>% as_tibble()

BiocManager::install("ReportingTools")
library("ReportingTools")
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./reports")

## Get a DESeqDataSet object
mockRna.dse <- DESeq(dds)
colData(mockRna.dse)$conditions <- dds$condition

# this will take a while
publish(mockRna.dse,des2Report, pvalueCutoff=0.05,
        annotation.db="fly.db0", factor = colData(mockRna.dse)$conditions,
        reportDir="./reports")
finish(des2Report)


############
#GO analysis
# Create list of significant gene names (padj < 0.05 and abs(log2FC) > 1)
sig_genes <- as.data.frame(resLFC) %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)


# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Preview
head(gene_entrez)
ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",       # Change to "MF" or "CC" if needed
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

dotplot(ego, showCategory = 40) +
  ggtitle("GO Biological Processes Enriched in DEGs")

# OR barplot
barplot(ego, showCategory = 40)


