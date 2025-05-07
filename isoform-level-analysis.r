library(GEOquery)

# Data Preprocessing ####

if (file.exists("data/GEO/GSE261172-columns.csv")) {
	coldata <- read.csv("data/GEO/GSE261172-columns.csv")
} else {
	gds <- getGEO("GSE261172", destdir="data/GEO")

	gse <- gds[[1L]]
	coldata <- pData(gse)
	write.csv(coldata,"data/GEO/GSE261172-columns.csv")
}



cts=as.matrix(read.csv("results/transcript_count_matrix_renamed.csv", row.names="transcript_id"))

coldata$sample_number <- gsub("[A-Za-z]","", coldata$title)

row.names(coldata) <- coldata$sample_number

coldata$BioSample <- gsub(".*biosample/", "", coldata$relation)

label_map <- rownames(coldata)
names(label_map) <- coldata$BioSample
count_samples <- unname(label_map[colnames(cts)])

colnames(cts) <- count_samples

coldata$condition <- gsub("Sarcomatoid|[0-9]{2}", "", coldata$title)

keepSamples <- coldata$condition %in% c("High", "Low")
coldata <- coldata[keepSamples,]
cts <- cts[,coldata$sample_number]


coldata$condition <- factor(coldata$condition, levels=c("Low", "High"))

all(rownames(coldata) == colnames(cts))

rs <- rowSums(cts)
keep <- rs >= 10 & !is.na(rs)
cts <- cts[keep,]

# Performing DESeq ####

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds, BPPARAM = BiocParallel::MulticoreParam())
res <- results(dds, alpha = 0.05)
summary(res)
# res <- res[!is.na(res$padj),]

# General Analysis & Visualization ####
resOrdered <- res[order(res$pvalue),]
head(resOrdered, 10)

res05 <- results(dds, alpha=0.05)
res05 <- res05[res05$padj < 0.05 & !is.na(res05$padj),]
write.csv(res05, "results/top-transcripts-results.csv")

summary(res05)

## MA and Count Plots ####
plotMA(res, ylim=c(-3,3))

resLFC <- lfcShrink(dds, coef="condition_High_vs_Low", type="apeglm")
resLFC

plotMA(resLFC, ylim=c(-3,3))

plotCounts(dds, gene=which(res$padj == min(res$padj, na.rm=TRUE)), intgroup="condition")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library(ggplot2)
ggplot(d, aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks=c(25,100,400))

## Volcano Plots ####

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

library(ggrepel)
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

resdf <- as.data.frame(res)
resdf$Significance <- resdf$padj < 0.05 & abs(resdf$log2FoldChange) > 1 & !is.na(resdf$padj)
sigIdx <- resdf$Significance
resdf$Significance[sigIdx] <- "Significant"
resdf$Significance[!sigIdx] <- "Not Significant"
write.csv(resdf, "isoform-de-results.csv")

g <- ggplot(resdf,aes(x=log2FoldChange, y=-log10(padj), color=Significance)) + geom_point()


g + scale_color_manual(values=c("gray", "red")) +
    geom_hline(yintercept=-log10(0.05),linetype="dashed") +
    geom_vline(xintercept = -1,linetype="dashed") +
    geom_vline(xintercept = 1,linetype="dashed") +
    theme_minimal()

# Looking at specific Genes ####
pdl1 <- c("ENST00000381577.4","ENST00000381573.8","ENST00000498261.1","ENST00000492923.1","ENST00000474218.1")
res[pdl1,]

fosl1 <- c("ENST00000312562.7","ENST00000448083.6","ENST00000531493.5","ENST00000532401.1","ENST00000534222.1")
res[fosl1,]

ctla4 <- c("ENST00000696479.1","ENST00000696049.1","ENST00000650075.1","ENST00000648405.2","ENST00000487393.1","ENST00000427473.3","ENST00000295854.10")
res[ctla4,]

df <- head(as.data.frame(resOrdered),20)
ids <- data.frame(transcript_id=rownames(df))
write.table(df, 'results/top-DE-transcripts.txt',quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

df2 <- read.table("results/top-de-isoform-info.txt",header=TRUE,sep="\t", na.strings = TRUE)
names(df2) <- gsub(".", " ",names(df2),fixed=TRUE)
df2 <- df2[match(rownames(df), df2$`Transcript stable ID version`),]

df3 <- data.frame(transcript_id=rownames(df), gene_id=df2$`Gene stable ID`, gene_name=df2$`Gene name`, gene_type = df2$`Gene type`, canonical=!is.na(df2$`Ensembl Canonical`), lfc=df$log2FoldChange, padj=df$padj)
write.table(df3, 'results/top-DE-summary.txt',quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")



topid <- rownames(res05)

write.table(data.frame(topid), 'results/top005-transcript-ids.txt',quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

library(data.table)
dt <- fread("results/gene-and-transcripts.txt")
dt2 <- dt[, .N, by=`Gene stable ID`]
dt3 <- dt2[N > 1]
genes <- dt3$`Gene stable ID`
dt4 <- dt[`Gene stable ID` %in% genes, .(`Transcript stable ID version`, `Gene stable ID`)]

dt4$log2FoldChange <- res[dt4$`Transcript stable ID version`,"log2FoldChange"]
dt4$padj <- res[dt4$`Transcript stable ID version`,"padj"]
dt4
