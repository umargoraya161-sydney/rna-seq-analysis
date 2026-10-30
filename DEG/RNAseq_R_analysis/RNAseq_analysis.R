# RNA-seq differential expression analysis (course practical)
# Run from the DEG directory containing counts.csv and metadata2.csv.
# Install packages once if needed (not on every run):
# install.packages(c("BiocManager", "ggplot2", "dplyr"))
# BiocManager::install(c("DESeq2", "org.Hs.eg.db", "ComplexHeatmap"))

library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

# 1. Load counts and sample information
count <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("metadata2.csv", row.names = 1, check.names = FALSE)

# The course metadata may store sample IDs in a column called X.
if (!all(colnames(count) %in% rownames(meta)) && "X" %in% names(meta)) {
  rownames(meta) <- meta$X
}
if (!all(colnames(count) %in% rownames(meta))) {
  stop("Count matrix sample names do not match metadata. Check the CSV headers.")
}
meta <- meta[colnames(count), , drop = FALSE]
stopifnot(identical(colnames(count), rownames(meta)))
stopifnot(all(vapply(count, is.numeric, logical(1))))
stopifnot("dexamethasone" %in% names(meta))

# 2. Differential expression using DESeq2
meta$dexamethasone <- factor(meta$dexamethasone)
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta,
                              design = ~ dexamethasone)
# Remove genes with fewer than 10 total reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
deg <- DESeq(dds)
res <- results(deg)
res0.05 <- results(deg, alpha = 0.05)
print(summary(res0.05))
write.csv(as.data.frame(res), "test_udemy2.csv")

# 3. Add human gene symbols
res0.05.df <- as.data.frame(res0.05)
res0.05.df$Symbol <- mapIds(org.Hs.eg.db,
                            keys = rownames(res0.05.df),
                            keytype = "ENSEMBL", column = "SYMBOL",
                            multiVals = "first")
write.csv(res0.05.df, "Finaltest.csv")

# 4. PCA plot
vsd <- vst(deg, blind = FALSE)
png("pcaplot.png", width = 1000, height = 750, res = 120)
print(plotPCA(vsd, intgroup = "dexamethasone"))
dev.off()

# 5. Inspect size factors and dispersion
print(sizeFactors(deg))
png("Dispersionplot.png", width = 1000, height = 750, res = 120)
plotDispEsts(deg)
dev.off()

# 6. MA plot
png("MAplot.png", width = 1000, height = 750, res = 120)
plotMA(res0.05)
dev.off()

# 7. Ten genes with smallest adjusted p-values
ranked <- res0.05.df %>% filter(!is.na(padj)) %>% arrange(padj)
best_genes <- head(ranked, 10)
best_genes$GeneID <- rownames(best_genes)
write.csv(best_genes, "best_genes.csv", row.names = FALSE)

# 8. Volcano plot: FDR < 0.05 and |log2 fold change| > 1
vol <- res0.05.df %>% filter(!is.na(padj), padj > 0)
vol$Significant <- with(vol, padj < 0.05 & abs(log2FoldChange) > 1)
vol$GeneID <- rownames(vol)
p <- ggplot(vol, aes(x = log2FoldChange, y = -log10(padj),
                     color = Significant)) +
  geom_point(alpha = 0.65, size = 1) +
  geom_text(data = subset(best_genes, !is.na(Symbol) & padj > 0),
            aes(label = Symbol), hjust = -0.1, check_overlap = TRUE,
            na.rm = TRUE, show.legend = FALSE) +
  labs(x = "log2 fold change", y = "-log10 adjusted p-value",
       title = "Volcano plot") +
  theme_minimal()
ggsave("Volcanoplot.png", plot = p, width = 9, height = 6, dpi = 150)

# 9. Heatmap of 30 genes with smallest adjusted p-values
# Use stable Ensembl IDs for selecting rows, not gene symbols.
top_genes <- head(ranked, 30)
mat <- counts(deg, normalized = TRUE)[rownames(top_genes), , drop = FALSE]
mat.z <- t(scale(t(mat)))
mat.z[!is.finite(mat.z)] <- 0
row_labels <- ifelse(is.na(top_genes$Symbol) | top_genes$Symbol == "",
                     rownames(top_genes), top_genes$Symbol)
png("heatmapplot.png", width = 1200, height = 1000, res = 120)
ComplexHeatmap::draw(Heatmap(mat.z, name = "Row Z-score",
                             cluster_rows = TRUE, cluster_columns = TRUE,
                             row_labels = row_labels,
                             column_labels = colnames(mat.z)))
dev.off()

message("Analysis complete. Review plots and CSV files in the DEG directory.")
