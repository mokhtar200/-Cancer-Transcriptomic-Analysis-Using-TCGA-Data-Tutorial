# TCGA Cancer Transcriptomic Analysis 

# Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "clusterProfiler", "org.Hs.eg.db", "pheatmap", "ggplot2"))

library(TCGAbiolinks)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)

#---------------------------#
# Step 1: Data Acquisition  #
#---------------------------#

# Query TCGA-BRCA Data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

# Download and Prepare Data
GDCdownload(query)
data <- GDCprepare(query)

#---------------------------#
# Step 2: Data Preprocessing#
#---------------------------#

# Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = assay(data),
                              colData = colData(data),
                              design = ~ sample_type)

# Filter Low-Count Genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Differential Expression Analysis
dds <- DESeq(dds)

#---------------------------#
# Step 3: DEG Analysis       #
#---------------------------#

# Extract DEGs
res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
res <- res[order(res$pvalue), ]

# Save DEGs
write.csv(as.data.frame(res), "DEGs.csv")

#---------------------------#
# Step 4: Enrichment Analysis#
#---------------------------#

# Select Significant DEGs
sig_genes <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])

# Perform GO Enrichment
ego <- enrichGO(gene = sig_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Save Enrichment Results
write.csv(as.data.frame(ego), "GO_enrichment.csv")

#---------------------------#
# Step 5: Visualization     #
#---------------------------#

# Volcano Plot
ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05)) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  ggsave("volcano_plot.png")

# Heatmap of Top 50 DEGs
top_genes <- head(rownames(res[order(res$padj), ]), 50)
pheatmap(assay(dds)[top_genes, ], scale = "row", clustering_distance_rows = "euclidean")

