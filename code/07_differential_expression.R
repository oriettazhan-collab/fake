##############################
# 第1步：安装必须的包（仅第一次需要）
##############################

if (!require("BiocManager")) install.packages("BiocManager")

BiocManager::install("DESeq2", update = FALSE)

install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("stringr")

##############################

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

setwd("C:/Users/orietta.zhan/Desktop/G A")


count_data <- read.table(
  "final_count_matrix_for_R.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

############################## SAMPLE INFO


sample_info <- data.frame(
  condition = factor(c(
    "BHI",
    "BHI",
    "BHI",
    "Serum",
    "Serum",
    "Serum"
  ))
)

rownames(sample_info) <- colnames(count_data)

##############################DESeq2 分析

dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData = sample_info,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Serum", "BHI")
)

##############################RESULT

write.table(
  res,
  "DESeq2_all_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

##############################PCA PLOT


pca_data <- rlog(dds)

pdf(
  "PCA_plot.pdf",
  width = 7,
  height = 5
)

plotPCA(
  pca_data,
  intgroup = "condition"
) +
  theme_bw() +
  ggtitle("PCA: BHI vs Serum")

dev.off()

##############################READ GFF 注释


gff <- readLines("E745_clean.gff")


gff <- gff[!grepl("^#", gff)]

#  CDS ONLY
gff_cds <- gff[grepl("\tCDS\t", gff)]
# WITHDRAW gene ID
extract_id <- function(x) {
  
  id <- str_match(
    x,
    "ID=([^;]+)"
  )[,2]
  
  return(id)
}

# gene symbol

extract_gene <- function(x) {
  
  gene <- str_match(
    x,
    "gene=([^;]+)"
  )[,2]
  
  return(gene)
}

extract_product <- function(x) {
  
  product <- str_match(
    x,
    "product=([^;]+)"
  )[,2]
  
  return(product)
}

# BUILD annotation 

annotation_df <- data.frame(
  gene_id = sapply(gff_cds, extract_id),
  gene_name = sapply(gff_cds, extract_gene),
  product = sapply(gff_cds, extract_product),
  stringsAsFactors = FALSE
)

annotation_df <- annotation_df %>%
  distinct(gene_id, .keep_all = TRUE)


res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

# CONBINE annotation
res_df <- merge(
  res_df,
  annotation_df,
  by = "gene_id",
  all.x = TRUE
)

#  locus tag

res_df$label <- ifelse(
  is.na(res_df$gene_name) |
    res_df$gene_name == "",
  res_df$gene_id,
  res_df$gene_name
)

# REMOVE NA


res_df <- na.omit(res_df)
res_df$status <- "Not significant"

res_df$status[
  res_df$log2FoldChange > 1 &
    res_df$padj < 0.05
] <- "Upregulated"

res_df$status[
  res_df$log2FoldChange < -1 &
    res_df$padj < 0.05
] <- "Downregulated"



# TWO SIDE
left_genes <- res_df %>%
  filter(
    log2FoldChange <= -3,
    padj < 0.05
  ) %>%
  arrange(log2FoldChange) %>%
  head(20)

right_genes <- res_df %>%
  filter(
    log2FoldChange >= 5,
    padj < 0.05
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

label_genes <- rbind(
  left_genes,
  right_genes
)

write.table(
  label_genes,
  "Extreme_volcano_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# VOLCANO
pdf(
  "Volcano_plot_publication.pdf",
  width = 10,
  height = 8
)

ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = status
  )
) +

geom_point(
  size = 2,
  alpha = 0.7
) +

scale_color_manual(
  values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "gray70"
  )
) +
  

geom_vline(
  xintercept = c(-1, 1),
  linetype = "dashed",
  color = "black"
) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  

geom_text_repel(
  data = label_genes,
  aes(label = label),
  
  size = 4,
  
  box.padding = 0.5,
  point.padding = 0.4,
  
  segment.color = "black",
  
  max.overlaps = 100
) +

ggtitle(
  "Volcano Plot: Serum vs BHI"
) +
  

xlab("log2(Fold Change)") +
  
  ylab("-log10(Adjusted p-value)") +

theme_bw(base_size = 15) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    legend.title = element_blank(),
    
    panel.grid.minor = element_blank()
  )

dev.off()

cat("====================================\n")
cat("✅ FINISH\n")
cat("====================================\n")
cat("OUTPUT：\n")
cat("1. PCA_plot.pdf\n")
cat("2. Volcano_plot_publication.pdf\n")
cat("3. DESeq2_all_genes.txt\n")
cat("4. Extreme_volcano_genes.txt\n")
cat("====================================\n")
