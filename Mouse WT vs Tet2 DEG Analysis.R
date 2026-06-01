# ==============================================================================
# DEG & Functional Enrichment Analysis: Mouse WT vs Tet2-KO (LSK Cells)
# WT  = Wild-Type (Reference Condition)
# KO  = Tet2 Knockout bone marrow stem cells
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES & INITIALIZATION
# ------------------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(forcats)
library(DOSE)


# Resolve namespace conflicts — must come AFTER all library() calls
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate


# Define global results directory
output_folder <- "DEG_Analysis_Results"
if (!dir.exists(output_folder)) dir.create(output_folder)

# ------------------------------------------------------------------------------
# 2. DATA LOADING, CLEANING, & EXPRESSION MATRIX PREPARATION
# ------------------------------------------------------------------------------

raw_data <- read_excel("Mouse WT vs Tet2 - LSK data.xlsx", sheet = "DE_analysis")

cat("BEFORE:", dim(raw_data), "\n")

cleaned_counts <- raw_data %>%
  dplyr::select(GeneSymbol, starts_with("WT"), starts_with("Tet2")) %>%
  dplyr::filter(!is.na(GeneSymbol)) %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()

cat("AFTER :", dim(cleaned_counts), "\n")

counts_matrix <- cleaned_counts %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

print(dim(counts_matrix))
print(head(counts_matrix))


# ------------------------------------------------------------------------------
# 3. EXPERIMENTAL INFO
# ------------------------------------------------------------------------------

sample_info <- data.frame(
  row.names = colnames(counts_matrix),
  condition = factor(
    c(rep("WT", 3), rep("Tet2", 4)),
    levels = c("WT", "Tet2")     # WT explicitly set as reference
  )
)

# Positive LFC = upregulated in Tet2 (relative to WT)
# Negative LFC = downregulated in Tet2 (relative to WT)
print(sample_info)


# ------------------------------------------------------------------------------
# 4. DIFFERENTIAL EXPRESSION ANALYSIS, LFC SHRINKAGE & DESeq2 OBJECT SETUP
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = sample_info,
  design    = ~ condition
)

# Keep genes with >= 10 counts in at least 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]

cat("\nGenes retained after filtering:", nrow(dds), "\n")


# RUN DESeq2 
dds <- DESeq(dds)

resultsNames(dds)   # confirm: should show "condition_Tet2_vs_WT"

# Raw results
res_raw <- results(dds, name = "condition_Tet2_vs_WT")

# Shrinkage — apeglm (best practice)
res <- lfcShrink(dds, coef = "condition_Tet2_vs_WT", type = "apeglm")

# Order by adjusted p-value — filter NAs
res_ordered <- res[order(res$padj), ]
res_ordered <- res_ordered[!is.na(res_ordered$padj) & res_ordered$padj < 0.05, ]

cat("\nTop 10 DE genes:\n")
print(head(res_ordered, 10))

# Save full results
write.csv(as.data.frame(res_ordered),
          file = file.path(output_folder, "DEG_results_Tet2_vs_WT.csv"),
          row.names = TRUE)
cat("DEG results saved.\n")


# ------------------------------------------------------------------------------
# 5.VISUALIZATIONS: MA & PCA PLOTS
# ------------------------------------------------------------------------------

graphics.off()

png(filename = file.path(output_folder, "MA_Plot.png"),
    width = 1400, height = 650, res = 140)

par(mfrow = c(1, 2))

plotMA(results(dds, name = "condition_Tet2_vs_WT"),
       ylim = c(-5, 5),
       main = "Before Shrinkage (Raw LFC)",
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change")
abline(h = 0, col = "red", lwd = 2)

plotMA(res,
       ylim = c(-5, 5),
       main = "After Shrinkage (apeglm)",
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change")
abline(h = 0, col = "red", lwd = 2)

dev.off()
cat("MA plot saved.\n")


# PCA PLOT
vst_data <- vst(dds, blind = FALSE)
plotPCA(vst_data, intgroup = "condition")

# ------------------------------------------------------------------------------
# 6. CORE EXPRESSION VISUALIZATION: VOLCANO PLOTS
# ------------------------------------------------------------------------------

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  dplyr::mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      padj < 0.05                        ~ "Significant (|LFC|<1)",
      TRUE                               ~ "NS"
    )
  )

# DEG counts
n_up    <- sum(res_df$significance == "Upregulated")
n_down  <- sum(res_df$significance == "Downregulated")
n_total <- n_up + n_down

cat("Upregulated  :", n_up,    "\n")
cat("Downregulated:", n_down,  "\n")
cat("Total DEGs   :", n_total, "\n")

# Target genes
genes_to_highlight <- c(
  "Il17ra", "Glp1r",  "Sik2",   "Sik3",  "Camp",   "Nfe2l2",
  "Creb1",  "Crtc2",  "Adcy1",  "Adcy2", "Atf4",   "Atf6",
  "Il6",    "Eif2s1", "Hspa5",  "Xbp1",  "Gclc",   "Nqo1"
)

res_df <- res_df %>%
  dplyr::mutate(
    label     = ifelse(gene %in% genes_to_highlight, gene, NA),
    highlight = gene %in% genes_to_highlight
  )

# Color palette — defined ONCE here, used in both plots
colors <- c(
  "Upregulated"           = "#C0392B",
  "Downregulated"         = "#2980B9",
  "Significant (|LFC|<1)" = "#E67E22",
  "NS"                    = "grey80"
)

# x_max — defined ONCE here, used in both plots
x_max <- max(abs(res_df$log2FoldChange), na.rm = TRUE)


#  function to build volcano
build_volcano <- function(df, y_max, y_annot) {
  
  df <- df %>% dplyr::mutate(y_plot = pmin(-log10(padj), y_max - 5))
  
  ggplot(df, aes(x = log2FoldChange, y = y_plot, color = significance)) +
    
    geom_point(size = 1.0, alpha = 0.5) +
    
    geom_point(data = dplyr::filter(df, highlight),
               size = 2.8, shape = 21, stroke = 0.7,
               aes(fill = significance), color = "black") +
    
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "black", linewidth = 0.4) +
    
    geom_label_repel(
      aes(label = label), na.rm = TRUE, size = 2.8, fontface = "bold",
      color = "black", fill = "white", box.padding = 0.4, point.padding = 0.3,
      segment.color = "black", segment.size = 0.3,
      max.overlaps = Inf, force = 6, min.segment.length = 0
    ) +
    
    annotate("text", x = -x_max * 0.75, y = y_annot,
             label = paste0("Down: ", n_down),
             color = "#2980B9", fontface = "bold", size = 4, hjust = 0) +
    annotate("text", x = 0, y = y_annot,
             label = paste0("Total DEGs: ", n_total),
             color = "black", fontface = "bold", size = 4, hjust = 0.5) +
    annotate("text", x = x_max * 0.75, y = y_annot,
             label = paste0("Up: ", n_up),
             color = "#C0392B", fontface = "bold", size = 4, hjust = 1) +
    
    scale_color_manual(values = colors,
                       breaks = c("Upregulated", "Downregulated",
                                  "Significant (|LFC|<1)", "NS")) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0.01, 0)) +
    
    labs(
      title    = "Tet2-KO vs WT: LSK Cells",
      subtitle = "Thresholds: padj < 0.05  |  |log\u2082FC| > 1",
      x        = expression(Log[2]~"Fold Change"),
      y        = expression(-Log[10]~"(adjusted p-value)"),
      color    = "Expression", fill = "Expression"
    ) +
    
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey40"),
      axis.title      = element_text(face = "bold", size = 10),
      axis.text       = element_text(size = 9),
      legend.title    = element_text(face = "bold", size = 9),
      legend.text     = element_text(size = 8),
      legend.position = "bottom",
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin     = margin(20, 20, 15, 15)
    ) +
    guides(fill = "none")
}

# Plot 1: Main 
p_main <- build_volcano(res_df, y_max = 310, y_annot = 300)
ggsave(file.path(output_folder, "Volcano_Main.png"),
       plot = p_main, width = 6, height = 6, dpi = 300)
cat("Main volcano saved.\n")

# Plot 2: Trimmed 
p_trimmed <- build_volcano(res_df, y_max = 100, y_annot = 93)
ggsave(file.path(output_folder, "Volcano_Trimmed.png"),
       plot = p_trimmed, width = 6, height = 6, dpi = 300)
cat("Trimmed volcano saved.\n")


# TARGET GENE CSV
target_gene_summary <- res_df %>%
  dplyr::filter(gene %in% genes_to_highlight) %>%
  dplyr::select(gene, log2FoldChange, pvalue, padj, significance) %>%
  dplyr::mutate(
    direction = case_when(
      significance == "Upregulated"           ~ "Up in Tet2",
      significance == "Downregulated"         ~ "Down in Tet2",
      significance == "Significant (|LFC|<1)" ~ "Significant (small LFC)",
      TRUE                                    ~ "Not Significant"
    ),
    log2FoldChange = round(log2FoldChange, 4),
    pvalue         = signif(pvalue, 4),
    padj           = signif(padj,   4)
  ) %>%
  arrange(significance, desc(abs(log2FoldChange)))

cat("\n===== TARGET GENE SUMMARY =====\n")
print(target_gene_summary, row.names = FALSE)

cat("\n===== BREAKDOWN =====\n")
print(table(target_gene_summary$significance))

write.csv(target_gene_summary,
          file = file.path(output_folder, "Target_Gene_Summary.csv"),
          row.names = FALSE)
cat("Target gene summary saved.\n")

# ------------------------------------------------------------------------------
# 7. TARGET GENE COMPILATION & CLUSTERED HEATMAP
# ------------------------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)

heatmap_genes <- genes_to_highlight[genes_to_highlight %in% rownames(norm_counts)]
cat("Genes in heatmap:", length(heatmap_genes), "\n")

heatmap_mat    <- log2(norm_counts[heatmap_genes, ] + 1)
heatmap_scaled <- t(scale(t(heatmap_mat)))

col_annotation <- data.frame(
  Condition = sample_info$condition,
  row.names = rownames(sample_info)
)

annotation_colors <- list(
  Condition = c(WT = "#4E79A7", Tet2 = "#F28E2B")
)

ggsave(filename = file.path(output_folder, "Heatmap_Target_Genes.png"),
    width = 900, height = 800, res = 130)

pheatmap(heatmap_scaled,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         annotation_col    = col_annotation,
         annotation_colors = annotation_colors,
         color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         border_color      = NA,
         fontsize_row      = 11,
         fontsize_col      = 10,
         main              = "Target Genes: Z-score Normalized Expression\nTet2-KO vs WT LSK Cells",
         show_colnames     = TRUE,
         angle_col         = 45)

dev.off()
cat("Heatmap saved.\n")


# ------------------------------------------------------------------------------
# 8. FUNCTIONAL ENRICHMENT: KEGG OVER-REPRESENTATION ANALYSIS (ORA)
# ------------------------------------------------------------------------------

# Significant DE genes for KEGG
sig_genes <- as.data.frame(res_ordered) %>%
  rownames_to_column("gene") %>%
  dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

cat("\nSignificant DE genes for KEGG:", nrow(sig_genes), "\n")

# Convert to Entrez IDs
entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys      = sig_genes$gene,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)
cat("Mapped to Entrez IDs:", length(entrez_ids), "\n")

# Run KEGG
kegg_result <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "mmu",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2
)

if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
  cat("No significant KEGG pathways found.\n")
} else {
  
  cat("\nTop 15 KEGG pathways:\n")
  print(head(kegg_result@result[, c("Description", "GeneRatio", "pvalue", "p.adjust", "Count")], 15))
  
  write.csv(kegg_result@result,
            file = file.path(output_folder, "KEGG_pathway_results.csv"),
            row.names = FALSE)
  
  # Standard dot plot
  png(filename = file.path(output_folder, "KEGG_dotplot.png"),
      width = 1100, height = 800, res = 130)
  print(
    dotplot(kegg_result, showCategory = 20, font.size = 10) +
      ggtitle("KEGG Pathway Enrichment\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()
  
  # Standard bar plot
  png(filename = file.path(output_folder, "KEGG_barplot.png"),
      width = 1100, height = 800, res = 130)
  print(
    barplot(kegg_result, showCategory = 20, font.size = 10) +
      ggtitle("KEGG Pathway Enrichment\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()
  
  cat("KEGG standard plots saved.\n")
}


# KEGG — 11 SELECTED PATHWAYS

pathways_11 <- c(
  "PI3K-Akt signaling pathway",
  "Hematopoietic cell lineage",
  "Calcium signaling pathway",
  "JAK-STAT signaling pathway",
  "Ras signaling pathway",
  "IL-17 signaling pathway",
  "Th17 cell differentiation",
  "cAMP signaling pathway",
  "MAPK signaling pathway",
  "Rap1 signaling pathway",
  "HIF-1 signaling pathway"
)

kegg_11 <- kegg_result@result %>%
  dplyr::filter(Description %in% pathways_11) %>%
  dplyr::mutate(
    GeneRatio_numeric = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    Significant = ifelse(p.adjust < 0.05, "Significant", "Not Significant"),
    Description = fct_reorder(Description, GeneRatio_numeric)
  )

cat("Pathways found:", nrow(kegg_11), "\n")
print(kegg_11[, c("Description", "p.adjust", "Count", "Significant")])

p_kegg_11 <- ggplot(kegg_11,
                    aes(x = GeneRatio_numeric,
                        y = Description,
                        size = Count,
                        color = p.adjust)) +
  
  geom_rect(
    data        = dplyr::filter(kegg_11, Significant == "Not Significant"),
    aes(ymin    = as.numeric(Description) - 0.45,
        ymax    = as.numeric(Description) + 0.45,
        xmin    = -Inf, xmax = Inf),
    fill        = "#F2F3F4",
    alpha       = 0.8,
    inherit.aes = FALSE
  ) +
  
  geom_point(alpha = 0.9) +
  
  geom_point(
    data   = dplyr::filter(kegg_11, Significant == "Significant"),
    shape  = 21, stroke = 1.5,
    color  = "#C0392B", fill = NA,
    aes(size = Count)
  ) +
  
  geom_point(
    data   = dplyr::filter(kegg_11, Significant == "Not Significant"),
    shape  = 21, stroke = 1.2,
    color  = "grey50", fill = NA,
    aes(size = Count)
  ) +
  
  geom_text(
    aes(label = ifelse(p.adjust < 0.001,
                       formatC(p.adjust, format = "e", digits = 1),
                       paste0("p=", round(p.adjust, 3)))),
    hjust = -0.15, size = 3, color = "grey30"
  ) +
  
  scale_y_discrete(
    labels = setNames(
      ifelse(levels(kegg_11$Description) %in%
               dplyr::filter(kegg_11, Significant == "Significant")$Description,
             paste0(levels(kegg_11$Description), "  *"),
             paste0(levels(kegg_11$Description), "  (ns)")),
      levels(kegg_11$Description)
    )
  ) +
  
  scale_color_gradient(
    low   = "#C0392B", high = "#AED6F1",
    name  = "Adjusted\np-value",
    guide = guide_colorbar(reverse = TRUE)
  ) +
  scale_size_continuous(name = "Gene Count", range = c(4, 13)) +
  
  labs(
    title    = "KEGG Pathway Enrichment \u2014 Selected Pathways",
    subtitle = "Tet2-KO vs WT LSK Cells  |  * Significant (padj < 0.05)  |  ns = not significant",
    x        = "Gene Ratio",
    y        = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
    axis.text.y   = element_text(size = 10, face = "bold"),
    axis.text.x   = element_text(size = 10),
    axis.title.x  = element_text(face = "bold", size = 11),
    legend.title  = element_text(face = "bold", size = 10),
    legend.text   = element_text(size = 9),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin   = margin(20, 100, 15, 15)
  )

ggsave(file.path(output_folder, "KEGG_11_Pathways_Final.png"),
       plot = p_kegg_11, width = 12, height = 7, dpi = 300)
cat("KEGG 11-pathway plot saved.\n")

# ==============================================================================
#                                 GSEA Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 9. ADVANCED ENRICHMENT: GENE ONTOLOGY (GO) GSEA SUITE
# ------------------------------------------------------------------------------
# Using DESeq2 dataframe  
# (Code Ref: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/)
df_gsea <- res_df

# We want the log2 fold change 
original_gene_list <- df_gsea$log2FoldChange

# Name the vector using your Gene Symbols (stored in the 'gene' column)
names(original_gene_list) <- df_gsea$gene

# Omit any NA values 
gene_list <- na.omit(original_gene_list)

# Sort the list in decreasing order (strictly required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Ranked gene list prepared! Total genes:", length(gene_list), "\n")

# RUN MOUSE GSEA GO 
gse <- gseGO(
  geneList      = gene_list,
  ont           = "ALL",             # Looks at all GO categories (BP, CC, MF)
  keyType       = "SYMBOL",          # Matches your mouse gene symbols (e.g., Tet2)
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  verbose       = TRUE,
  OrgDb         = org.Mm.eg.db,      # This explicitly tells R to use MOUSE data
  pAdjustMethod = "BH"
)

require(DOSE)
dotplot(gse, showCategory=08, split=".sign") + facet_grid(.~.sign)

# Open a clean PNG file in your results folder
ggsave(filename = file.path(output_folder, "GSEA_GO_Dotplot.png"),
    width = 1200, height = 900, res = 130)

# Draw the plot
print(dotplot(gse, showCategory=08, split=".sign") + facet_grid(.~.sign))

# Close the file connection
dev.off()

# Calculate the pairwise similarity matrix (Crucial Step!)
gse_pairwise <- pairwise_termsim(gse)

# Open a clean PNG file in your results folder
png(filename = file.path(output_folder, "GSEA_Emapplot.png"),
    width = 1200, height = 900, res = 130)

# Draw the plot using the new pairwise object
print(emapplot(gse_pairwise, showCategory = 10))

# Close the file connection
dev.off()


# ------------------------------------------------------------------------------
# 10. ADVANCED ENRICHMENT: MSIGDB HALLMARK GSEA SUITE
# ------------------------------------------------------------------------------
# Download/retrieve the official Hallmark gene sets specifically for MOUSE
mouse_hallmarks <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# Run GSEA with relaxed cutoffs to see the top trends
gsea_hallmark <- GSEA(
  geneList      = gene_list,         
  TERM2GENE     = mouse_hallmarks,   
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 1,                 # Relaxed to 1 to capture ALL pathways
  pAdjustMethod = "BH",              
  verbose       = TRUE
)

# Check and Save Results
if (is.null(gsea_hallmark) || nrow(gsea_hallmark@result) == 0) {
  cat("No Hallmark pathways mapped at all.\n")
} else {
  
  # CORRECTED SORT: Sorts by absolute NES strength using base R
  gsea_hallmark@result <- gsea_hallmark@result[order(-abs(gsea_hallmark@result$NES)), ]
  
  cat("\nTop 10 Enriched Hallmark Pathways (Ranked by NES Strength):\n")
  print(head(gsea_hallmark@result[, c("Description", "NES", "pvalue", "p.adjust")], 10))
  
  # Save the full table to inspect
  write.csv(as.data.frame(gsea_hallmark),
            file = file.path(output_folder, "GSEA_Hallmark_all_results.csv"),
            row.names = FALSE)
  
  # Save the HALLMARK RIDGEPLOT
  png(filename = file.path(output_folder, "GSEA_Hallmark_Ridgeplot.png"),
      width = 1200, height = 900, res = 130)
  
  print(
    ridgeplot(gsea_hallmark, showCategory = 15) + 
      labs(x = "Log2 Fold Change") +
      ggtitle("MSigDB Hallmark GSEA: Gene Distributions\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  
  dev.off()
  
  cat("Hallmark GSEA Ridgeplot completed and saved successfully!\n")
}


# ------------------------------------------------------------------------------
# 11. ADVANCED ENRICHMENT: KEGG PATHWAY GSEA SUITE
# ------------------------------------------------------------------------------
# Map Gene Symbols to ENTREZIDs for KEGG compatibility
entrez_df <- bitr(names(gene_list), 
                  fromType = "SYMBOL", 
                  toType   = "ENTREZID", 
                  OrgDb    = org.Mm.eg.db)

# Re-align names and recreate your ranked gene list with Entrez IDs
kegg_gene_list <- gene_list[entrez_df$SYMBOL]
names(kegg_gene_list) <- entrez_df$ENTREZID
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

# Run GSEA KEGG using the mouse organism code 'mmu'
gse_kegg <- gseKEGG(
  geneList      = kegg_gene_list,
  organism      = "mmu",             # 'mmu' is the explicit KEGG code for mouse
  keyType       = "ncbi-geneid",     # Tells R we are using Entrez IDs
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 1,                 # Ensures all pathways are captured
  pAdjustMethod = "BH",
  verbose       = TRUE
)

kegg_results_df <- as.data.frame(gse_kegg)

write.csv(kegg_results_df, 
          file = file.path(output_folder, "GSEA_KEGG_Full_Results.csv"), 
          row.names = FALSE)

cat("Success: Full KEGG GSEA results table saved!\n")


# FILTER & SAVE 11 TARGETS PATHWAYS
target_pathways <- c(
  "PI3K-Akt signaling pathway", 
  "Hematopoietic cell lineage", 
  "Calcium signaling pathway", 
  "JAK-STAT signaling pathway", 
  "Ras signaling pathway", 
  "IL-17 signaling pathway", 
  "Th17 cell differentiation", 
  "cAMP signaling pathway", 
  "MAPK signaling pathway", 
  "Rap1 signaling pathway", 
  "HIF-1 signaling pathway"
)

# Extract only your 11 pathways from the master dataframe
highlighted_kegg <- kegg_results_df %>% 
  filter(Description %in% target_pathways)

# Save the subsetted 11 pathway table
write.csv(highlighted_kegg, 
          file = file.path(output_folder, "Highlighted_11_KEGG_GSEA.csv"), 
          row.names = FALSE)

# Print a preview of your 11 pathways to the console right now
cat("\n--- Preview of Your 11 Highlighted Target Pathways ---\n")
print(highlighted_kegg[, c("ID", "Description", "NES", "pvalue", "p.adjust")])


# SAVE THE HIGHLIGHTED KEGG DOTPLOT

# Ensure there are actually matching pathways to plot
if (nrow(highlighted_kegg) == 0) {
  cat("Warning: None of your 11 target pathways were found in the dataset to plot.\n")
} else {
  
  # Open a clean PNG file in your results folder
  png(filename = file.path(output_folder, "GSEA_KEGG_Highlighted_Dotplot.png"),
      width = 1200, height = 900, res = 130)
  
  # Draw the dotplot specifically using your filtered 11 pathways dataset
  print(
    dotplot(gse_kegg, showCategory = target_pathways, split = ".sign") + 
      facet_grid(.~.sign) +
      ggtitle("Target KEGG Pathways Enrichment\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
            axis.text.y = element_text(size = 10)) # Ensures pathway names remain easy to read
  )
  
  # Close the file connection
  dev.off()
  
  cat("Success: Highlighted KEGG GSEA Dotplot saved!\n")
}

# Find the unique KEGG internal code IDs for your two priority sets
il17_id <- highlighted_kegg$ID[highlighted_kegg$Description == "IL-17 signaling pathway"]
th17_id <- highlighted_kegg$ID[highlighted_kegg$Description == "Th17 cell differentiation"]

# Save IL-17 Signaling GSEA Plot
if(length(il17_id) > 0) {
  png(filename = file.path(output_folder, "GSEA_Plot_IL17_Signaling.png"), width = 900, height = 700, res = 130)
  print(gseaplot2(gse_kegg, geneSetID = il17_id, title = "KEGG: IL-17 Signaling Pathway"))
  dev.off()
}

# Save Th17 Differentiation GSEA Plot
if(length(th17_id) > 0) {
  png(filename = file.path(output_folder, "GSEA_Plot_Th17_Differentiation.png"), width = 900, height = 700, res = 130)
  print(gseaplot2(gse_kegg, geneSetID = th17_id, title = "KEGG: Th17 Cell Differentiation"))
  dev.off()
}


# ============================================================================ #






