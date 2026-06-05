# ==============================================================================
# DEG & Functional Enrichment Analysis: Mouse WT vs Tet2-KO (LSK Cells)
# WT  = Wild-Type (Reference Condition)
# KO  = Tet2 Knockout bone marrow stem cells
# ==============================================================================


# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
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

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

# Output directory
output_folder <- "DEG_Analysis_Results"
if (!dir.exists(output_folder)) dir.create(output_folder)


# ==============================================================================
# 2. DATA LOADING & EXPRESSION MATRIX
# ==============================================================================
raw_data <- read_excel("Mouse WT vs Tet2 - LSK data.xlsx", sheet = "DE_analysis")

cat("Raw dimensions:", dim(raw_data), "\n")

cleaned_counts <- raw_data %>%
  dplyr::select(GeneSymbol, starts_with("WT"), starts_with("Tet2")) %>%
  dplyr::filter(!is.na(GeneSymbol)) %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()

cat("After cleaning:", dim(cleaned_counts), "\n")

counts_matrix <- cleaned_counts %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

cat("Matrix dimensions:", dim(counts_matrix), "\n")
print(head(counts_matrix))


# ==============================================================================
# 3. SAMPLE METADATA
# ==============================================================================
sample_info <- data.frame(
  row.names = colnames(counts_matrix),
  condition = factor(
    c(rep("WT", 3), rep("Tet2", 4)),
    levels = c("WT", "Tet2")   # WT is the reference level
  )
)
# Positive LFC = upregulated in Tet2 vs WT
# Negative LFC = downregulated in Tet2 vs WT
print(sample_info)


# ==============================================================================
# 4. DESEQ2: MODEL FIT, FILTERING, RAW + SHRUNK RESULTS
# ==============================================================================
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = sample_info,
  design    = ~ condition
)

# Prefilter: keep genes with >= 10 counts in at least 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]
cat("\nGenes retained after filtering:", nrow(dds), "\n")

# Fit DESeq2 model
dds <- DESeq(dds)
resultsNames(dds)   # confirm: "condition_Tet2_vs_WT"


# ------------------------------------------------------------------------------
# SHRINKAGE NOTE:
#    apeglm-shrunk LFC used for:  Volcano plot, MA plot, heatmap, tables
#    NOT USED shrunk LFC for:  GSEA ranking (used Wald stat from raw results)
#
#   apeglm shrinks noisy, low-count genes toward LFC = 0.
#   Ideal for visualisation and reporting (removes noise artifacts),
#   For GSEA the full ranked signal across ALL genes — shrinkage
#   flattens the middle of the rank list and buries pathway-level signal.
# ------------------------------------------------------------------------------

# --- RAW results (used for GSEA ranking) ---
res_raw <- results(dds, name = "condition_Tet2_vs_WT")

# --- SHRUNK results (used for volcano, MA, heatmap, tables) ---
res_shrunk <- lfcShrink(dds, coef = "condition_Tet2_vs_WT", type = "apeglm")

# Significant DEGs table (shrunk LFC, padj < 0.05) — saved as data frame
res_sig <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  dplyr::filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

cat("\nSignificant DEGs (padj < 0.05):", nrow(res_sig), "\n")
cat("Top 10:\n")
print(head(res_sig, 10))

write.csv(res_sig,
          file = file.path(output_folder, "DEG_results_Tet2_vs_WT.csv"),
          row.names = FALSE)
cat("DEG results saved.\n")


# ==============================================================================
# 5. QC PLOTS: MA PLOT & PCA
# ==============================================================================
graphics.off()

# MA Plot — before vs after shrinkage (shrunk LFC used for display)
png(filename = file.path(output_folder, "MA_Plot.png"),
    width = 1400, height = 650, res = 140)
par(mfrow = c(1, 2))

plotMA(res_raw,
       ylim  = c(-5, 5),
       main  = "Before Shrinkage (Raw LFC)",
       xlab  = "Mean of Normalized Counts",
       ylab  = "Log2 Fold Change")
abline(h = 0, col = "red", lwd = 2)

plotMA(res_shrunk,
       ylim  = c(-5, 5),
       main  = "After Shrinkage (apeglm)",
       xlab  = "Mean of Normalized Counts",
       ylab  = "Log2 Fold Change")
abline(h = 0, col = "red", lwd = 2)

dev.off()
cat("MA plot saved.\n")

# PCA
vst_data <- vst(dds, blind = FALSE)
p_pca <- plotPCA(vst_data, intgroup = "condition") +
  theme_classic(base_size = 12) +
  ggtitle("PCA: Tet2-KO vs WT LSK Cells") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_folder, "PCA_Plot.png"),
       plot = p_pca, width = 6, height = 5, dpi = 300)
cat("PCA plot saved.\n")


# ==============================================================================
# 6. VOLCANO PLOT  (uses shrunk LFC — correct)
# ==============================================================================
res_df <- as.data.frame(res_shrunk) %>%
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

n_up    <- sum(res_df$significance == "Upregulated")
n_down  <- sum(res_df$significance == "Downregulated")
n_total <- n_up + n_down
cat("Upregulated  :", n_up,    "\n")
cat("Downregulated:", n_down,  "\n")
cat("Total DEGs   :", n_total, "\n")

# Genes to highlight on volcano
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

colors <- c(
  "Upregulated"           = "#C0392B",
  "Downregulated"         = "#2980B9",
  "Significant (|LFC|<1)" = "#E67E22",
  "NS"                    = "grey80"
)
x_max <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

build_volcano <- function(df, y_max, y_annot) {
  df <- df %>% dplyr::mutate(y_plot = pmin(-log10(padj), y_max - 5))

  ggplot(df, aes(x = log2FoldChange, y = y_plot, color = significance)) +
    geom_point(size = 1.0, alpha = 0.5) +
    geom_point(data    = dplyr::filter(df, highlight),
               size    = 2.8, shape = 21, stroke = 0.7,
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


p_main    <- build_volcano(res_df, y_max = 310, y_annot = 300)
p_trimmed <- build_volcano(res_df, y_max = 100, y_annot = 93)

ggsave(file.path(output_folder, "Volcano_Main.png"),
       plot = p_main, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_folder, "Volcano_Trimmed.png"),
       plot = p_trimmed, width = 6, height = 6, dpi = 300)
cat("Volcano plots saved.\n")

# Target gene summary table
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
print(table(target_gene_summary$significance))

write.csv(target_gene_summary,
          file = file.path(output_folder, "Target_Gene_Summary.csv"),
          row.names = FALSE)
cat("Target gene summary saved.\n")


# ==============================================================================
# 7. HEATMAP — TARGET GENES  (shrunk-normalized counts, z-scored)
# ==============================================================================
norm_counts  <- counts(dds, normalized = TRUE)
heatmap_genes <- genes_to_highlight[genes_to_highlight %in% rownames(norm_counts)]
cat("Genes in heatmap:", length(heatmap_genes), "\n")

heatmap_mat    <- log2(norm_counts[heatmap_genes, ] + 1)
heatmap_scaled <- t(scale(t(heatmap_mat)))

col_annotation <- data.frame(
  Condition = sample_info$condition,
  row.names = rownames(sample_info)
)
annotation_colors <- list(Condition = c(WT = "#4E79A7", Tet2 = "#F28E2B"))

# CORRECT device call for pheatmap
png(filename = file.path(output_folder, "Heatmap_Target_Genes.png"),
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


# ==============================================================================
# 8. KEGG ORA  (Over-Representation Analysis)
#    Input: significant DEGs (padj < 0.05, |LFC| > 1) — used shrunk LFC 
#    here because ORA only cares about significant genes, not rank.
# ==============================================================================
sig_genes_ora <- res_df %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1)

cat("\nSignificant DE genes for ORA:", nrow(sig_genes_ora), "\n")

entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys      = sig_genes_ora$gene,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)
cat("Mapped to Entrez IDs:", length(entrez_ids), "\n")

kegg_ora <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "mmu",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2
)

if (is.null(kegg_ora) || nrow(kegg_ora@result) == 0) {
  cat("No significant KEGG ORA pathways found.\n")
} else {
  cat("\nTop 15 KEGG ORA pathways:\n")
  print(head(kegg_ora@result[, c("Description", "GeneRatio", "pvalue", "p.adjust", "Count")], 15))

  write.csv(kegg_ora@result,
            file = file.path(output_folder, "KEGG_ORA_results.csv"),
            row.names = FALSE)

  png(filename = file.path(output_folder, "KEGG_ORA_dotplot.png"),
      width = 1100, height = 800, res = 130)
  print(
    dotplot(kegg_ora, showCategory = 20, font.size = 10) +
      ggtitle("KEGG ORA Enrichment\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()

  png(filename = file.path(output_folder, "KEGG_ORA_barplot.png"),
      width = 1100, height = 800, res = 130)
  print(
    barplot(kegg_ora, showCategory = 20, font.size = 10) +
      ggtitle("KEGG ORA Enrichment\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()
  cat("KEGG ORA plots saved.\n")
}


# --- 11 selected pathways dot plot (ORA) ---
pathways_11 <- c(
  "PI3K-Akt signaling pathway", "Hematopoietic cell lineage",
  "Calcium signaling pathway",  "JAK-STAT signaling pathway",
  "Ras signaling pathway",       "IL-17 signaling pathway",
  "Th17 cell differentiation",  "cAMP signaling pathway",
  "MAPK signaling pathway",     "Rap1 signaling pathway",
  "HIF-1 signaling pathway"
)

kegg_11 <- kegg_ora@result %>%
  dplyr::filter(Description %in% pathways_11) %>%
  dplyr::mutate(
    GeneRatio_numeric = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    Significant = ifelse(p.adjust < 0.05, "Significant", "Not Significant"),
    Description = fct_reorder(Description, GeneRatio_numeric)
  )

cat("\nORA — 11 pathway summary:\n")
print(kegg_11[, c("Description", "p.adjust", "Count", "Significant")])

p_kegg_11 <- ggplot(kegg_11,
                    aes(x = GeneRatio_numeric, y = Description,
                        size = Count, color = p.adjust)) +
  geom_rect(
    data        = dplyr::filter(kegg_11, Significant == "Not Significant"),
    aes(ymin = as.numeric(Description) - 0.45,
        ymax = as.numeric(Description) + 0.45,
        xmin = -Inf, xmax = Inf),
    fill = "#F2F3F4", alpha = 0.8, inherit.aes = FALSE
  ) +
  geom_point(alpha = 0.9) +
  geom_point(data  = dplyr::filter(kegg_11, Significant == "Significant"),
             shape = 21, stroke = 1.5, color = "#C0392B", fill = NA,
             aes(size = Count)) +
  geom_point(data  = dplyr::filter(kegg_11, Significant == "Not Significant"),
             shape = 21, stroke = 1.2, color = "grey50", fill = NA,
             aes(size = Count)) +
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
  scale_color_gradient(low = "#C0392B", high = "#AED6F1",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_size_continuous(name = "Gene Count", range = c(4, 13)) +
  labs(
    title    = "KEGG ORA \u2014 11 Selected Pathways",
    subtitle = "Tet2-KO vs WT LSK Cells  |  * padj < 0.05  |  ns = not significant",
    x = "Gene Ratio", y = NULL
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

ggsave(file.path(output_folder, "KEGG_ORA_11_Pathways.png"),
       plot = p_kegg_11, width = 12, height = 7, dpi = 300)
cat("KEGG ORA 11-pathway plot saved.\n")


# ==============================================================================
# 9. RANKED GENE LIST FOR GSEA
#
# KEY CORRECTION: Use DESeq2 Wald statistic from RAW results (res_raw).
#
#   Wald stat = log2FC / SE  → encodes direction + statistical confidence
#   This is the community-standard metric for GSEA ranking with DESeq2.
#
#   Do NOT use apeglm-shrunk LFC here. Shrinkage is designed for
#   visualisation/reporting — it deliberately moves genes toward zero,
#   which destroys rank-order information at the tails needed by GSEA.
#
#   Reference: Love et al. 2014 (DESeq2); Subramanian et al. 2005 (GSEA)
# ==============================================================================
res_raw_df <- as.data.frame(res_raw) %>%
  rownames_to_column("gene") %>%
  dplyr::filter(!is.na(stat))     # stat = Wald statistic

gene_list_gsea <- res_raw_df$stat
names(gene_list_gsea) <- res_raw_df$gene
gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)

cat("\nGSEA ranked gene list (Wald stat) — total genes:", length(gene_list_gsea), "\n")
cat("Range:", round(range(gene_list_gsea), 2), "\n")


# ==============================================================================
# 10. GO GSEA
# ==============================================================================
gse_go <- gseGO(
  geneList      = gene_list_gsea,
  ont           = "ALL",
  keyType       = "SYMBOL",
  minGSSize     = 10,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = TRUE,
  OrgDb         = org.Mm.eg.db
)

if (is.null(gse_go) || nrow(gse_go@result) == 0) {
  cat("No significant GO GSEA terms at padj < 0.05.\n")
} else {
  cat("\nTop GO GSEA terms:\n")
  print(head(gse_go@result[, c("Description", "NES", "pvalue", "p.adjust")], 10))

  write.csv(as.data.frame(gse_go),
            file = file.path(output_folder, "GSEA_GO_results.csv"),
            row.names = FALSE)

  # Dotplot split by activation direction
  png(filename = file.path(output_folder, "GSEA_GO_Dotplot.png"),
      width = 1200, height = 900, res = 130)
  print(
    dotplot(gse_go, showCategory = 8, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle("GO GSEA: Tet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()

  # Enrichment map
  gse_go_pw <- pairwise_termsim(gse_go)
  png(filename = file.path(output_folder, "GSEA_GO_Emapplot.png"),
      width = 1200, height = 900, res = 130)
  print(emapplot(gse_go_pw, showCategory = 10))
  dev.off()

  cat("GO GSEA plots saved.\n")
}


# ==============================================================================
# 11. HALLMARK GSEA
# ==============================================================================
mouse_hallmarks <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmark <- GSEA(
  geneList      = gene_list_gsea,
  TERM2GENE     = mouse_hallmarks,
  minGSSize     = 10,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = TRUE
)

if (is.null(gsea_hallmark) || nrow(gsea_hallmark@result) == 0) {
  cat("No significant Hallmark pathways at padj < 0.05.\n")

  # Run with pvalueCutoff = 1 to capture all trends for reporting
  cat("Re-running with pvalueCutoff = 1 to inspect trend-level results...\n")
  gsea_hallmark_all <- GSEA(
    geneList      = gene_list_gsea,
    TERM2GENE     = mouse_hallmarks,
    minGSSize     = 10,
    maxGSSize     = 800,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  gsea_hallmark_all@result <- gsea_hallmark_all@result[
    order(-abs(gsea_hallmark_all@result$NES)), ]

  cat("\nTop 10 Hallmark trends (unsupervised, all pathways):\n")
  print(head(gsea_hallmark_all@result[, c("Description", "NES", "pvalue", "p.adjust")], 10))

  write.csv(as.data.frame(gsea_hallmark_all),
            file = file.path(output_folder, "GSEA_Hallmark_all_trends.csv"),
            row.names = FALSE)

  png(filename = file.path(output_folder, "GSEA_Hallmark_Ridgeplot.png"),
      width = 1200, height = 900, res = 130)
  print(
    ridgeplot(gsea_hallmark_all, showCategory = 15) +
      labs(x = "Wald Statistic (ranked)") +
      ggtitle("MSigDB Hallmark GSEA Trends (all, unsupervised)\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()

} else {
  gsea_hallmark@result <- gsea_hallmark@result[order(-abs(gsea_hallmark@result$NES)), ]
  cat("\nTop 10 significant Hallmark pathways:\n")
  print(head(gsea_hallmark@result[, c("Description", "NES", "pvalue", "p.adjust")], 10))

  write.csv(as.data.frame(gsea_hallmark),
            file = file.path(output_folder, "GSEA_Hallmark_significant.csv"),
            row.names = FALSE)

  png(filename = file.path(output_folder, "GSEA_Hallmark_Ridgeplot.png"),
      width = 1200, height = 900, res = 130)
  print(
    ridgeplot(gsea_hallmark, showCategory = 15) +
      labs(x = "Wald Statistic (ranked)") +
      ggtitle("MSigDB Hallmark GSEA\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
  dev.off()
  cat("Hallmark GSEA plots saved.\n")
}


# ==============================================================================
# 12. KEGG GSEA
#     Input: Wald stat ranked list, mapped to Entrez IDs for KEGG
# ==============================================================================

# Map gene symbols → Entrez IDs
entrez_map <- bitr(names(gene_list_gsea),
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Mm.eg.db)

# Rebuild ranked list with Entrez IDs (preserves Wald stat rank order)
kegg_gene_list <- gene_list_gsea[entrez_map$SYMBOL]
names(kegg_gene_list) <- entrez_map$ENTREZID
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

cat("\nKEGG GSEA input gene list length:", length(kegg_gene_list), "\n")

# --- Run KEGG GSEA (standard cutoff) ---
gse_kegg <- gseKEGG(
  geneList      = kegg_gene_list,
  organism      = "mmu",
  keyType       = "ncbi-geneid",
  minGSSize     = 10,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = TRUE
)

# Save significant results
kegg_gsea_sig <- as.data.frame(gse_kegg)
write.csv(kegg_gsea_sig,
          file = file.path(output_folder, "GSEA_KEGG_Significant.csv"),
          row.names = FALSE)
cat("\nSignificant KEGG GSEA pathways (padj < 0.05):", nrow(kegg_gsea_sig), "\n")
if (nrow(kegg_gsea_sig) > 0) {
  print(kegg_gsea_sig[, c("Description", "NES", "pvalue", "p.adjust")])
}

# --- Also capture all pathways (pvalueCutoff = 1) for 11-pathway inspection ---
gse_kegg_all <- gseKEGG(
  geneList      = kegg_gene_list,
  organism      = "mmu",
  keyType       = "ncbi-geneid",
  minGSSize     = 10,
  maxGSSize     = 800,
  pvalueCutoff  = 1,          # capture everything to inspect your 11 targets
  pAdjustMethod = "BH",
  verbose       = FALSE
)

kegg_gsea_all_df <- as.data.frame(gse_kegg_all)
write.csv(kegg_gsea_all_df,
          file = file.path(output_folder, "GSEA_KEGG_All_Pathways.csv"),
          row.names = FALSE)
cat("Full KEGG GSEA results saved.\n")


# --- Extract & display your 11 target pathways regardless of significance ---
target_pathways <- c(
  "PI3K-Akt signaling pathway", "Hematopoietic cell lineage",
  "Calcium signaling pathway",  "JAK-STAT signaling pathway",
  "Ras signaling pathway",       "IL-17 signaling pathway",
  "Th17 cell differentiation",  "cAMP signaling pathway",
  "MAPK signaling pathway",     "Rap1 signaling pathway",
  "HIF-1 signaling pathway"
)

highlighted_11 <- kegg_gsea_all_df %>%
  dplyr::filter(Description %in% target_pathways) %>%
  dplyr::mutate(
    Significant_GSEA = ifelse(p.adjust < 0.05, "Yes", "No")
  ) %>%
  dplyr::select(ID, Description, setSize, NES, pvalue, p.adjust, Significant_GSEA) %>%
  arrange(pvalue)

cat("\n===== 11 TARGET PATHWAYS — KEGG GSEA (Wald stat ranked) =====\n")
print(highlighted_11, row.names = FALSE)

write.csv(highlighted_11,
          file = file.path(output_folder, "GSEA_KEGG_11_Targets.csv"),
          row.names = FALSE)
cat("11-pathway GSEA table saved.\n")


# --- Dotplot for significant KEGG GSEA pathways ---
if (nrow(kegg_gsea_sig) > 0) {
  png(filename = file.path(output_folder, "GSEA_KEGG_Dotplot.png"),
      width = 1200, height = 900, res = 130)
  print(
    dotplot(gse_kegg, showCategory = 15, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle("KEGG GSEA — Significant Pathways\nTet2-KO vs WT LSK Cells") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
            axis.text.y = element_text(size = 10))
  )
  dev.off()
  cat("KEGG GSEA dotplot saved.\n")
}


# --- Individual GSEA enrichment plots for IL-17 and Th17 ---
#     Only runs if the pathway appears in the full (pvalueCutoff=1) results.
#     This guards against the geneSetID error in the original code.

il17_id <- kegg_gsea_all_df$ID[kegg_gsea_all_df$Description == "IL-17 signaling pathway"]
th17_id <- kegg_gsea_all_df$ID[kegg_gsea_all_df$Description == "Th17 cell differentiation"]

if (length(il17_id) > 0) {
  png(filename = file.path(output_folder, "GSEA_Plot_IL17_Signaling.png"),
      width = 900, height = 700, res = 130)
  print(gseaplot2(gse_kegg_all, geneSetID = il17_id,
                  title = paste0("KEGG GSEA: IL-17 Signaling Pathway\n",
                                 "NES = ", round(kegg_gsea_all_df$NES[kegg_gsea_all_df$ID == il17_id], 3),
                                 "  padj = ", signif(kegg_gsea_all_df$p.adjust[kegg_gsea_all_df$ID == il17_id], 3))))
  dev.off()
  cat("IL-17 GSEA plot saved.\n")
} else {
  cat("IL-17 pathway not found in GSEA results.\n")
}

if (length(th17_id) > 0) {
  png(filename = file.path(output_folder, "GSEA_Plot_Th17_Differentiation.png"),
      width = 900, height = 700, res = 130)
  print(gseaplot2(gse_kegg_all, geneSetID = th17_id,
                  title = paste0("KEGG GSEA: Th17 Cell Differentiation\n",
                                 "NES = ", round(kegg_gsea_all_df$NES[kegg_gsea_all_df$ID == th17_id], 3),
                                 "  padj = ", signif(kegg_gsea_all_df$p.adjust[kegg_gsea_all_df$ID == th17_id], 3))))
  dev.off()
  cat("Th17 GSEA plot saved.\n")
} else {
  cat("Th17 pathway not found in GSEA results.\n")
}


# ==============================================================================
# 13. DIRECTION-SPLIT ORA (UP vs DOWN separately)
#     Reveals WHICH arm of each pathway is driving ORA signal.
#     Helps interpret ORA vs GSEA discordance for your PI.
# ==============================================================================
up_genes   <- res_df %>% dplyr::filter(padj < 0.05, log2FoldChange >  1) %>% pull(gene)
down_genes <- res_df %>% dplyr::filter(padj < 0.05, log2FoldChange < -1) %>% pull(gene)

map_entrez <- function(genes) {
  ids <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID",
                keytype = "SYMBOL", multiVals = "first")
  na.omit(ids)
}

up_entrez   <- map_entrez(up_genes)
down_entrez <- map_entrez(down_genes)

run_kegg_ora <- function(entrez, label) {
  res <- enrichKEGG(gene = entrez, organism = "mmu",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
  if (!is.null(res) && nrow(res@result) > 0) {
    df <- res@result %>%
      dplyr::filter(Description %in% target_pathways) %>%
      dplyr::select(Description, Count, p.adjust) %>%
      dplyr::mutate(Direction = label)
    return(df)
  }
  return(NULL)
}

ora_up   <- run_kegg_ora(up_entrez,   "Upregulated in Tet2")
ora_down <- run_kegg_ora(down_entrez, "Downregulated in Tet2")

direction_split <- bind_rows(ora_up, ora_down)

if (nrow(direction_split) > 0) {
  cat("\n===== DIRECTION-SPLIT ORA — 11 TARGET PATHWAYS =====\n")
  print(direction_split, row.names = FALSE)
  write.csv(direction_split,
            file = file.path(output_folder, "ORA_DirectionSplit_11_Pathways.csv"),
            row.names = FALSE)
  cat("Direction-split ORA saved.\n")
} else {
  cat("No target pathways significant in direction-split ORA.\n")
}



# ==============================================================================
# KEGG GSEA PLOTS — 11 TARGET PATHWAYS
# Uses gse_kegg_all (pvalueCutoff = 1) so all 11 targets are always present
# ==============================================================================

# Step 1: Subset the gse_kegg_all object to your 11 pathways only
#         This creates a valid GSEA object that enrichplot can handle
gse_kegg_11 <- gse_kegg_all

# Keep only rows matching your 11 targets in the result slot
gse_kegg_11@result <- gse_kegg_all@result %>%
  dplyr::filter(Description %in% target_pathways) %>%
  dplyr::mutate(
    # Add a significance flag column for manual annotation
    sig_label = ifelse(p.adjust < 0.05,
                       paste0(Description, " *"),
                       Description)
  ) %>%
  arrange(NES)   # sort by NES so suppressed pathways appear at top

# ------------------------------------------------------------------------------
# Plot 1: Dotplot — all 11 pathways, split by activation direction
# ------------------------------------------------------------------------------
png(filename = file.path(output_folder, "GSEA_KEGG_11_Dotplot.png"),
    width = 1300, height = 900, res = 130)

print(
  dotplot(gse_kegg_11,
          showCategory = 11,
          split        = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("KEGG GSEA — 11 Target Pathways\nTet2-KO vs WT LSK Cells") +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.text.y  = element_text(size = 10),
      strip.text   = element_text(size = 11, face = "bold")
    )
)

dev.off()
cat("11-pathway GSEA dotplot saved.\n")


# ------------------------------------------------------------------------------
# Plot 2: Custom NES bar chart — cleaner for publication
#         Shows all 11 pathways, colored by direction, sized by setSize
# ------------------------------------------------------------------------------

plot_df <- gse_kegg_11@result %>%
  dplyr::mutate(
    Direction   = ifelse(NES > 0, "Activated in Tet2-KO", "Suppressed in Tet2-KO"),
    Significant = ifelse(p.adjust < 0.05, "padj < 0.05", "ns"),
    # Star marker for significant pathways
    label       = ifelse(p.adjust < 0.05,
                         paste0(Description, "  (padj=",
                                formatC(p.adjust, format = "e", digits = 1), ") *"),
                         paste0(Description, "  (ns)")),
    Description = fct_reorder(Description, NES)
  )

p_nes <- ggplot(plot_df,
                aes(x = NES,
                    y = Description,
                    fill = Direction,
                    alpha = Significant,
                    size  = setSize)) +
  
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey40") +
  
  geom_point(shape = 21, color = "white", stroke = 0.3) +
  
  geom_text(aes(label = label,
                hjust = ifelse(NES < 0, 1.05, -0.05)),
            size = 3, color = "grey30") +
  
  scale_fill_manual(
    values = c("Activated in Tet2-KO"  = "#C0392B",
               "Suppressed in Tet2-KO" = "#2980B9"),
    name = "Direction"
  ) +
  scale_alpha_manual(
    values = c("padj < 0.05" = 0.95, "ns" = 0.35),
    name   = "Significance"
  ) +
  scale_size_continuous(
    name   = "Gene set size",
    range  = c(4, 14)
  ) +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = seq(-2.5, 2.5, by = 0.5)
  ) +
  
  labs(
    title    = "KEGG GSEA — 11 Target Pathways",
    subtitle = "Tet2-KO vs WT LSK Cells  |  * padj < 0.05  |  Ranked by Wald statistic",
    x        = "Normalized Enrichment Score (NES)",
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
    plot.margin   = margin(20, 160, 15, 20)   # right margin for labels
  )

ggsave(file.path(output_folder, "GSEA_KEGG_11_NES_Barplot.png"),
       plot = p_nes, width = 13, height = 7, dpi = 300)
cat("11-pathway NES barplot saved.\n")


# ------------------------------------------------------------------------------
# Plot 3: Individual GSEA enrichment curves for the 4 significant pathways
# ------------------------------------------------------------------------------

sig_4_ids <- gse_kegg_11@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  pull(ID)

if (length(sig_4_ids) > 0) {
  png(filename = file.path(output_folder, "GSEA_KEGG_11_EnrichmentCurves.png"),
      width = 1400, height = 900, res = 130)
  
  print(
    gseaplot2(
      gse_kegg_all,
      geneSetID = sig_4_ids,
      title     = "KEGG GSEA — Significant Target Pathways\nTet2-KO vs WT LSK Cells",
      pvalue_table = TRUE,
      ES_geom   = "line"
    )
  )
  
  dev.off()
  cat("Enrichment curves for", length(sig_4_ids), "significant pathways saved.\n")
}









# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================

# ==============================================================================
# UPDATED GENE TARGET LIST
# ==============================================================================
# Combining the 25 genes from the dot plot and the 9 added chemokines
genes_to_highlight <- c(
  # From Image (Stemness, Exhaustion, Inflammation)
  "Selp", "Clu", "Vwf", "Slamf1", "Ly6a", "Kit", "Mecom", "Junb", "Egr1", 
  "Ndn", "Cdkn1b", "Hoxb5", "Fgd5", "Ctnna1", "Rptor", "Eif4e", "Rps6", 
  "G0s2", "Cdkn1c", "Hmgb1", "S100a8", "S100a9", "Ifit1", "Isg15", "Cdkn1a",
  # Added Chemokines
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl12", "Cxcl16", "Ccl2", "Ccl3", "Ccl4", "Ccl5"
)

# ==============================================================================
# VOLCANO PLOT PREPARATION & EXECUTION
# ==============================================================================
res_df <- as.data.frame(res_shrunk) %>%
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

n_up    <- sum(res_df$significance == "Upregulated")
n_down  <- sum(res_df$significance == "Downregulated")
n_total <- n_up + n_down
cat("Upregulated  :", n_up,    "\n")
cat("Downregulated:", n_down,  "\n")
cat("Total DEGs   :", n_total, "\n")

res_df <- res_df %>%
  dplyr::mutate(
    label     = ifelse(gene %in% genes_to_highlight, gene, NA),
    highlight = gene %in% genes_to_highlight
  )

colors <- c(
  "Upregulated"           = "#C0392B",
  "Downregulated"         = "#2980B9",
  "Significant (|LFC|<1)" = "#E67E22",
  "NS"                    = "grey80"
)
x_max <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

build_volcano <- function(df, y_max, y_annot) {
  df <- df %>% dplyr::mutate(y_plot = pmin(-log10(padj), y_max - 5))
  
  ggplot(df, aes(x = log2FoldChange, y = y_plot, color = significance)) +
    geom_point(size = 1.0, alpha = 0.5) +
    geom_point(data    = dplyr::filter(df, highlight),
               size    = 2.8, shape = 21, stroke = 0.7,
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

p_main    <- build_volcano(res_df, y_max = 310, y_annot = 300)
p_trimmed <- build_volcano(res_df, y_max = 100, y_annot = 93)

ggsave(file.path(output_folder, "Volcano_Cytotoxins.png"),
       plot = p_main, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_folder, "Volcano_Trimmed_CT.png"),
       plot = p_trimmed, width = 6, height = 6, dpi = 300)
cat("Volcano plots saved.\n")

# Target gene summary table
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
print(table(target_gene_summary$significance))

write.csv(target_gene_summary,
          file = file.path(output_folder, "CT Target_Gene_Summary.csv"),
          row.names = FALSE)
cat("Target gene summary saved.\n")


# ==============================================================================
# HEATMAP — NEW TARGET GENE LIST (shrunk-normalized counts, z-scored)
# ==============================================================================
norm_counts   <- counts(dds, normalized = TRUE)
heatmap_genes <- genes_to_highlight[genes_to_highlight %in% rownames(norm_counts)]
cat("Genes in heatmap:", length(heatmap_genes), "\n")

# Check if any genes from your list were entirely missing from the dds object
missing_genes <- setdiff(genes_to_highlight, rownames(norm_counts))
if(length(missing_genes) > 0) {
  cat("Warning: The following genes were not found in the dataset expression matrix:\n", 
      paste(missing_genes, collapse = ", "), "\n")
}

heatmap_mat    <- log2(norm_counts[heatmap_genes, ] + 1)
heatmap_scaled <- t(scale(t(heatmap_mat)))

col_annotation <- data.frame(
  Condition = sample_info$condition,
  row.names = rownames(sample_info)
)
annotation_colors <- list(Condition = c(WT = "#4E79A7", Tet2 = "#F28E2B"))

# Adjusted output image height slightly since the list grew from 18 to 34 genes
png(filename = file.path(output_folder, "Heatmap_Target_Genes_CT.png"),
    width = 900, height = 1000, res = 130)

pheatmap(heatmap_scaled,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         annotation_col    = col_annotation,
         annotation_colors = annotation_colors,
         color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         border_color      = NA,
         fontsize_row      = 10,  # Adjusted for better readability of more rows
         fontsize_col      = 10,
         main              = "Target Genes: Z-score Normalized Expression\nTet2-KO vs WT LSK Cells",
         show_colnames     = TRUE,
         angle_col         = 45)
dev.off()
cat("Heatmap saved.\n")



# ==============================================================================
# HEATMAPS BY FUNCTIONAL CATEGORY
# ==============================================================================
norm_counts <- counts(dds, normalized = TRUE)

# Define the separate functional gene lists
gene_categories <- list(
  "Stemness_Quiescence" = c("Selp", "Clu", "Vwf", "Slamf1", "Ly6a", "Kit", 
                            "Mecom", "Junb", "Egr1", "Ndn", "Cdkn1b", "Hoxb5", "Fgd5"),
  "Cell_Exhaustion"     = c("Ctnna1", "Rptor", "Eif4e", "Rps6", "G0s2", "Cdkn1c"),
  "Inflammatory"        = c("Hmgb1", "S100a8", "S100a9", "Ifit1", "Isg15", "Cdkn1a"),
  "Cytokines_Chemokines"= c("Cxcl1", "Cxcl2", "Cxcl3", "Cxcl12", "Cxcl16", 
                            "Ccl2", "Ccl3", "Ccl4", "Ccl5")
)

# Shared heatmap annotations
col_annotation <- data.frame(
  Condition = sample_info$condition,
  row.names = rownames(sample_info)
)
annotation_colors <- list(Condition = c(WT = "#4E79A7", Tet2 = "#F28E2B"))

# Dynamic adjustment for row spacing based on gene counts per plot
get_heatmap_height <- function(num_genes) {
  # Base height 400px + 20px per gene (minimum 500px total)
  return(max(500, 400 + (num_genes * 20)))
}

# Loop through each category and save an individual heatmap
for (category_name in names(gene_categories)) {
  
  # Extract target genes for this specific loop iteration
  target_genes <- gene_categories[[category_name]]
  
  # Filter to keep only genes present in the dataset matrix
  present_genes <- target_genes[target_genes %in% rownames(norm_counts)]
  missing_genes <- setdiff(target_genes, rownames(norm_counts))
  
  cat("\nProcessing Category:", category_name, "\n")
  cat("Genes found:", length(present_genes), "/", length(target_genes), "\n")
  
  if(length(missing_genes) > 0) {
    cat("Warning: Missing from data matrix:", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Skip generation if no genes from the category exist in your matrix
  if(length(present_genes) == 0) {
    cat("Skipping heatmap for", category_name, "due to zero matching genes.\n")
    next
  }
  
  # Prepare the expression matrix (log2 transformed & Z-score scaled)
  heatmap_mat    <- log2(norm_counts[present_genes, , drop = FALSE] + 1)
  heatmap_scaled <- t(scale(t(heatmap_mat)))
  
  # Format clean title text (replace underscores with spaces)
  clean_title <- gsub("_", " ", category_name)
  file_name   <- paste0("Heatmap_", category_name, ".png")
  
  # Open graphics device with a dynamically adjusted layout height
  png(filename = file.path(output_folder, file_name),
      width = 850, 
      height = get_heatmap_height(length(present_genes)), 
      res = 130)
  
  # Draw the heatmap
  pheatmap(heatmap_scaled,
           cluster_rows      = TRUE,
           cluster_cols      = FALSE,
           annotation_col    = col_annotation,
           annotation_colors = annotation_colors,
           color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           border_color      = NA,
           fontsize_row      = 11,
           fontsize_col      = 10,
           main              = paste0(clean_title, " Pathway: Z-score Expression\nTet2-KO vs WT LSK Cells"),
           show_colnames     = TRUE,
           angle_col         = 45)
  
  dev.off()
  cat("Saved:", file_name, "\n")
}


# ---------------------------------------------------------------------------- #

