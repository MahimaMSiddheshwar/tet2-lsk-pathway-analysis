# Tet2-KO vs WT LSK Cells: Differential Expression & Pathway Enrichment Workflow

A robust computational biology workflow utilizing R to analyze Differential Gene Expression (DGE), Over-Representation Analysis (ORA), and Gene Set Enrichment Analysis (GSEA) on hematopoietic stem and progenitor (LSK) cells following **Tet2 Knockout (KO)**.

---

## 📦 Packages & Technologies Used

### Core Language & Environment
![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-8A2BE2?style=for-the-badge&logo=bioconductor&logoColor=white)

### Bioinformatics & Enrichment Packages
| Package | Badge | Purpose |
| :--- | :--- | :--- |
| **clusterProfiler** | ![clusterProfiler](https://img.shields.io/badge/clusterProfiler-v4.0+-blue?style=flat-square) | Core GSEA & KEGG engine execution |
| **pathview** | ![pathview](https://img.shields.io/badge/pathview-Integration-orange?style=flat-square) | Mapping expression data onto native KEGG maps |
| **org.Mm.eg.db** | ![org.Mm.eg.db](https://img.shields.io/badge/org.Mm.eg.db-Mouse-green?style=flat-square) | Alphanumeric Entrez ID & Symbol mapping |

### Data Wrangling & Visualization
| Package/Tool | Badge | Purpose |
| :--- | :--- | :--- |
| **tidyverse** | ![tidyverse](https://img.shields.io/badge/tidyverse-Data_Wrangling-blueviolet?style=flat-square) | Data manipulation (`dplyr`, `tidyr`, `purrr`) |
| **ggplot2** | ![ggplot2](https://img.shields.io/badge/ggplot2-Visualization-blue?style=flat-square) | Base canvas rendering engine |
| **enrichplot** | ![enrichplot](https://img.shields.io/badge/enrichplot-GSEA_Plots-red?style=flat-square) | Generating Ridgeplots, Dotplots, and Emapplots |
| **pheatmap** | ![pheatmap](https://img.shields.io/badge/pheatmap-Heatmaps-yellowgreen?style=flat-square) | Clustering and plotting Z-score expressions |


## 🧬 Project Context & Biological Insights (JMML Pathogenesis)

This analysis evaluates the transcriptional landscape shifts in LSK cells upon the loss of Tet2, uncovering a distinct pre-leukemic phenotype that closely mirrors the molecular priming phase of **Juvenile Myelomonocytic Leukemia (JMML)**.

### Core Core Biological Findings:
* **Epigenetic De-repression:** A massive asymmetry in differential expression (over 5,300 upregulated genes vs. ~680 downregulated genes) serves as a classic transcriptomic footprint of Tet2-loss-mediated chromatin de-repression, priming LSK cells to escape normal lineage constraints.
* **Stemness Collapse & Mobilization:** Marked downregulation of classic stem/quiescence factors (`Kit`, `Slamf1`, `Ly6a`, `Mecom`, `Egr1`) shows that Tet2-KO cells are abandoning homeostasis. Concurrently, the downregulation of `Cxcl12` indicates a disruption in bone marrow niche anchoring, a prerequisite for extramedullary migration in JMML.
* **Hyper-Metabolic Exhaustion Engine:** Upregulation of translational machinery (`Eif4e`, `Rptor`) paired with the loss of cell cycle brakes (`Cdkn1c`) demonstrates a hyper-proliferative, metabolically demanding state characteristic of expanding leukemic myelomonocytic progenitors.
* **Autonomous Myeloid-Recruiting Niche:** Robust activation of monocyte/myeloid chemoattractants (`Ccl2`, `Ccl5`, `Cxcl1`) establishes a self-sustaining pro-inflammatory autocrine loop that drives myeloid skewing.
* **Lymphoid Lineage Block:** Programmatic comparison highlights a massive coordinate suppression of alternative fates, such as the *Th17 Cell Differentiation* program ($NES = -1.935$, $padj = 0.000675$), effectively trapping progenitors in a myeloid-dominant lineage pathway.

---

## 📊 Core Workflow Outputs

The pipeline automatically processes and outputs data frames and plots to the `DEG_Analysis_Results/` directory:

1. **`Volcano_Trimmed.jpg`** - Showcases the asymmetric total DEG landscape (6,047 total DEGs; thresholds: $padj < 0.05$, $|log_2FC| > 1$) demonstrating global transcriptional de-repression.
2. **`GSEA_KEGG_11_NES_Barplot.jpg` & `KEGG_ORA_11_Pathways.jpg`** - Side-by-side methodological comparisons highlighting threshold-free rank-based GSEA performance vs. threshold-bound ORA across 11 targeted signaling cascades.
3. **`GSEA_Plot_Th17_Differentiation.png`** - Running enrichment curve demonstrating the highly significant coordinate suppression ($padj = 0.000675$) that blocks alternative lymphoid lineage potential.
4. **`GSEA_Plot_IL17_Signaling.png`** - Serving as an internal control, confirming no true alteration in the broader IL-17 pathway ($padj = 0.671$).
5. **`GSEA_Hallmark_Ridgeplot.png`** - Distribution shifts highlighting global downstream pathway activation trends across the hallmark gene sets.
6. **Targeted Condition Heatmaps (`pheatmap` Z-scores):**
    * `Heatmap_Stemness_Quiescence.png` - Captures the uniform exhaustion of self-renewal capacity across replicates.
    * `Heatmap_Cell_Exhaustion.png` - Visualizes the activation of the translational/mTOR engine driving hyper-proliferation.
    * `Heatmap_Cytokines_Chemokines.png` - Maps the autonomous production of the myeloid chemoattractant profile (*Ccl2*, *Ccl5*, *Cxcl1*).
    * `Heatmap_Inflammatory.png` - Illustrates downstream transcriptional evidence of chronic progenitor cell stress.

---

## 🚀 Getting Started

### Prerequisites
Ensure you have **R (version 4.0 or higher)** installed.

### Installation
Run the following inside your R environment to sync dependencies:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "pathview", "enrichplot"))
install.packages(c("tidyverse", "pheatmap"))

---

