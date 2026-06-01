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

---

## 📊 Key Workflow Outputs

The pipeline automatically outputs clean visual plots and detailed tabular datasets to the `DEG_Analysis_Results/` output folder:

1. **`GSEA_Hallmark_Ridgeplot.png`** - Distribution shifts highlighting global phenotype trends (e.g., KRAS signaling, EMT).
2. **`GSEA_KEGG_Highlighted_Dotplot.png`** - Targeted evaluation screen of 11 critical signaling cascades (IL-17, Th17, Calcium, cAMP).
3. **`GSEA_Plot_IL17_Signaling.png` & `GSEA_Plot_Th17_Differentiation.png`** - Detailed traditional running enrichment curves with highlighted leading-edge lists.
4. **`Heatmap_Target_Genes.png`** - High-contrast sample-level expression metrics across specific multi-gene cascades.
5. **`GSEA_Emapplot.png`** - Functional overlap network connecting ion-exchange channels and distinct signaling components.

---

## 🚀 Getting Started

### Prerequisites
Ensure you have **R (version 4.0 or higher)** and the necessary Bioconductor manager framework installed.

### Installation
Run the following script inside your R terminal to sync all required environments:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "pathview", "enrichplot"))
install.packages(c("tidyverse", "pheatmap", "ggplot2"))
