# Single-Cell RNA-seq Analysis Pipeline

A comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline for atrial fibrillation (AF) research.

## Overview

This pipeline provides end-to-end workflows for single-cell RNA-seq data analysis, including:

- **Data Import**: Support for 10X Genomics, H5, and TXT/TSV formats
- **QC & Clustering**: Quality control, normalization, dimensionality reduction (PCA/UMAP), and clustering
- **Cell Annotation**: Both manual and automatic cell type annotation based on marker genes
- **Subclustering**: Detailed subclustering of specific cell types
- **Differential Expression**: Differential gene expression analysis between conditions
- **Trajectory Analysis**: Monocle3-based pseudotime trajectory analysis
- **Cell-Cell Communication**: CellChat-based intercellular communication analysis

## Quick Start

```bash
# Install dependencies
Rscript install_packages.R

# Run analysis pipeline
Rscript 数据处理降维聚类.R     # Step 1: QC & Clustering
Rscript cell_annotation.R     # Step 2: Cell Annotation (or: Rscript 细胞注释.R)
Rscript 差异分析.R            # Step 3: Differential Expression
```

## Directory Structure

```
.
├── 10X数据读取.R           # Read 10X Genomics data
├── h5数据读取.R            # Read H5 format data
├── txt或tsv数据读取.R       # Read TXT/TSV format data
├── 数据处理降维聚类.R      # QC & Clustering
├── 细胞注释.R              # Manual cell annotation
├── cell_annotation.R     # Automatic cell annotation
├── 细胞亚型再聚类.R        # Subclustering
├── 差异分析.R            # Differential expression
├── continue_analysis.R   # Visualization
├── 拟时序分析.R          # Trajectory analysis (Monocle3)
├── 细胞通讯及可视化.R      # Cell communication (CellChat)
├── data/                 # Differential expression results (*_diff.csv)
└── graph/                # Visualization outputs (*.pdf, *.tiff)
```

## Data Requirements

### Input Data Format
- **10X**: `matrix.mtx`, `features.tsv`, `barcodes.tsv`
- **H5**: `.h5` file
- **TXT**: Expression matrix (rows=genes, columns=cells)

### Sample Metadata
Required columns in metadata:
- `sample`: Sample ID
- `condition`: Group info (e.g., EarlyAF, PermanentAF)

## Analysis Results

### Differential Expression Output (data/)
- `{cell_type}_diff.csv` - Differential genes for each cell type
- `CAFs_DEGs.csv` - DEGs with up/down regulation
- `allmarkers.csv` - All marker genes

### Visualization Output (graph/)
- `01_QCbefore.*` - QC plots before filtering
- `02_FeatureScatter.*` - Feature scatter plots
- `03_QCafter.*` - QC plots after filtering
- `04_RidgePlot.*` - Cell cycle scores
- `09_CellType.*` - Cell type UMAP
- `Volcano_plot.*` - Volcano plots

## Dependencies

### R Version
- R >= 4.2.0

### Required R Packages
- Seurat (>= 4.0)
- tidyverse
- ggplot2
- dplyr
- monocle3
- CellChat
- clusterProfiler

### Installation
```bash
# Using conda (recommended)
conda env create -f environment.yml

# Or install R packages manually
Rscript install_packages.R
```

## Cell Type Markers

| Cell Type | Marker Genes |
|----------|-------------|
| Osteoblastic OS cells | ALPL, RUNX2, IBSP |
| Myeloid cells | LYZ, CD68 |
| Osteoclasts OCs | ACP5, CTSK |
| CAFs | COL1A1, FAP, VIM |
| NK/T cells | CD2, CD3D, GNLY, NKG7 |
| Endothelial cells | EGFL7, PLVAP |
| B cells | MS4A1, CD79A |
| Plasma cells | IGHG1, MZB1 |

## License

MIT License

## Citation

If you use this pipeline, please cite:
- Seurat: Hao Y, et al. (2021)
- monocle3: Trapnell C, et al.
- CellChat: Jin S, et al. (2021)

## Author

MikeOutlook

## Contact

For questions or issues, please open an issue on GitHub.