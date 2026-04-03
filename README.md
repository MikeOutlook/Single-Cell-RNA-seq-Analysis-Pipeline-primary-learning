# 🔬 Single-Cell RNA-seq Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Platform: Linux](https://img.shields.io/badge/Platform-Linux-green.svg)](https://www.linux.org/)
[![Language: R](https://img.shields.io/badge/Language-R-red.svg)](https://www.r-project.org/)

> A powerful and user-friendly pipeline for analyzing single-cell RNA sequencing data in atrial fibrillation research.

## 📊 Overview

This pipeline provides **end-to-end single-cell RNA-seq analysis** for understanding cellular heterogeneity in atrial fibrillation (AF). It covers everything from raw data processing to differential expression analysis, featuring automated cell-type annotation and publication-ready visualizations.

### ✨ Features

- 🚀 **Fast & Efficient** - Optimized for large-scale scRNA-seq data
- 📱 **Multiple Data Formats** - 10X Genomics, H5, TXT/TSV support
- 🧬 **Smart Annotation** - Automated cell-type identification
- 📈 **Publication-Ready** - High-quality visualizations
- 🔄 **Reproducible** - Well-documented and automated workflow

---

## 🏁 Quick Start

```bash
# 1️⃣ QC & Clustering
Rscript 数据处理降维聚类.R

# 2️⃣ Cell Annotation  
Rscript cell_annotation.R

# 3️⃣ Differential Expression
Rscript 差异分析.R
```

---

## 📋 Analysis Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    scRNA-seq Analysis Workflow                              │
└─────────────────────────────────────────────────────────────────────────────┘

  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
  │   Raw Data  │───▶│     QC     │───▶│ Clustering │───▶│  UMAP/tSNE │
  │  (10X/H5)  │    │ Filtering  │    │  (Seurat)  │    │  Visualize │
  └─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
                                                        
  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
  │    Cell   │───▶│   Marker   │───▶│   Diff.    │───▶│  Enrich.   │
  │Annotation │    │  Finding   │    │   Express  │    │  (GO/KEGG)│
  └─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
```

---

## 📁 File Structure

| File | Description |
|------|-------------|
| `数据处理降维聚类.R` | QC, normalization, clustering |
| `cell_annotation.R` | Automated cell annotation |
| `细胞注释.R` | Manual marker-based annotation |
| `细胞亚型再聚类.R` | Subclustering analysis |
| `差异分析.R` | Differential expression |
| `拟时序分析.R` | Monocle3 trajectory |
| `细胞通讯及可视化.R` | CellChat communication |

---

## 🧬 Cell Types & Markers

| Cell Type | Markers | Colors |
|-----------|--------|--------|
| Osteoblastic OS cells | ALPL, RUNX2 | ████ |
| Myeloid cells | LYZ, CD68 | ████ |
| Osteoclasts OCs | ACP5, CTSK | ████ |
| CAFs | COL1A1, FAP | ████ |
| NK/T cells | CD3D, GNLY, NKG7 | ████ |
| B cells | MS4A1, CD79A | ████ |
| Plasma cells | IGHG1, MZB1 | ████ |

---

## 📈 Output Examples

### Differential Expression Results
```
data/
├── CAFs_diff.csv               # 1,735 DEGs
├── B_cells_diff.csv
├── Myeloid_cells_diff.csv
├── NK_T_cells_diff.csv
└── CAFs_DEGs.csv            # Significant DEGs (Up/Down)
```

### Visualization Outputs
```
graph/
├── 01_QCbefore.pdf          # Pre-QC plots
├── 03_QCafter.pdf          # Post-QC plots  
├── 09_CellType.pdf         # Cell type UMAP
└── Volcano_plot.pdf       # Volcano plots
```

---

## 🔧 Installation

```bash
# Option 1: Install R packages manually
Rscript install_packages.R

# Option 2: Using Conda (recommended)
conda env create -f environment.yml
conda activate scRNAseq
```

### Requirements
- R >= 4.2.0
- Seurat >= 4.0
- 16GB+ RAM recommended

---

## 📖 Usage Examples

### 1. Process 10X Data
```bash
Rscript 10X数据读取.R
Rscript 数据处理降维聚类.R
```

### 2. Annotate Cells
```bash
Rscript cell_annotation.R
```

### 3. Find DEGs
```bash
Rscript 差异分析.R
```

---

## 📊 Data Source

This pipeline was developed and tested with **GSE261170** dataset from atrial fibrillation research:

- **Samples**: 9 atrial tissue samples
- **Cells**: 66,434 single cells
- **Conditions**: EarlyAF vs PermanentAF
- **Cell Types**: 9 major types identified

---

## 📚 Citation

If you use this pipeline, please cite:

> **Seurat**: Hao Y, et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*.

> **CellChat**: Jin S, et al. (2021). Inference and analysis of cell-cell communication. *Molecular Systems Biology*.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## 📝 License

This project is licensed under the **MIT License** - see the LICENSE file for details.

---

## 👤 Author

**MikeOutlook**

- GitHub: [@MikeOutlook](https://github.com/MikeOutlook)

---

<p align="center">
  <strong>⭐ Star this repo if you find it useful! ⭐</strong>
</p>