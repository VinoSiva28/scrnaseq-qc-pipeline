# Single-Cell RNA-seq QC Pipeline

A Python pipeline for performing quality control (QC) on single-cell RNA-seq count data. It computes standard per-cell QC metrics, flags low-quality cells using both fixed thresholds and an Isolation Forest anomaly detection model, and outputs summary tables, diagnostic plots, and a text report.

---

## Overview

```
Raw gene x cell count matrix (.txt)
        ↓
Load & transpose (cells as rows)
        ↓
Compute per-cell QC metrics:
  - total counts
  - genes detected
  - mitochondrial percentage
        ↓
Flag low-quality cells:
  - fixed-threshold QC
  - Isolation Forest anomaly detection
        ↓
Save QC tables + plots + text report
```

## Features

- Loads a tab-separated gene expression count matrix and reshapes it to one row per cell.
- Calculates total counts per cell.
- Calculates the number of genes detected per cell.
- Identifies mitochondrial genes (columns prefixed with `MT-`) and calculates the percentage of each cell's counts coming from mitochondrial genes.
- Flags cells as pass/fail QC using configurable thresholds (mitochondrial %, genes detected, total counts).
- Alternative anomaly-based QC flagging using an Isolation Forest model on total counts, genes detected, and mitochondrial percentage.
- Generates distribution plots (with cutoff lines) for mitochondrial percentage, total counts, and genes detected.
- Generates a bar plot of QC status counts (pass vs. fail).
- Saves QC summary tables, a filtered "passed cells" table, and a plain-text summary report.

## Repository Contents

| File/Folder | Description |
|---|---|
| `main.py` | Entry point — runs the full QC pipeline end to end. |
| `data/CRC_subset.txt` | Input gene x cell count matrix (tab-separated, genes as rows). |
| `results/` | Output directory for QC tables, plots, and the text report. |

## Software Requirements

- Python 3.x
- Python packages: `pandas`, `matplotlib`, `seaborn`, `scikit-learn`

## Usage

Place your input count matrix at `data/CRC_subset.txt` (tab-separated, genes as rows and cells as columns), then run:

```bash
python main.py
```

This will create a `results/` folder (if it doesn't already exist) containing:

- `QC_summary.csv` — per-cell QC metrics and pass/fail status.
- `passed_cells.csv` — only the cells that passed QC.
- `QC_passed_summary.csv` — full per-cell metrics table after QC assignment.
- `qc_report.txt` — plain-text summary (total cells, passed/failed counts, average mitochondrial %, average total counts).
- `mitochondrial_percentage_plot.png`, `total_counts_plot.png`, `genes_detected_plot.png` — distribution plots with QC cutoff lines.
- `qc_status_counts.png` — bar plot of the number of cells that passed vs. failed QC.

## QC Approach

Two QC flagging methods are available:

1. **Threshold-based QC** (`assign_qc_status`): flags a cell as failing QC if its mitochondrial percentage exceeds a threshold (default 20%), its genes detected falls below a threshold (default 200), or its total counts falls below a threshold (default 500).
2. **Isolation Forest QC** (`assign_qc_status_isolation_forest`): fits an Isolation Forest on total counts, genes detected, and mitochondrial percentage to flag anomalous cells (default contamination rate 0.1), without relying on fixed cutoffs.

The pipeline currently runs the Isolation Forest method by default in `main()`.
