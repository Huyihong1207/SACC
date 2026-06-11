# GSE294017 Cleanup Report

Run date: 2026-06-11

## Input

- File: `GSE294017_harmony.h5ad`
- Cells: `105,723`

## Filters applied

1. Basic QC
   - `n_genes_by_counts >= 500`
   - `total_counts >= 1000`
   - `pct_counts_mt < 20`
   - `pct_counts_in_top_50_genes < 45`
2. Doublet detection
   - `Scrublet` run separately for each `sample_id`
   - because the provided object is already processed and auto-thresholding was unstable, a fallback rule marked the top `6%` highest-scoring cells per sample as predicted doublets
3. Ambient RNA proxy filter
   - conservative flag for likely low-complexity, contamination-dominated nuclei:
   - `n_genes_by_counts < 800`
   - `total_counts < 1200`
   - `pct_counts_in_top_50_genes > 45`

## Summary

- Input cells: `105,723`
- Basic QC pass: `95,722`
- Predicted doublets: `6,344`
- Ambient-like cells: `3,540`
- Final kept cells: `89,436`

## Important limitation

This is not a true ambient RNA correction workflow. The input object already contains only `3000` genes and is not the original unfiltered raw droplet matrix, so methods like `CellBender` or `SoupX` cannot be applied faithfully here. The current workflow therefore uses a conservative ambient-like filtering proxy instead of full decontamination.
