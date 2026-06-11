# Paper-like Reclustering Notes

Reference article:

- [Nature Communications 2025](https://www.nature.com/articles/s41467-025-60421-0)

Targeted paper steps adapted here:

- use count matrix as input
- normalize and log-transform
- use a 3000-gene feature space
- PCA with `50` components
- Harmony integration across `sample_id`
- build neighbors on Harmony embedding
- Leiden clustering with resolution `0.5`
- compute UMAP from the Harmony neighbor graph

Important differences from the published workflow:

- the paper used `CellBender` and `scDblFinder` upstream; our local run starts from the already cleaned `GSE294017_harmony.cleaned.h5ad`
- the paper selected `3000-5000` HVGs from raw data; here the input object already contains a curated `3000`-gene set, so we use that directly
- the paper used Seurat-based preprocessing conventions; this local reproduction uses the equivalent Scanpy/Harmony implementation
