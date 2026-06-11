# Paper-like Reproduction Notes

Reference article:

- [Nature Communications 2025](https://www.nature.com/articles/s41467-025-60421-0)

This local workflow follows the paper as closely as possible on the available input:

- cell filtering: `>=1000 UMIs`, `>=300 genes`
- remove cells with `>20%` mitochondrial counts
- remove cells with `>10%` ribosomal counts
- doublet removal: paper used `scDblFinder`; here `Scrublet` is used as a practical local proxy because the current machine does not have `R/scDblFinder`
- normalization + log transform
- exclude `MALAT1`, mitochondrial genes, and ribosomal genes before final dimensional reduction
- `50 PCs`
- `Harmony` integration on `sample_id`
- `Leiden resolution = 0.5`

Important non-identical parts:

- the paper used `CellBender` upstream, which cannot be reproduced from the current input because the provided file is an already processed `h5ad` rather than the raw droplet matrix
- the paper selected HVGs from a richer raw feature space; here the input already contains only `3000` genes, so this workflow uses that restricted feature set directly
