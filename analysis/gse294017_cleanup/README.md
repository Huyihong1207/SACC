# GSE294017 Cleanup Workflow

This folder contains a reproducible cleanup workflow for the user-provided `GSE294017_harmony.h5ad`.

## What this workflow does

- reads the provided Harmony-integrated AnnData object
- joins sample-level metadata from GEO
- applies conservative QC filters
- runs Scrublet per sample to flag likely doublets
- applies an additional conservative `ambient_like` filter for low-complexity, likely contamination-dominated nuclei

## Important limitation

This input file already contains only `3000` genes and is not the original raw droplet matrix. Because of that:

- true ambient RNA correction methods such as CellBender or SoupX cannot be run faithfully here
- the `ambient_like` step in this project is a conservative filtering proxy, not a full decontamination algorithm

If you later provide the raw unfiltered 10x matrix, we can upgrade this workflow to a true ambient RNA correction pipeline.

## Main script

```bash
source /Users/niubi/Desktop/SACC/.venv-scanpy/bin/activate
python /Users/niubi/Desktop/SACC/analysis/gse294017_cleanup/run_gse294017_cleanup.py
```

## Outputs

- `analysis/gse294017_cleanup/output/GSE294017_harmony.annotated_qc.h5ad`
- `analysis/gse294017_cleanup/output/GSE294017_harmony.cleaned.h5ad`
- `analysis/gse294017_cleanup/output/qc_summary.tsv`
- `analysis/gse294017_cleanup/output/sample_doublet_summary.tsv`

