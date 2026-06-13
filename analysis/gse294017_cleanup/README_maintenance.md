Workspace maintenance rules for `gse294017_cleanup`:

1. Large matrices and cNMF working directories stay in `analysis/gse294017_cleanup/output/` and are ignored by Git.
2. Small summaries that are worth versioning are copied into `analysis/gse294017_cleanup/reports/`.
3. Safe junk files are removed by `maintain_workspace.py`:
   - `.DS_Store`
   - `__pycache__/`
   - `*.pyc`
   - temporary `*.zero_feature_filtered.h5ad` files after a successful cNMF prepare

Recommended usage:

```bash
source /Users/niubi/Desktop/SACC/.venv-scanpy/bin/activate
python /Users/niubi/Desktop/SACC/analysis/gse294017_cleanup/maintain_workspace.py
```

Suggested Git policy:

- Track scripts in `analysis/gse294017_cleanup/`
- Track reports in `analysis/gse294017_cleanup/reports/`
- Do not track `.h5ad`, `.h5`, or `output/` contents
