#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from pathlib import Path
import json

import harmonypy as hm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scrublet as scr


ROOT = Path("/Users/niubi/Desktop/SACC")
INPUT_H5AD = ROOT / "GSE294017_harmony.h5ad"
META_CSV = ROOT / "GSE294017" / "GSE294017_Dataframe_samples_single_nuclei_GEO.csv.gz"
OUTDIR = ROOT / "analysis" / "gse294017_cleanup" / "output" / "paper_like"
OUTPUT_H5AD = OUTDIR / "GSE294017_paper_like.h5ad"

MIN_COUNTS = 1000
MIN_GENES = 300
MAX_MT = 20.0
MAX_RIBO = 10.0
EXPECTED_DOUBLET_RATE = 0.10
LEIDEN_RESOLUTION = 0.5


def parse_pid(sample_id: str) -> str:
    parts = str(sample_id).split("_")
    return "_".join(parts[1:]) if len(parts) > 1 else str(sample_id)


def add_qc_metrics(adata):
    meta = pd.read_csv(META_CSV, sep="\t").rename(
        columns={"PID": "pid", "Biopsy site": "biopsy_site", "Tumor entity": "tumor_entity"}
    )
    adata.obs["pid"] = adata.obs["sample_id"].astype(str).map(parse_pid)
    adata.obs = adata.obs.join(meta.set_index("pid"), on="pid")

    adata.X = adata.layers["counts"].copy()
    adata.var["mt_gene"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo_gene"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["malat1_gene"] = adata.var_names.str.upper() == "MALAT1"
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt_gene", "ribo_gene"], inplace=True)
    return adata


def run_scrublet_by_sample(adata):
    scores = pd.Series(index=adata.obs_names, dtype=float)
    preds = pd.Series(index=adata.obs_names, dtype=object)
    summary = []
    for sample_id, idx in adata.obs.groupby("sample_id", observed=False).groups.items():
        subset = adata[idx].copy()
        matrix = subset.layers["counts"]
        if not sp.issparse(matrix):
            matrix = sp.csr_matrix(matrix)
        scrub = scr.Scrublet(matrix, expected_doublet_rate=EXPECTED_DOUBLET_RATE)
        sample_scores, sample_preds = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
        )
        threshold = getattr(scrub, "threshold_", None)
        source = "auto"
        if threshold is None or pd.isna(threshold):
            n_doublets = max(1, int(round(len(sample_scores) * EXPECTED_DOUBLET_RATE)))
            order = np.argsort(sample_scores)
            fallback = np.zeros(len(sample_scores), dtype=bool)
            fallback[order[-n_doublets:]] = True
            sample_preds = fallback
            threshold = float(np.min(np.asarray(sample_scores)[fallback]))
            source = "fallback_top_rate"
        scores.loc[subset.obs_names] = sample_scores
        preds.loc[subset.obs_names] = sample_preds
        summary.append(
            {
                "sample_id": sample_id,
                "n_cells": subset.n_obs,
                "doublet_threshold": threshold,
                "threshold_source": source,
                "predicted_doublets": int(np.sum(sample_preds)),
            }
        )
    adata.obs["doublet_score_paper_like"] = scores
    adata.obs["predicted_doublet_paper_like"] = preds.fillna(False).astype(bool)
    return adata, pd.DataFrame(summary)


def filter_cells(adata):
    adata.obs["paper_like_qc_pass"] = (
        (adata.obs["total_counts"] >= MIN_COUNTS)
        & (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["pct_counts_mt_gene"] < MAX_MT)
        & (adata.obs["pct_counts_ribo_gene"] < MAX_RIBO)
        & (~adata.obs["predicted_doublet_paper_like"])
    )
    return adata[adata.obs["paper_like_qc_pass"]].copy()


def preprocess_like_paper(adata):
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    if "log1p" in adata.uns:
        del adata.uns["log1p"]
    sc.pp.log1p(adata)
    adata.raw = adata

    keep_genes = ~(adata.var["mt_gene"] | adata.var["ribo_gene"] | adata.var["malat1_gene"])
    adata = adata[:, keep_genes].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample_id", verbose=True)
    adata.obsm["X_pca_harmony_paper_like"] = ho.Z_corr
    sc.pp.neighbors(adata, use_rep="X_pca_harmony_paper_like", n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added="leiden_paper_like")
    sc.tl.rank_genes_groups(adata, "leiden_paper_like", method="wilcoxon", key_added="rank_genes_paper_like")
    return adata


def save_outputs(adata, dbl_summary, input_n):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata.write(OUTPUT_H5AD)
    dbl_summary.to_csv(OUTDIR / "paper_like_doublet_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "metric": [
                "input_cells",
                "kept_after_paper_like_qc",
                "n_clusters",
                "n_genes_after_feature_exclusion",
            ],
            "value": [
                int(input_n),
                int(adata.n_obs),
                int(adata.obs["leiden_paper_like"].nunique()),
                int(adata.n_vars),
            ],
        }
    ).to_csv(OUTDIR / "paper_like_summary.tsv", sep="\t", index=False)
    settings = {
        "MIN_COUNTS": MIN_COUNTS,
        "MIN_GENES": MIN_GENES,
        "MAX_MT": MAX_MT,
        "MAX_RIBO": MAX_RIBO,
        "EXPECTED_DOUBLET_RATE": EXPECTED_DOUBLET_RATE,
        "LEIDEN_RESOLUTION": LEIDEN_RESOLUTION,
    }
    (OUTDIR / "paper_like_settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")

    sc.pl.umap(adata, color=["sample_id", "tumor_entity", "leiden_paper_like"], show=False, use_raw=False)
    plt.savefig(OUTDIR / "umap_paper_like_overview.png", bbox_inches="tight", dpi=180)
    plt.close("all")

    acc = adata[adata.obs["tumor_entity"].astype(str).str.contains("Adenoid cystic carcinoma", case=False, na=False)].copy()
    if acc.n_obs > 0:
        sc.pl.umap(acc, color=["sample_id", "leiden_paper_like"], show=False, use_raw=False)
        plt.savefig(OUTDIR / "umap_paper_like_acc_only.png", bbox_inches="tight", dpi=180)
        plt.close("all")

    sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_paper_like").to_csv(
        OUTDIR / "markers_leiden_paper_like.csv", index=False
    )


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    adata = sc.read_h5ad(INPUT_H5AD)
    input_n = adata.n_obs
    adata = add_qc_metrics(adata)
    adata, dbl_summary = run_scrublet_by_sample(adata)
    adata = filter_cells(adata)
    adata = preprocess_like_paper(adata)
    save_outputs(adata, dbl_summary, input_n)
    print("input_cells", input_n)
    print("kept_cells", adata.n_obs)
    print("n_clusters", adata.obs["leiden_paper_like"].nunique())


if __name__ == "__main__":
    main()
