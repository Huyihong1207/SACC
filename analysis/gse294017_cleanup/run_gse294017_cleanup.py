#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from __future__ import annotations

from pathlib import Path
import json
import re

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
OUTDIR = ROOT / "analysis" / "gse294017_cleanup" / "output"

MIN_GENES = 500
MIN_COUNTS = 1000
MAX_MT = 20.0
MAX_TOP50 = 45.0
AMBIENT_MAX_GENES = 800
AMBIENT_MAX_COUNTS = 1200
AMBIENT_MIN_TOP50 = 45.0
EXPECTED_DOUBLET_RATE = 0.06


def parse_pid(sample_id: str) -> str:
    m = re.match(r"GSM\d+_(.+)$", sample_id)
    return m.group(1) if m else sample_id


def load_data():
    adata = sc.read_h5ad(INPUT_H5AD)
    if "counts" not in adata.layers:
        raise ValueError("Input h5ad does not contain a counts layer")
    meta = pd.read_csv(META_CSV, sep="\t").rename(
        columns={"PID": "pid", "Biopsy site": "biopsy_site", "Tumor entity": "tumor_entity"}
    )
    adata.obs["pid"] = adata.obs["sample_id"].astype(str).map(parse_pid)
    adata.obs = adata.obs.join(meta.set_index("pid"), on="pid")
    return adata


def add_qc_metrics(adata):
    # Preserve existing QC fields from the input object; these are more reliable than
    # recomputing percentages from an already HVG-restricted matrix.
    original_pct_mt = adata.obs["pct_counts_mt"].copy() if "pct_counts_mt" in adata.obs.columns else None
    original_n_genes = adata.obs["n_genes_by_counts"].copy() if "n_genes_by_counts" in adata.obs.columns else None
    original_total_counts = adata.obs["total_counts"].copy() if "total_counts" in adata.obs.columns else None
    original_top50 = adata.obs["pct_counts_in_top_50_genes"].copy() if "pct_counts_in_top_50_genes" in adata.obs.columns else None

    counts = adata.layers["counts"].copy()
    adata.X = counts
    adata.var["mt_gene"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo_gene"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb_gene"] = adata.var_names.str.upper().str.match(r"^HB[ABDEGMQZ]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt_gene", "ribo_gene", "hb_gene"], inplace=True)

    if original_pct_mt is not None:
        adata.obs["pct_counts_mt_used"] = original_pct_mt
    else:
        adata.obs["pct_counts_mt_used"] = adata.obs["pct_counts_mt_gene"]
    if original_n_genes is not None:
        adata.obs["n_genes_by_counts_used"] = original_n_genes
    else:
        adata.obs["n_genes_by_counts_used"] = adata.obs["n_genes_by_counts"]
    if original_total_counts is not None:
        adata.obs["total_counts_used"] = original_total_counts
    else:
        adata.obs["total_counts_used"] = adata.obs["total_counts"]
    if original_top50 is not None:
        adata.obs["pct_counts_in_top_50_genes_used"] = original_top50
    else:
        adata.obs["pct_counts_in_top_50_genes_used"] = adata.obs["pct_counts_in_top_50_genes"]
    return adata


def initial_qc_flags(adata):
    adata.obs["qc_pass_basic"] = (
        (adata.obs["n_genes_by_counts_used"] >= MIN_GENES)
        & (adata.obs["total_counts_used"] >= MIN_COUNTS)
        & (adata.obs["pct_counts_mt_used"] < MAX_MT)
        & (adata.obs["pct_counts_in_top_50_genes_used"] < MAX_TOP50)
    )
    adata.obs["ambient_like"] = (
        (adata.obs["n_genes_by_counts_used"] < AMBIENT_MAX_GENES)
        & (adata.obs["total_counts_used"] < AMBIENT_MAX_COUNTS)
        & (adata.obs["pct_counts_in_top_50_genes_used"] > AMBIENT_MIN_TOP50)
    )
    return adata


def run_scrublet_by_sample(adata):
    scores = pd.Series(index=adata.obs_names, dtype=float)
    preds = pd.Series(index=adata.obs_names, dtype=object)
    thresholds = []
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
        threshold_source = "auto"
        if threshold is None or pd.isna(threshold):
            n_doublets = max(1, int(round(len(sample_scores) * EXPECTED_DOUBLET_RATE)))
            order = np.argsort(sample_scores)
            fallback_preds = np.zeros(len(sample_scores), dtype=bool)
            fallback_preds[order[-n_doublets:]] = True
            sample_preds = fallback_preds
            threshold = float(np.min(np.asarray(sample_scores)[fallback_preds]))
            threshold_source = "fallback_top_rate"
        scores.loc[subset.obs_names] = sample_scores
        preds.loc[subset.obs_names] = sample_preds
        thresholds.append(
            {
                "sample_id": sample_id,
                "n_cells": subset.n_obs,
                "scrublet_threshold": threshold,
                "threshold_source": threshold_source,
                "predicted_doublets": int(np.sum(sample_preds)),
            }
        )
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = preds.fillna(False).astype(bool)
    return adata, pd.DataFrame(thresholds)


def finalize_filters(adata):
    adata.obs["keep_cell"] = (
        adata.obs["qc_pass_basic"]
        & (~adata.obs["predicted_doublet"])
        & (~adata.obs["ambient_like"])
    )
    return adata


def save_summary(adata, sample_summary):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    qc_summary = pd.DataFrame(
        {
            "metric": [
                "input_cells",
                "basic_qc_pass",
                "predicted_doublets",
                "ambient_like",
                "final_kept",
            ],
            "value": [
                int(adata.n_obs),
                int(adata.obs["qc_pass_basic"].sum()),
                int(adata.obs["predicted_doublet"].sum()),
                int(adata.obs["ambient_like"].sum()),
                int(adata.obs["keep_cell"].sum()),
            ],
        }
    )
    qc_summary.to_csv(OUTDIR / "qc_summary.tsv", sep="\t", index=False)
    sample_summary.to_csv(OUTDIR / "sample_doublet_summary.tsv", sep="\t", index=False)

    by_sample = (
        adata.obs.groupby("sample_id")[["qc_pass_basic", "predicted_doublet", "ambient_like", "keep_cell"]]
        .sum()
        .reset_index()
    )
    by_sample.to_csv(OUTDIR / "sample_filter_summary.tsv", sep="\t", index=False)

    settings = {
        "MIN_GENES": MIN_GENES,
        "MIN_COUNTS": MIN_COUNTS,
        "MAX_MT": MAX_MT,
        "MAX_TOP50": MAX_TOP50,
        "AMBIENT_MAX_GENES": AMBIENT_MAX_GENES,
        "AMBIENT_MAX_COUNTS": AMBIENT_MAX_COUNTS,
        "AMBIENT_MIN_TOP50": AMBIENT_MIN_TOP50,
        "EXPECTED_DOUBLET_RATE": EXPECTED_DOUBLET_RATE,
    }
    (OUTDIR / "run_settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")


def save_plots(adata):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    metrics = [
        "n_genes_by_counts_used",
        "total_counts_used",
        "pct_counts_mt_used",
        "pct_counts_in_top_50_genes_used",
    ]
    for ax, metric in zip(axes, metrics):
        ax.hist(adata.obs[metric].to_numpy(), bins=100, color="#4C78A8")
        ax.set_title(metric)
    fig.tight_layout()
    fig.savefig(OUTDIR / "qc_histograms.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(adata.obs["doublet_score"].dropna().to_numpy(), bins=80, color="#F58518")
    ax.set_title("Scrublet Doublet Score")
    fig.tight_layout()
    fig.savefig(OUTDIR / "scrublet_scores.png", dpi=180)
    plt.close(fig)

    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color=["sample_id", "predicted_doublet", "ambient_like", "keep_cell"], show=False)
        plt.savefig(OUTDIR / "umap_cleanup_flags.png", bbox_inches="tight", dpi=180)
        plt.close("all")


def save_h5ad(adata):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    annotated = OUTDIR / "GSE294017_harmony.annotated_qc.h5ad"
    cleaned = OUTDIR / "GSE294017_harmony.cleaned.h5ad"
    adata.write(annotated)
    adata[adata.obs["keep_cell"]].copy().write(cleaned)


def main():
    sc.settings.verbosity = 2
    adata = load_data()
    adata = add_qc_metrics(adata)
    adata = initial_qc_flags(adata)
    adata, sample_summary = run_scrublet_by_sample(adata)
    adata = finalize_filters(adata)
    save_summary(adata, sample_summary)
    save_plots(adata)
    save_h5ad(adata)
    print("input_cells", adata.n_obs)
    print("kept_cells", int(adata.obs["keep_cell"].sum()))


if __name__ == "__main__":
    main()
