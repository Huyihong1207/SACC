#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from __future__ import annotations

import json
from pathlib import Path

import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("/Users/niubi/Desktop/SACC")
INPUT_H5AD = ROOT / "analysis" / "gse294017_cleanup" / "output" / "paper_85142" / "GSE294017_paper_85142.h5ad"
OUTDIR = ROOT / "analysis" / "gse294017_cleanup" / "output" / "paper_85142_reclustered"

N_PCS = 50
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.5
RANDOM_STATE = 0
BATCH_KEY = "sample_id"
CLUSTER_KEY = "paper_85142_leiden"


def keep_feature(gene_name: str) -> bool:
    gene = str(gene_name).upper()
    if gene == "MALAT1":
        return False
    if gene.startswith("MT-"):
        return False
    if gene.startswith("RPL") or gene.startswith("RPS"):
        return False
    return True


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT_H5AD)
    adata_dr = adata.raw.to_adata()

    feature_mask = [keep_feature(gene) for gene in adata_dr.var_names]
    adata_dr = adata_dr[:, feature_mask].copy()
    sc.pp.highly_variable_genes(adata_dr, n_top_genes=3000, flavor="seurat")
    adata_dr = adata_dr[:, adata_dr.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(adata_dr, max_value=10)

    sc.tl.pca(
        adata_dr,
        n_comps=N_PCS,
        svd_solver="arpack",
        random_state=RANDOM_STATE,
    )
    harmony_out = hm.run_harmony(
        adata_dr.obsm["X_pca"],
        adata_dr.obs,
        BATCH_KEY,
        random_state=RANDOM_STATE,
    )
    adata_dr.obsm["X_pca_harmony"] = np.asarray(harmony_out.Z_corr)
    sc.pp.neighbors(
        adata_dr,
        n_neighbors=N_NEIGHBORS,
        n_pcs=N_PCS,
        use_rep="X_pca_harmony",
        random_state=RANDOM_STATE,
    )
    sc.tl.umap(adata_dr, random_state=RANDOM_STATE)
    sc.tl.leiden(
        adata_dr,
        resolution=LEIDEN_RESOLUTION,
        key_added=CLUSTER_KEY,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    adata.obsm["X_pca_paper_85142"] = adata_dr.obsm["X_pca"].copy()
    adata.obsm["X_pca_harmony_paper_85142"] = adata_dr.obsm["X_pca_harmony"].copy()
    adata.obsm["X_umap_paper_85142"] = adata_dr.obsm["X_umap"].copy()
    adata.obs[CLUSTER_KEY] = adata_dr.obs[CLUSTER_KEY].astype(str).values
    adata.uns[f"{CLUSTER_KEY}_colors"] = adata_dr.uns.get(f"{CLUSTER_KEY}_colors", [])

    sc.tl.rank_genes_groups(
        adata,
        groupby=CLUSTER_KEY,
        method="wilcoxon",
        use_raw=False,
    )

    adata.write(OUTDIR / "GSE294017_paper_85142.reclustered.h5ad")

    sc.pl.umap(
        adata_dr,
        color=[CLUSTER_KEY, BATCH_KEY],
        wspace=0.35,
        show=False,
    )
    plt.savefig(OUTDIR / "umap_overview.png", dpi=200, bbox_inches="tight")
    plt.close("all")

    sc.pl.umap(
        adata_dr,
        color=CLUSTER_KEY,
        legend_loc="on data",
        show=False,
    )
    plt.savefig(OUTDIR / "umap_clusters_ondata.png", dpi=220, bbox_inches="tight")
    plt.close("all")

    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_df.to_csv(OUTDIR / "markers_paper_85142_leiden.csv", index=False)

    cluster_counts = (
        adata.obs.groupby([BATCH_KEY, CLUSTER_KEY], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    cluster_counts.to_csv(OUTDIR / "cluster_counts_by_sample.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "input_cells",
                "input_genes_main_matrix",
                "input_genes_raw_matrix",
                "genes_after_mt_ribo_malat1_removal",
                "hvg_genes_used",
                "dimension_reduction_genes_removed",
                "n_pcs",
                "n_neighbors",
                "leiden_resolution",
                "n_clusters",
            ],
            "value": [
                adata.n_obs,
                adata.n_vars,
                adata.raw.shape[1],
                int(sum(feature_mask)),
                adata_dr.n_vars,
                int(len(feature_mask) - sum(feature_mask)),
                N_PCS,
                N_NEIGHBORS,
                LEIDEN_RESOLUTION,
                adata.obs[CLUSTER_KEY].nunique(),
            ],
        }
    )
    summary.to_csv(OUTDIR / "summary.tsv", sep="\t", index=False)

    settings = {
        "input_h5ad": str(INPUT_H5AD),
        "batch_key": BATCH_KEY,
        "cluster_key": CLUSTER_KEY,
        "source_matrix_for_hvg": "adata.raw",
        "n_pcs": N_PCS,
        "n_neighbors": N_NEIGHBORS,
        "leiden_resolution": LEIDEN_RESOLUTION,
        "random_state": RANDOM_STATE,
        "hvg_n_top_genes": 3000,
        "hvg_flavor": "seurat",
        "removed_feature_prefixes": ["MT-", "RPL", "RPS"],
        "removed_feature_names": ["MALAT1"],
    }
    (OUTDIR / "settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
