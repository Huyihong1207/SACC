#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from pathlib import Path

import harmonypy as hm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


ROOT = Path("/Users/niubi/Desktop/SACC")
INPUT_H5AD = ROOT / "analysis" / "gse294017_cleanup" / "output" / "GSE294017_harmony.cleaned.h5ad"
OUTDIR = ROOT / "analysis" / "gse294017_cleanup" / "output" / "reclustered"
OUTPUT_H5AD = OUTDIR / "GSE294017_harmony.cleaned.reclustered.h5ad"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    adata = sc.read_h5ad(INPUT_H5AD)
    adata.X = adata.layers["counts"].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # The cleaned object already contains only a curated 3000-gene feature set,
    # which we use as a proxy for the paper's 3000-5000 HVGs.
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample_id", verbose=True)
    adata.obsm["X_pca_harmony_recluster"] = ho.Z_corr
    sc.pp.neighbors(adata, use_rep="X_pca_harmony_recluster", n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_recluster")
    sc.tl.rank_genes_groups(adata, "leiden_recluster", method="wilcoxon", key_added="rank_genes_recluster")

    adata.write(OUTPUT_H5AD)

    sc.pl.umap(adata, color=["sample_id", "tumor_entity", "leiden_recluster"], show=False)
    plt.savefig(OUTDIR / "umap_recluster_overview.png", bbox_inches="tight", dpi=180)
    plt.close("all")

    acc = adata[adata.obs["tumor_entity"].astype(str).str.contains("Adenoid cystic carcinoma", case=False, na=False)].copy()
    if acc.n_obs > 0:
        sc.pl.umap(acc, color=["sample_id", "leiden_recluster"], show=False)
        plt.savefig(OUTDIR / "umap_recluster_acc_only.png", bbox_inches="tight", dpi=180)
        plt.close("all")

    sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_recluster").to_csv(
        OUTDIR / "markers_leiden_recluster.csv", index=False
    )
    adata.obs.groupby(["sample_id", "tumor_entity", "leiden_recluster"], observed=True).size().rename("n_cells").reset_index().to_csv(
        OUTDIR / "cluster_counts_by_sample.tsv", sep="\t", index=False
    )

    print("output_cells", adata.n_obs)
    print("n_clusters", adata.obs["leiden_recluster"].nunique())


if __name__ == "__main__":
    main()
