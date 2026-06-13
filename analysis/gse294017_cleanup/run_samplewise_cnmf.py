#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz


ROOT = Path("/Users/niubi/Desktop/SACC")
INPUT_DIR = (
    ROOT
    / "analysis"
    / "gse294017_cleanup"
    / "output"
    / "malignant_counts_by_sample"
)
OUTDIR = (
    ROOT
    / "analysis"
    / "gse294017_cleanup"
    / "output"
    / "malignant_cnmf_per_sample"
)
SUMMARY_TSV = INPUT_DIR / "extracted_samples.tsv"
GENE_MIN_CELLS = 10
N_TOP_GENES = 10000
MIN_CELLS_FOR_CNMF = 500
N_ITER = 50
MAX_NMF_ITER = 1000
SEED = 0
DENSITY_THRESHOLD = 0.5
LOCAL_NEIGHBORHOOD_SIZE = 0.30


def keep_gene(gene_name: str) -> bool:
    gene = str(gene_name).upper()
    return not (
        gene == "MALAT1"
        or gene.startswith("MT-")
        or gene.startswith("RPL")
        or gene.startswith("RPS")
    )


def k_grid_for_sample(n_cells: int) -> list[int]:
    if n_cells >= 5000:
        return [4, 5, 6, 7, 8, 9]
    if n_cells >= 3000:
        return [4, 5, 6, 7, 8]
    if n_cells >= 1000:
        return [3, 4, 5, 6, 7]
    return [2, 3, 4, 5]


def load_sample_table() -> pd.DataFrame:
    table = pd.read_csv(SUMMARY_TSV, sep="\t")
    table = table[table["n_malignant_cells_written"] >= MIN_CELLS_FOR_CNMF].copy()
    table["input_h5ad"] = table["output_h5ad"].astype(str)
    table["k_grid"] = table["n_malignant_cells_written"].map(k_grid_for_sample)
    return table.sort_values("n_malignant_cells_written", ascending=False)


def completed_k_values(sample_id: str) -> list[int]:
    sample_dir = OUTDIR / sample_id / sample_id / "cnmf_tmp"
    run_params = load_df_from_npz(
        str(OUTDIR / sample_id / sample_id / "cnmf_tmp" / f"{sample_id}.nmf_params.df.npz")
    )
    expected = (
        run_params.groupby("n_components", observed=True)
        .size()
        .rename("expected")
        .reset_index()
        .rename(columns={"n_components": "k"})
    )
    rows = []
    for k in sorted(expected["k"].astype(int).tolist()):
        observed = len(list(sample_dir.glob(f"{sample_id}.spectra.k_{k}.iter_*.df.npz")))
        rows.append({"k": k, "observed": observed})
    observed_df = pd.DataFrame(rows)
    merged = expected.merge(observed_df, on="k", how="left")
    merged["observed"] = merged["observed"].fillna(0).astype(int)
    return merged.loc[merged["observed"] >= merged["expected"], "k"].astype(int).tolist()


def build_gene_list(adata: sc.AnnData, out_path: Path) -> list[str]:
    work = adata.copy()
    work.var_names_make_unique()
    feature_mask = [keep_gene(g) for g in work.var_names]
    work = work[:, feature_mask].copy()
    sc.pp.filter_genes(work, min_cells=GENE_MIN_CELLS)

    X = work.X
    if hasattr(X, "tocsr"):
        X = X.tocsr()
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
    else:
        X = np.asarray(X)
        mean = X.mean(axis=0)
        mean_sq = np.square(X).mean(axis=0)

    var = np.maximum(mean_sq - np.square(mean), 0.0)
    dispersion = np.divide(var, mean + 1e-12, out=np.zeros_like(var), where=mean > 0)
    order = np.argsort(dispersion)[::-1]
    top_n = min(N_TOP_GENES, work.n_vars)
    genes = work.var_names[order[:top_n]].astype(str).tolist()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(genes) + "\n")
    return genes


def prepare_one(sample_id: str, input_h5ad: Path, n_cells: int, k_grid: list[int]) -> dict[str, object]:
    sample_outdir = OUTDIR / sample_id
    sample_outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_h5ad)
    genes_file = sample_outdir / f"{sample_id}.overdispersed_10000_genes.txt"
    selected_genes = build_gene_list(adata, genes_file)
    selected_mask = adata.var_names.isin(selected_genes)
    selected_view = adata[:, selected_mask]
    if hasattr(selected_view.X, "sum"):
        selected_sums = np.asarray(selected_view.X.sum(axis=1)).ravel()
    else:
        selected_sums = np.asarray(selected_view.X).sum(axis=1)
    keep_cells = selected_sums > 0
    filtered = adata[keep_cells].copy()
    filtered_h5ad = sample_outdir / f"{sample_id}.cnmf_input.h5ad"
    filtered.write(filtered_h5ad)

    runner = cNMF(output_dir=str(sample_outdir), name=sample_id)
    runner.prepare(
        counts_fn=str(filtered_h5ad),
        components=k_grid,
        n_iter=N_ITER,
        seed=SEED,
        num_highvar_genes=None,
        genes_file=str(genes_file),
        init="random",
        max_NMF_iter=MAX_NMF_ITER,
    )

    return {
        "sample_id": sample_id,
        "n_cells": n_cells,
        "n_cells_after_zero_feature_filter": int(filtered.n_obs),
        "k_grid": ",".join(map(str, k_grid)),
        "n_selected_genes": len(selected_genes),
        "genes_file": str(genes_file),
        "cnmf_input_h5ad": str(filtered_h5ad),
        "cnmf_dir": str(sample_outdir / sample_id),
        "status": "prepared",
    }


def factorize_one(sample_id: str, total_workers: int, worker_index: int) -> dict[str, object]:
    runner = cNMF(output_dir=str(OUTDIR / sample_id), name=sample_id)
    runner.factorize(worker_i=worker_index, total_workers=total_workers)
    return {"sample_id": sample_id, "status": f"factorized_worker_{worker_index}_of_{total_workers}"}


def finalize_one(sample_id: str, k_grid: list[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    runner = cNMF(output_dir=str(OUTDIR / sample_id), name=sample_id)
    runner.combine(components=k_grid)
    stats = []
    norm_counts = sc.read(runner.paths["normalized_counts"])
    for k in k_grid:
        stats.append(
            runner.consensus(
                k,
                skip_density_and_return_after_stats=True,
                show_clustering=False,
                close_clustergram_fig=True,
                norm_counts=norm_counts,
            ).stats
        )
    stats = pd.DataFrame(stats)
    stats.reset_index(drop=True, inplace=True)
    stats.to_csv(OUTDIR / sample_id / f"{sample_id}.k_selection_stats.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(stats.k, stats.silhouette, "o-", color="b")
    ax1.set_ylabel("Stability", color="b")
    for tick in ax1.get_yticklabels():
        tick.set_color("b")
    ax2.plot(stats.k, stats.prediction_error, "o-", color="r")
    ax2.set_ylabel("Error", color="r")
    for tick in ax2.get_yticklabels():
        tick.set_color("r")
    ax1.set_xlabel("Number of Components")
    ax1.grid("on")
    plt.tight_layout()
    fig.savefig(OUTDIR / sample_id / f"{sample_id}.k_selection.png", dpi=250)
    plt.close(fig)

    stats["sample_id"] = sample_id
    stats = stats.sort_values("k").reset_index(drop=True)

    max_silhouette = float(stats["silhouette"].max())
    candidate = stats[stats["silhouette"] >= max_silhouette - 0.02].copy()
    candidate = candidate.sort_values(["prediction_error", "k"], ascending=[True, True])
    selected_k = int(candidate.iloc[0]["k"])

    runner.consensus(
        selected_k,
        density_threshold=DENSITY_THRESHOLD,
        local_neighborhood_size=LOCAL_NEIGHBORHOOD_SIZE,
        show_clustering=False,
        close_clustergram_fig=True,
    )
    selected_stats = stats.loc[stats["k"] == selected_k].copy()
    selected_stats["sample_id"] = sample_id
    selected_stats["selected_k"] = selected_k
    selected_stats["k_grid"] = ",".join(map(str, k_grid))
    selected_stats["density_threshold"] = DENSITY_THRESHOLD
    selected_stats["local_neighborhood_size"] = LOCAL_NEIGHBORHOOD_SIZE
    selected_stats["consensus_spectra_txt"] = runner.paths["consensus_spectra__txt"] % (
        selected_k,
        str(DENSITY_THRESHOLD).replace(".", "_"),
    )
    selected_stats["consensus_usages_txt"] = runner.paths["consensus_usages__txt"] % (
        selected_k,
        str(DENSITY_THRESHOLD).replace(".", "_"),
    )
    return stats, selected_stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run sample-wise cNMF on malignant cells.")
    parser.add_argument(
        "--stage",
        choices=["prepare", "factorize", "finalize"],
        required=True,
    )
    parser.add_argument(
        "--sample",
        default="all",
        help="One sample_id or 'all'.",
    )
    parser.add_argument(
        "--worker-index",
        type=int,
        default=0,
        help="Worker index for factorize stage.",
    )
    parser.add_argument(
        "--total-workers",
        type=int,
        default=1,
        help="Total workers for factorize stage.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    table = load_sample_table()
    if args.sample != "all":
        table = table[table["sample_id"] == args.sample].copy()
    if table.empty:
        raise SystemExit("No eligible sample found for the requested stage.")

    if args.stage == "prepare":
        rows = []
        for row in table.itertuples(index=False):
            rows.append(
                prepare_one(
                    sample_id=str(row.sample_id),
                    input_h5ad=Path(row.input_h5ad),
                    n_cells=int(row.n_malignant_cells_written),
                    k_grid=list(row.k_grid),
                )
            )
        pd.DataFrame(rows).to_csv(OUTDIR / "prepare_summary.tsv", sep="\t", index=False)
        print(pd.DataFrame(rows).to_string(index=False))
        return

    if args.stage == "factorize":
        rows = []
        for row in table.itertuples(index=False):
            rows.append(
                factorize_one(
                    sample_id=str(row.sample_id),
                    total_workers=args.total_workers,
                    worker_index=args.worker_index,
                )
            )
        print(pd.DataFrame(rows).to_string(index=False))
        return

    all_stats = []
    selected_rows = []
    for row in table.itertuples(index=False):
        k_grid = completed_k_values(str(row.sample_id))
        if not k_grid:
            continue
        stats, selected = finalize_one(
            sample_id=str(row.sample_id),
            k_grid=k_grid,
        )
        all_stats.append(stats)
        selected_rows.append(selected)

    pd.concat(all_stats, ignore_index=True).to_csv(
        OUTDIR / "k_selection_stats.tsv",
        sep="\t",
        index=False,
    )
    pd.concat(selected_rows, ignore_index=True).to_csv(
        OUTDIR / "selected_programs.tsv",
        sep="\t",
        index=False,
    )
    print(pd.concat(selected_rows, ignore_index=True).to_string(index=False))


if __name__ == "__main__":
    main()
