#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd
import scanpy as sc


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
MIN_CELLS_FOR_CNMF = 500
N_ITER = 50
NUM_GENES = 3000
MAX_NMF_ITER = 1000
SEED = 0


def k_grid_for_sample(n_cells: int) -> list[int]:
    if n_cells >= 5000:
        return [4, 5, 6]
    if n_cells >= 3000:
        return [4, 5, 6]
    if n_cells >= 1000:
        return [3, 4, 5]
    return [2, 3, 4]


def load_sample_table() -> pd.DataFrame:
    table = pd.read_csv(SUMMARY_TSV, sep="\t")
    table = table[table["n_malignant_cells_written"] >= MIN_CELLS_FOR_CNMF].copy()
    table["k_grid"] = table["n_malignant_cells_written"].map(k_grid_for_sample)
    return table.sort_values("n_malignant_cells_written", ascending=False)


def run_prepare(sample_id: str, counts_h5ad: str, k_grid: list[int]) -> dict[str, object]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sample_root = OUTDIR / sample_id / sample_id
    norm_counts = sample_root / "cnmf_tmp" / f"{sample_id}.norm_counts.h5ad"
    genes_file = sample_root / f"{sample_id}.overdispersed_genes.txt"
    if norm_counts.exists():
        return {
            "sample_id": sample_id,
            "counts_h5ad": counts_h5ad,
            "k_grid": ",".join(map(str, k_grid)),
            "n_iter": N_ITER,
            "num_genes": NUM_GENES,
            "status": "already_prepared",
        }

    cmd = [
        str(ROOT / ".venv-scanpy" / "bin" / "cnmf"),
        "prepare",
        "--output-dir",
        str(OUTDIR / sample_id),
        "--name",
        sample_id,
        "--counts",
        counts_h5ad,
        "--components",
        *[str(k) for k in k_grid],
        "--n-iter",
        str(N_ITER),
        "--seed",
        str(SEED),
        "--numgenes",
        str(NUM_GENES),
        "--max-nmf-iter",
        str(MAX_NMF_ITER),
        "--init",
        "random",
    ]
    try:
        subprocess.run(cmd, check=True)
        status = "prepared"
    except subprocess.CalledProcessError:
        if not genes_file.exists():
            raise
        filtered_h5ad = filter_zero_feature_cells(
            sample_id=sample_id,
            counts_h5ad=Path(counts_h5ad),
            genes_file=genes_file,
        )
        fallback_cmd = [
            str(ROOT / ".venv-scanpy" / "bin" / "cnmf"),
            "prepare",
            "--output-dir",
            str(OUTDIR / sample_id),
            "--name",
            sample_id,
            "--counts",
            str(filtered_h5ad),
            "--components",
            *[str(k) for k in k_grid],
            "--n-iter",
            str(N_ITER),
            "--seed",
            str(SEED),
            "--genes-file",
            str(genes_file),
            "--max-nmf-iter",
            str(MAX_NMF_ITER),
            "--init",
            "random",
        ]
        subprocess.run(fallback_cmd, check=True)
        status = "prepared_after_zero_feature_filter"
    return {
        "sample_id": sample_id,
        "counts_h5ad": counts_h5ad,
        "k_grid": ",".join(map(str, k_grid)),
        "n_iter": N_ITER,
        "num_genes": NUM_GENES,
        "status": status,
    }


def filter_zero_feature_cells(sample_id: str, counts_h5ad: Path, genes_file: Path) -> Path:
    adata = sc.read_h5ad(counts_h5ad)
    genes = [line.strip() for line in genes_file.read_text().splitlines() if line.strip()]
    keep_genes = [g for g in genes if g in adata.var_names]
    view = adata[:, keep_genes]
    cell_sums = view.X.sum(axis=1)
    if hasattr(cell_sums, "A1"):
        cell_sums = cell_sums.A1
    kept = adata[cell_sums > 0].copy()
    out = OUTDIR / sample_id / f"{sample_id}.zero_feature_filtered.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    kept.write(out)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare official dylkot/cNMF runs by sample.")
    parser.add_argument("--sample", default="all", help="One sample_id or 'all'.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    table = load_sample_table()
    if args.sample != "all":
        table = table[table["sample_id"] == args.sample].copy()
    if table.empty:
        raise SystemExit("No eligible sample found.")

    rows = []
    for row in table.itertuples(index=False):
        rows.append(
            run_prepare(
                sample_id=str(row.sample_id),
                counts_h5ad=str(row.output_h5ad),
                k_grid=list(row.k_grid),
            )
        )

    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "prepare_summary.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
