#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path("/Users/niubi/Desktop/SACC")
ANNOTATED_H5AD = (
    ROOT
    / "analysis"
    / "gse294017_cleanup"
    / "output"
    / "paper_85142_annotated"
    / "GSE294017_paper_85142.annotated.h5ad"
)
RAW_H5_DIR = ROOT / "GSE294017" / "raw_h5"
TEST_H5_DIR = ROOT / "GSE294017" / "raw_h5_test"
OUTDIR = (
    ROOT
    / "analysis"
    / "gse294017_cleanup"
    / "output"
    / "malignant_counts_by_sample"
)
CELL_LABEL = "Malignant_cells"
MIN_MALIGNANT_CELLS = 200


def barcode_from_obs_name(obs_name: str) -> str:
    return str(obs_name).rsplit("_", 1)[-1]


def resolve_h5(sample_id: str) -> Path:
    suffix = sample_id.split("GSM", 1)[-1]
    gsm = f"GSM{suffix.split('_', 1)[0]}"
    patient = suffix.split("_", 1)[1]
    filename = f"{gsm}_{patient}_output_filtered.h5"
    for base in (RAW_H5_DIR, TEST_H5_DIR):
        candidate = base / filename
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Missing raw h5 for {sample_id}: {filename}")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(ANNOTATED_H5AD)
    malignant = adata[adata.obs["cluster_annotation"].astype(str) == CELL_LABEL].copy()
    malignant.obs["cell_barcode"] = [barcode_from_obs_name(x) for x in malignant.obs_names]

    summary_rows = []
    skipped_rows = []

    for sample_id, sample_obs in malignant.obs.groupby("sample_id", observed=True):
        h5_path = resolve_h5(str(sample_id))
        sample_barcodes = sorted(sample_obs["cell_barcode"].astype(str).unique())
        if len(sample_barcodes) < MIN_MALIGNANT_CELLS:
            skipped_rows.append(
                {
                    "sample_id": str(sample_id),
                    "n_malignant_cells": len(sample_barcodes),
                    "reason": f"below_min_malignant_cells_{MIN_MALIGNANT_CELLS}",
                }
            )
            continue

        raw = sc.read_10x_h5(h5_path)
        raw.var_names_make_unique()
        common = sorted(set(raw.obs_names).intersection(sample_barcodes))
        if len(common) == 0:
            skipped_rows.append(
                {
                    "sample_id": str(sample_id),
                    "n_malignant_cells": len(sample_barcodes),
                    "reason": "no_barcode_overlap_with_raw_h5",
                }
            )
            continue

        out = raw[common].copy()
        meta = sample_obs.set_index("cell_barcode").loc[common].copy()
        out.obs = meta
        out.obs["sample_id"] = str(sample_id)
        out.obs["raw_h5_path"] = str(h5_path)

        out_path = OUTDIR / f"{sample_id}.malignant_counts.h5ad"
        out.write(out_path)
        summary_rows.append(
            {
                "sample_id": str(sample_id),
                "raw_h5_path": str(h5_path),
                "n_malignant_cells_requested": len(sample_barcodes),
                "n_malignant_cells_written": out.n_obs,
                "n_genes": out.n_vars,
                "output_h5ad": str(out_path),
            }
        )

    pd.DataFrame(summary_rows).sort_values(
        "n_malignant_cells_written", ascending=False
    ).to_csv(OUTDIR / "extracted_samples.tsv", sep="\t", index=False)
    pd.DataFrame(skipped_rows).sort_values(
        "n_malignant_cells", ascending=False
    ).to_csv(OUTDIR / "skipped_samples.tsv", sep="\t", index=False)

    print(pd.DataFrame(summary_rows).to_string(index=False))
    if skipped_rows:
        print("\nSkipped:")
        print(pd.DataFrame(skipped_rows).to_string(index=False))


if __name__ == "__main__":
    main()
