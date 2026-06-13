#!/Users/niubi/Desktop/SACC/.venv-scanpy/bin/python
from __future__ import annotations

from pathlib import Path
import shutil

import pandas as pd


ROOT = Path("/Users/niubi/Desktop/SACC")
ANALYSIS_DIR = ROOT / "analysis" / "gse294017_cleanup"
OUTPUT_DIR = ANALYSIS_DIR / "output"
REPORTS_DIR = ANALYSIS_DIR / "reports"
CNMF_DIR = OUTPUT_DIR / "malignant_cnmf_per_sample"
COUNTS_DIR = OUTPUT_DIR / "malignant_counts_by_sample"


def remove_junk_files() -> list[str]:
    removed: list[str] = []
    patterns = [".DS_Store", "*.pyc"]
    for pattern in patterns:
        for path in ROOT.rglob(pattern):
            if path.is_file():
                path.unlink()
                removed.append(str(path))
    for path in ROOT.rglob("__pycache__"):
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(str(path))
    return removed


def remove_retry_intermediates() -> list[str]:
    removed: list[str] = []
    if not CNMF_DIR.exists():
        return removed

    for filtered_h5ad in CNMF_DIR.glob("*/*.zero_feature_filtered.h5ad"):
        sample_id = filtered_h5ad.parent.name
        norm_counts = filtered_h5ad.parent / sample_id / "cnmf_tmp" / f"{sample_id}.norm_counts.h5ad"
        if norm_counts.exists():
            filtered_h5ad.unlink()
            removed.append(str(filtered_h5ad))
    return removed


def copy_if_exists(src: Path, dest: Path) -> bool:
    if not src.exists():
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return True


def build_disk_usage_report() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path in [OUTPUT_DIR, COUNTS_DIR, CNMF_DIR]:
        if path.exists():
            size_bytes = sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
            rows.append(
                {
                    "path": str(path),
                    "size_gb": round(size_bytes / (1024 ** 3), 3),
                }
            )
    return pd.DataFrame(rows)


def sync_reports() -> list[str]:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []

    summary_map = {
        COUNTS_DIR / "extracted_samples.tsv": REPORTS_DIR / "malignant_counts_extracted.tsv",
        COUNTS_DIR / "skipped_samples.tsv": REPORTS_DIR / "malignant_counts_skipped.tsv",
        CNMF_DIR / "prepare_summary.tsv": REPORTS_DIR / "cnmf_prepare_summary.tsv",
        CNMF_DIR / "selected_programs.tsv": REPORTS_DIR / "cnmf_selected_programs.tsv",
        ROOT / ".gitignore": REPORTS_DIR / "gitignore_snapshot.txt",
    }
    for src, dest in summary_map.items():
        if copy_if_exists(src, dest):
            copied.append(str(dest))

    disk_usage = build_disk_usage_report()
    if not disk_usage.empty:
        disk_usage.to_csv(REPORTS_DIR / "disk_usage.tsv", sep="\t", index=False)
        copied.append(str(REPORTS_DIR / "disk_usage.tsv"))

    return copied


def main() -> None:
    removed_junk = remove_junk_files()
    removed_retry = remove_retry_intermediates()
    copied = sync_reports()

    print("Removed junk files:", len(removed_junk))
    for path in removed_junk[:20]:
        print(path)

    print("\nRemoved retry intermediates:", len(removed_retry))
    for path in removed_retry[:20]:
        print(path)

    print("\nSynced report files:", len(copied))
    for path in copied:
        print(path)


if __name__ == "__main__":
    main()
