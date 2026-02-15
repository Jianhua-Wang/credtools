"""Shared fixtures for plot unit tests."""

import gzip

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def qc_summary_gz(tmp_path):
    """Generate a gzipped QC summary TSV with enough rows for groupby stats."""
    rows = []
    for locus in ["chr1_100_200", "chr2_300_400", "chr3_500_600"]:
        for popu, cohort in [("EUR", "UKB"), ("EAS", "BBJ"), ("AFR", "APCDR")]:
            rows.append(
                {
                    "locus_id": locus,
                    "popu": popu,
                    "cohort": cohort,
                    "lambda_s": np.random.default_rng(42).uniform(0.8, 1.5),
                    "maf_corr": np.random.default_rng(42).uniform(0.5, 1.0),
                    "n_lambda_s_outlier": np.random.default_rng(42).integers(0, 5),
                    "n_dentist_s_outlier": np.random.default_rng(42).integers(0, 5),
                }
            )
    df = pd.DataFrame(rows)
    filepath = tmp_path / "qc.txt.gz"
    df.to_csv(filepath, sep="\t", index=False, compression="gzip")
    return filepath


@pytest.fixture
def locus_dir_with_pips(tmp_path):
    """Generate a locus directory containing a pips.txt.gz file."""
    locus_dir = tmp_path / "chr1_49782265_50282265"
    locus_dir.mkdir()

    rng = np.random.default_rng(42)
    n = 8
    df = pd.DataFrame(
        {
            "BP": np.arange(49800000, 49800000 + n * 1000, 1000),
            "CRED": [0, 1, 1, 0, 2, 0, 0, 1],
            "PIP": rng.uniform(0, 1, n).round(3),
            "UKB_P": rng.uniform(1e-10, 0.05, n),
            "UKB_R2": rng.uniform(0, 1, n).round(3),
            "BBJ_P": rng.uniform(1e-10, 0.05, n),
            "BBJ_R2": rng.uniform(0, 1, n).round(3),
        }
    )

    pip_path = locus_dir / "pips.txt.gz"
    df.to_csv(pip_path, sep="\t", index=False, compression="gzip")
    return locus_dir
