"""Shared fixtures for preprocessing unit tests."""

import gzip
import json

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def minimal_gwas_df():
    """Minimal valid GWAS DataFrame with 7 mandatory columns."""
    return pd.DataFrame(
        {
            "CHR": [1, 1, 2, 2, 3],
            "BP": [1000, 2000, 3000, 4000, 5000],
            "EA": ["A", "C", "G", "T", "A"],
            "NEA": ["G", "T", "C", "A", "C"],
            "BETA": [0.1, -0.2, 0.3, -0.4, 0.5],
            "SE": [0.05, 0.06, 0.07, 0.08, 0.09],
            "P": [0.01, 0.001, 0.0001, 0.00001, 0.05],
        }
    )


@pytest.fixture
def full_gwas_df(minimal_gwas_df):
    """Full GWAS DataFrame with optional columns."""
    df = minimal_gwas_df.copy()
    df["SNPID"] = ["1-1000-A-G", "1-2000-C-T", "2-3000-C-G", "2-4000-A-T", "3-5000-A-C"]
    df["EAF"] = [0.3, 0.45, 0.1, 0.5, 0.25]
    df["RSID"] = ["rs1", "rs2", "rs3", "rs4", "rs5"]
    df["N"] = [10000, 10000, 10000, 10000, 10000]
    df["MAF"] = [0.3, 0.45, 0.1, 0.5, 0.25]
    return df


@pytest.fixture
def tab_sumstats_file(tmp_path, minimal_gwas_df):
    """Tab-separated sumstats file."""
    filepath = tmp_path / "sumstats.tsv"
    minimal_gwas_df.to_csv(filepath, sep="\t", index=False)
    return str(filepath)


@pytest.fixture
def comma_sumstats_file(tmp_path, minimal_gwas_df):
    """Comma-separated sumstats file."""
    filepath = tmp_path / "sumstats.csv"
    minimal_gwas_df.to_csv(filepath, sep=",", index=False)
    return str(filepath)


@pytest.fixture
def gzipped_sumstats_file(tmp_path, minimal_gwas_df):
    """Gzipped tab-separated sumstats file."""
    filepath = tmp_path / "sumstats.tsv.gz"
    minimal_gwas_df.to_csv(filepath, sep="\t", index=False, compression="gzip")
    return str(filepath)


@pytest.fixture
def space_sumstats_file(tmp_path, minimal_gwas_df):
    """Space-separated sumstats file."""
    filepath = tmp_path / "sumstats.txt"
    minimal_gwas_df.to_csv(filepath, sep=" ", index=False)
    return str(filepath)


@pytest.fixture
def sample_config_dict():
    """Valid config dictionary for munging."""
    return {
        "column_mapping": {
            "CHROM": "CHR",
            "POS": "BP",
            "A1": "EA",
            "A2": "NEA",
            "BETA": "BETA",
            "SE": "SE",
            "PVAL": "P",
        }
    }


@pytest.fixture
def sample_config_file(tmp_path, sample_config_dict):
    """Valid config JSON file."""
    filepath = tmp_path / "config.json"
    with open(filepath, "w") as f:
        json.dump(sample_config_dict, f, indent=2)
    return str(filepath)
