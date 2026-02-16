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


# ---------- chunk-related fixtures ----------


def _make_munged_df(chr_list, n_per_chr, sig_indices=None, base_p=0.5):
    """Helper to build a munged-style DataFrame.

    Parameters
    ----------
    chr_list : list[int]
        Chromosomes to include.
    n_per_chr : int
        Number of SNPs per chromosome.
    sig_indices : dict[int, list[int]] | None
        Mapping of chr -> list of intra-chr indices that should be genome-wide
        significant (P = 1e-9).  If *None*, the first SNP of each chr is
        significant.
    base_p : float
        Default (non-significant) p-value.
    """
    rows = []
    if sig_indices is None:
        sig_indices = {c: [0] for c in chr_list}
    for c in chr_list:
        for j in range(n_per_chr):
            bp = 100_000 + j * 10_000
            p = 1e-9 if j in sig_indices.get(c, []) else base_p
            rows.append(
                {
                    "CHR": c,
                    "BP": bp,
                    "SNPID": f"{c}-{bp}-A-G",
                    "EA": "A",
                    "NEA": "G",
                    "EAF": 0.3,
                    "BETA": 0.1,
                    "SE": 0.05,
                    "P": p,
                    "N": 10000,
                    "RSID": f"rs{c}_{j}",
                    "MAF": 0.3,
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture
def munged_sumstats_df():
    """Munged DataFrame: 3 chromosomes × 20 SNPs, first SNP per chr significant."""
    return _make_munged_df([1, 2, 3], 20)


@pytest.fixture
def munged_sumstats_no_sig_df():
    """Munged DataFrame where all P > 5e-8 (no genome-wide significant SNPs)."""
    return _make_munged_df([1, 2], 10, sig_indices={}, base_p=0.1)


@pytest.fixture
def munged_sumstats_gz_file(tmp_path, munged_sumstats_df):
    """Write munged_sumstats_df to a gzip TSV file."""
    filepath = tmp_path / "EUR.munged.txt.gz"
    munged_sumstats_df.to_csv(filepath, sep="\t", index=False, compression="gzip")
    return str(filepath)


@pytest.fixture
def two_ancestry_gz_files(tmp_path):
    """EUR + ASN gzip files with partially overlapping loci.

    EUR chr1 significant at index 0 (BP=100000), chr2 significant at index 0.
    ASN chr1 significant at index 0 (same region → overlap), chr3 significant at index 0.
    """
    eur_df = _make_munged_df([1, 2], 20, sig_indices={1: [0], 2: [0]})
    asn_df = _make_munged_df([1, 3], 20, sig_indices={1: [0], 3: [0]})

    eur_path = tmp_path / "EUR.munged.txt.gz"
    asn_path = tmp_path / "ASN.munged.txt.gz"

    eur_df.to_csv(eur_path, sep="\t", index=False, compression="gzip")
    asn_df.to_csv(asn_path, sep="\t", index=False, compression="gzip")

    return {"EUR": str(eur_path), "ASN": str(asn_path)}


@pytest.fixture
def sample_loci_df():
    """Simulated output of identify_independent_loci."""
    return pd.DataFrame(
        [
            {
                "chr": 1,
                "start": 1,
                "end": 350000,
                "lead_snp": "1-100000-A-G",
                "lead_bp": 100000,
                "lead_p": 1e-9,
                "ancestry": "EUR",
                "n_variants": 20,
                "locus_id": "chr1_1_350000",
            },
            {
                "chr": 2,
                "start": 1,
                "end": 350000,
                "lead_snp": "2-100000-A-G",
                "lead_bp": 100000,
                "lead_p": 1e-9,
                "ancestry": "EUR",
                "n_variants": 20,
                "locus_id": "chr2_1_350000",
            },
        ]
    )


@pytest.fixture
def sample_chunk_info_df(tmp_path):
    """Simulated output of chunk_sumstats."""
    return pd.DataFrame(
        [
            {
                "locus_id": "chr1_1_350000",
                "ancestry": "EUR",
                "chr": 1,
                "start": 1,
                "end": 350000,
                "n_variants": 20,
                "sumstats_file": str(
                    tmp_path / "chunks" / "EUR.chr1_1_350000.sumstats.gz"
                ),
            },
            {
                "locus_id": "chr2_1_350000",
                "ancestry": "EUR",
                "chr": 2,
                "start": 1,
                "end": 350000,
                "n_variants": 20,
                "sumstats_file": str(
                    tmp_path / "chunks" / "EUR.chr2_1_350000.sumstats.gz"
                ),
            },
        ]
    )
