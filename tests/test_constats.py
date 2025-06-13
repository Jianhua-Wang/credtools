"""Tests for the constants module."""

from credtools.constants import CHROM_LIST, ColName


def test_colname_mandatory_columns():
    """
    Test that the mandatory columns are correctly defined.

    Returns
    -------
    None
    """
    expected_cols = ["CHR", "BP", "EA", "NEA", "EAF", "BETA", "SE", "P"]
    assert ColName.mandatory_cols == expected_cols


def test_colname_sumstat_columns():
    """
    Test that the sumstat columns are correctly defined.

    Returns
    -------
    None
    """
    expected_cols = [
        "SNPID",
        "CHR",
        "BP",
        "rsID",
        "EA",
        "NEA",
        "EAF",
        "MAF",
        "BETA",
        "SE",
        "P",
    ]
    assert ColName.sumstat_cols == expected_cols


def test_colname_loci_columns():
    """
    Test that the loci columns are correctly defined.

    Returns
    -------
    None
    """
    expected_cols = ["CHR", "START", "END", "LEAD_SNP", "LEAD_SNP_P", "LEAD_SNP_BP"]
    assert ColName.loci_cols == expected_cols


def test_chrom_list():
    """
    Test that the chromosome list is correctly defined.

    Returns
    -------
    None
    """
    expected_chrom_list = [i for i in range(1, 24)]
    assert CHROM_LIST == expected_chrom_list
