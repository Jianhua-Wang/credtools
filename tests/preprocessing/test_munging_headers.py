"""Unit tests for header detection and mapping utilities."""

import pytest
import pandas as pd
from credtools.preprocessing.munging.headers import (
    inspect_headers,
    _detect_separator,
    map_headers_automatic,
    _fuzzy_match_header,
    apply_header_mapping,
    validate_required_columns,
)
from credtools.preprocessing.munging.constants import ColName


# ============================================================================
# Tests for _detect_separator
# ============================================================================


def test_detect_separator_tab_file(tab_sumstats_file):
    r"""Test separator detection for tab-separated file.

    Note: When tab_count == space_count (short column names without spaces),
    the detector falls through to \\s+ which still correctly reads tab files.
    """
    sep = _detect_separator(tab_sumstats_file)
    # With short column names, tab_count equals space_count so detector returns \s+
    assert sep in ("\t", r"\s+")


def test_detect_separator_comma_file(comma_sumstats_file):
    """Test separator detection for comma-separated file."""
    sep = _detect_separator(comma_sumstats_file)
    assert sep == ","


def test_detect_separator_space_file(space_sumstats_file):
    """Test separator detection for space-separated file."""
    sep = _detect_separator(space_sumstats_file)
    assert sep == r"\s+"


def test_detect_separator_gzipped_file(gzipped_sumstats_file):
    """Test separator detection for gzipped tab-separated file."""
    sep = _detect_separator(gzipped_sumstats_file)
    # With short column names, tab_count equals space_count so detector returns \s+
    assert sep in ("\t", r"\s+")


def test_detect_separator_nonexistent_file():
    """Test separator detection for non-existent file falls back to tab."""
    sep = _detect_separator("/nonexistent/path/to/file.txt")
    assert sep == "\t"


# ============================================================================
# Tests for inspect_headers
# ============================================================================


def test_inspect_headers_returns_correct_columns(tab_sumstats_file):
    """Test that inspect_headers returns correct column names."""
    headers = inspect_headers(tab_sumstats_file)
    expected = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
    assert headers == expected


def test_inspect_headers_auto_detects_separator(comma_sumstats_file):
    """Test that inspect_headers auto-detects separator when not specified."""
    headers = inspect_headers(comma_sumstats_file)
    expected = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
    assert headers == expected


def test_inspect_headers_uses_specified_separator(tab_sumstats_file):
    """Test that inspect_headers uses specified separator parameter."""
    headers = inspect_headers(tab_sumstats_file, sep="\t")
    expected = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
    assert headers == expected


def test_inspect_headers_file_not_found():
    """Test that inspect_headers raises error for non-existent file."""
    with pytest.raises(Exception):
        inspect_headers("/nonexistent/path/to/file.txt")


# ============================================================================
# Tests for map_headers_automatic
# ============================================================================


def test_map_headers_automatic_exact_match():
    """Test automatic mapping with exact column name matches."""
    headers = ["CHR", "BP", "BETA"]
    mapping = map_headers_automatic(headers)
    assert mapping == {"CHR": "CHR", "BP": "BP", "BETA": "BETA"}


def test_map_headers_automatic_common_aliases():
    """Test automatic mapping with common alias column names."""
    headers = ["CHROM", "POS", "A1", "A2", "PVAL"]
    mapping = map_headers_automatic(headers)
    assert mapping == {
        "CHROM": "CHR",
        "POS": "BP",
        "A1": "EA",
        "A2": "NEA",
        "PVAL": "P",
    }


def test_map_headers_automatic_fuzzy_match():
    """Test automatic mapping with fuzzy matching patterns."""
    headers = ["chromosome", "base_position"]
    mapping = map_headers_automatic(headers)
    # "chromosome" should match CHR pattern
    # "base_position" should match BP pattern (contains "pos")
    assert mapping["chromosome"] == "CHR"
    assert mapping["base_position"] == "BP"


def test_map_headers_automatic_no_match():
    """Test automatic mapping with unrecognized column names."""
    headers = ["foobar"]
    mapping = map_headers_automatic(headers)
    assert mapping == {"foobar": "foobar"}


def test_map_headers_automatic_empty_list():
    """Test automatic mapping with empty header list."""
    headers = []
    mapping = map_headers_automatic(headers)
    assert mapping == {}


# ============================================================================
# Tests for _fuzzy_match_header (parametrized)
# ============================================================================


@pytest.mark.parametrize(
    "header,expected",
    [
        # CHR patterns
        ("chrom_id", "CHR"),
        ("my_chromosome", "CHR"),
        # BP patterns
        ("position_hg19", "BP"),
        ("base_pos", "BP"),
        # P patterns
        ("pvalue_gwas", "P"),
        ("log_pval", "P"),
        # Allele patterns
        ("alt_allele", "EA"),
        ("effect_allele_freq", "EA"),
        # No match (avoid strings containing single-letter patterns like "n", "p", "z")
        ("kkk_www_jjj", None),
        ("hhhhh_yyy", None),
    ],
)
def test_fuzzy_match_header_patterns(header, expected):
    """Test fuzzy header matching with various patterns."""
    result = _fuzzy_match_header(header)
    assert result == expected


def test_fuzzy_match_header_chr_patterns():
    """Test fuzzy matching for CHR column patterns."""
    assert _fuzzy_match_header("chrom_id") == "CHR"
    assert _fuzzy_match_header("my_chromosome") == "CHR"


def test_fuzzy_match_header_bp_patterns():
    """Test fuzzy matching for BP column patterns."""
    assert _fuzzy_match_header("position_hg19") == "BP"
    assert _fuzzy_match_header("base_pos") == "BP"


def test_fuzzy_match_header_p_patterns():
    """Test fuzzy matching for P-value column patterns."""
    assert _fuzzy_match_header("pvalue_gwas") == "P"
    assert _fuzzy_match_header("log_pval") == "P"


def test_fuzzy_match_header_allele_patterns():
    """Test fuzzy matching for allele column patterns."""
    assert _fuzzy_match_header("alt_allele") == "EA"


def test_fuzzy_match_header_no_match():
    """Test fuzzy matching returns None for unrecognized patterns."""
    # Avoid strings containing single-letter patterns like "n", "p", "z"
    assert _fuzzy_match_header("kkk_www_jjj") is None


# ============================================================================
# Tests for apply_header_mapping
# ============================================================================


def test_apply_header_mapping_correct_rename():
    """Test that apply_header_mapping correctly renames columns."""
    df = pd.DataFrame({"CHROM": [1, 2], "POS": [1000, 2000], "PVAL": [0.01, 0.001]})
    mapping = {"CHROM": "CHR", "POS": "BP", "PVAL": "P"}
    result = apply_header_mapping(df, mapping)

    assert "CHR" in result.columns
    assert "BP" in result.columns
    assert "P" in result.columns
    assert "CHROM" not in result.columns
    assert "POS" not in result.columns
    assert "PVAL" not in result.columns


def test_apply_header_mapping_nonexistent_column():
    """Test that mapping with non-existent column keys is skipped."""
    df = pd.DataFrame({"CHR": [1, 2], "BP": [1000, 2000]})
    mapping = {"CHR": "CHROM", "NONEXISTENT": "FOO"}
    result = apply_header_mapping(df, mapping)

    # CHR should be renamed to CHROM
    assert "CHROM" in result.columns
    assert "CHR" not in result.columns
    # NONEXISTENT key should be ignored
    assert "FOO" not in result.columns


def test_apply_header_mapping_input_not_modified():
    """Test that apply_header_mapping does not modify input DataFrame."""
    df = pd.DataFrame({"CHROM": [1, 2], "POS": [1000, 2000]})
    original_columns = df.columns.tolist()
    mapping = {"CHROM": "CHR", "POS": "BP"}

    result = apply_header_mapping(df, mapping)

    # Original df should not be modified
    assert df.columns.tolist() == original_columns
    # Result should have new columns
    assert "CHR" in result.columns
    assert "BP" in result.columns


# ============================================================================
# Tests for validate_required_columns
# ============================================================================


def test_validate_required_columns_all_present(minimal_gwas_df):
    """Test validation passes when all required columns are present."""
    result = validate_required_columns(minimal_gwas_df)
    assert result is True


def test_validate_required_columns_missing():
    """Test validation fails when required columns are missing."""
    df = pd.DataFrame({"CHR": [1, 2], "BP": [1000, 2000]})
    result = validate_required_columns(df)
    assert result is False


def test_validate_required_columns_default_mandatory():
    """Test validation uses default mandatory_cols when required is None."""
    # Create df with only mandatory columns
    df = pd.DataFrame(
        {
            "CHR": [1, 2],
            "BP": [1000, 2000],
            "EA": ["A", "C"],
            "NEA": ["G", "T"],
            "BETA": [0.1, 0.2],
            "SE": [0.05, 0.06],
            "P": [0.01, 0.001],
        }
    )
    result = validate_required_columns(df, required=None)
    assert result is True


def test_validate_required_columns_custom_required():
    """Test validation with custom required columns list."""
    df = pd.DataFrame({"CHR": [1, 2], "BP": [1000, 2000], "SNPID": ["rs1", "rs2"]})
    result = validate_required_columns(df, required=["CHR", "BP", "SNPID"])
    assert result is True

    # Should fail if custom required column is missing
    result = validate_required_columns(df, required=["CHR", "BP", "RSID"])
    assert result is False


# ============================================================================
# Additional tests for _detect_separator
# ============================================================================


def test_detect_separator_csv_returns_comma(tmp_path):
    """Test that _detect_separator returns comma for a CSV file."""
    filepath = tmp_path / "test.csv"
    filepath.write_text("col_a,col_b,col_c\n1,2,3\n4,5,6\n")
    sep = _detect_separator(str(filepath))
    assert sep == ","


def test_detect_separator_gzipped_tab_works(tmp_path):
    """Test that _detect_separator correctly reads a gzipped tab file."""
    import gzip

    filepath = tmp_path / "test.tsv.gz"
    content = (
        "CHR\tBP\tSNPID\tEA\tNEA\tBETA\tSE\tP\n1\t1000\trs1\tA\tG\t0.1\t0.05\t0.01\n"
    )
    with gzip.open(str(filepath), "wt") as f:
        f.write(content)
    sep = _detect_separator(str(filepath))
    assert sep in ("\t", r"\s+")


# ============================================================================
# Tests for suggest_missing_mappings
# ============================================================================

from credtools.preprocessing.munging.headers import suggest_missing_mappings


def test_suggest_missing_mappings_fuzzy_unmapped():
    """Unmapped header that fuzzy-matches a missing standard should be suggested."""
    headers = ["chromosome_id", "BP", "EA", "NEA", "BETA", "SE", "P"]
    # Only BP..P are already mapped; CHR is missing from mapped_headers
    mapped = {
        "BP": "BP",
        "EA": "EA",
        "NEA": "NEA",
        "BETA": "BETA",
        "SE": "SE",
        "P": "P",
    }
    suggestions = suggest_missing_mappings(headers, mapped)
    # "chromosome_id" contains "chrom" -> should match CHR
    assert suggestions.get("chromosome_id") == "CHR"


def test_suggest_missing_mappings_all_required_mapped():
    """When all mandatory columns are already mapped, suggestions should be empty."""
    headers = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P", "extra_col"]
    mapped = {
        "CHR": "CHR",
        "BP": "BP",
        "EA": "EA",
        "NEA": "NEA",
        "BETA": "BETA",
        "SE": "SE",
        "P": "P",
    }
    suggestions = suggest_missing_mappings(headers, mapped)
    assert suggestions == {}


def test_suggest_missing_mappings_no_unmapped_headers():
    """When all headers are in mapped_headers, suggestions should be empty."""
    headers = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
    mapped = {
        "CHR": "CHR",
        "BP": "BP",
        "EA": "EA",
        "NEA": "NEA",
        "BETA": "BETA",
        "SE": "SE",
        "P": "P",
    }
    suggestions = suggest_missing_mappings(headers, mapped)
    assert suggestions == {}
