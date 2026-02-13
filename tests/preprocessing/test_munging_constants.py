"""
Unit tests for munging constants module.

Tests the column definitions, data types, validation ranges, and column
mapping suggestion functionality used in GWAS summary statistics munging.
"""

import pytest
import numpy as np

from credtools.preprocessing.munging.constants import (
    ColName,
    ColType,
    ColRange,
    ColAllowNA,
    COMMON_COLNAMES,
    CHROM_LENGTHS,
    suggest_column_mapping,
)


class TestColNameStructure:
    """Tests for ColName class structure and consistency."""

    def test_mandatory_cols_subset_of_sumstat_cols(self):
        """Test that mandatory_cols is a proper subset of sumstat_cols.

        All mandatory columns must be present in the full set of summary
        statistic columns.
        """
        assert set(ColName.mandatory_cols).issubset(set(ColName.sumstat_cols))

        # Also verify mandatory columns are actually mandatory
        expected_mandatory = {
            ColName.CHR,
            ColName.BP,
            ColName.EA,
            ColName.NEA,
            ColName.BETA,
            ColName.SE,
            ColName.P,
        }
        assert set(ColName.mandatory_cols) == expected_mandatory

    def test_output_cols_has_11_columns(self):
        """Test that output_cols contains exactly 11 specific columns.

        The output columns for credtools munge command must include:
        CHR, BP, SNPID, EA, NEA, EAF, BETA, SE, P, N, RSID
        """
        assert len(ColName.output_cols) == 11

        expected_output_cols = [
            ColName.CHR,
            ColName.BP,
            ColName.SNPID,
            ColName.EA,
            ColName.NEA,
            ColName.EAF,
            ColName.BETA,
            ColName.SE,
            ColName.P,
            ColName.N,
            ColName.RSID,
        ]
        assert ColName.output_cols == expected_output_cols

        # Verify all output columns are in sumstat_cols
        assert set(ColName.output_cols).issubset(set(ColName.sumstat_cols))


class TestColRange:
    """Tests for ColRange validation constraints."""

    def test_colrange_min_less_than_max(self):
        """Test that all MIN values are less than corresponding MAX values.

        This ensures range validation constraints are logically consistent.
        """
        # Test all MIN/MAX pairs
        assert ColRange.CHR_MIN < ColRange.CHR_MAX
        assert ColRange.BP_MIN < ColRange.BP_MAX
        assert ColRange.P_MIN < ColRange.P_MAX
        assert ColRange.EAF_MIN < ColRange.EAF_MAX
        assert ColRange.MAF_MIN < ColRange.MAF_MAX
        assert ColRange.INFO_MIN < ColRange.INFO_MAX

        # Test specific expected values
        assert ColRange.CHR_MIN == 1
        assert ColRange.CHR_MAX == 23  # Including X chromosome
        assert ColRange.BP_MIN == 0
        assert ColRange.BP_MAX == 300_000_000
        assert ColRange.P_MIN == 0.0
        assert ColRange.P_MAX == 1.0
        assert ColRange.EAF_MIN == 0.0
        assert ColRange.EAF_MAX == 1.0
        assert ColRange.MAF_MIN == 0.0
        assert ColRange.MAF_MAX == 0.5
        assert ColRange.INFO_MIN == 0.0
        assert ColRange.INFO_MAX == 1.0

        # Test minimum-only values are positive
        assert ColRange.SE_MIN == 0.0
        assert ColRange.OR_MIN == 1e-10
        assert ColRange.OR_MIN > 0
        assert ColRange.N_MIN == 1
        assert ColRange.N_MIN > 0


class TestCommonColnames:
    """Tests for COMMON_COLNAMES mapping dictionary."""

    def test_common_colnames_values_are_valid(self):
        """Test that all COMMON_COLNAMES values are valid ColName attributes.

        All mapped values must be present in sumstat_cols to ensure consistent
        column naming throughout the pipeline.
        """
        for input_name, mapped_name in COMMON_COLNAMES.items():
            assert (
                mapped_name in ColName.sumstat_cols
            ), f"{input_name} maps to {mapped_name}, which is not in sumstat_cols"

        # Verify some key mappings exist
        assert "CHR" in COMMON_COLNAMES
        assert "BP" in COMMON_COLNAMES
        assert "BETA" in COMMON_COLNAMES
        assert "P" in COMMON_COLNAMES
        assert "EA" in COMMON_COLNAMES
        assert "NEA" in COMMON_COLNAMES

        # Verify alternative names map correctly
        assert COMMON_COLNAMES["CHROM"] == ColName.CHR
        assert COMMON_COLNAMES["POS"] == ColName.BP
        assert COMMON_COLNAMES["A1"] == ColName.EA
        assert COMMON_COLNAMES["A2"] == ColName.NEA
        assert COMMON_COLNAMES["PVAL"] == ColName.P


class TestChromLengths:
    """Tests for CHROM_LENGTHS chromosome size data."""

    def test_chrom_lengths_has_23_entries(self):
        """Test that CHROM_LENGTHS contains all 23 chromosomes (1-22 + X).

        Keys should be integers 1-23 (23 represents X chromosome).
        All values should be positive integers representing base pair lengths.
        """
        assert len(CHROM_LENGTHS) == 23

        # Verify all chromosomes 1-23 are present
        expected_chroms = set(range(1, 24))
        assert set(CHROM_LENGTHS.keys()) == expected_chroms

        # Verify all lengths are positive integers
        for chrom, length in CHROM_LENGTHS.items():
            assert isinstance(chrom, int)
            assert isinstance(length, int)
            assert length > 0

        # Verify chromosome lengths are reasonable (GRCh37/hg19)
        # Chromosome 1 is the longest
        assert CHROM_LENGTHS[1] == 249250621
        # Chromosome 21 is one of the shortest autosomes
        assert CHROM_LENGTHS[21] == 46709983
        # X chromosome (23) should be reasonable
        assert CHROM_LENGTHS[23] == 156040895

        # All lengths should be less than 300Mb (ColRange.BP_MAX)
        for length in CHROM_LENGTHS.values():
            assert length < ColRange.BP_MAX


class TestSuggestColumnMapping:
    """Tests for suggest_column_mapping function."""

    def test_suggest_column_mapping_exact_match(self):
        """Test exact matches produce high confidence suggestions.

        When column names exactly match COMMON_COLNAMES keys, the function
        should return high confidence mappings.
        """
        headers = ["CHR", "BP", "BETA", "P", "EA", "NEA"]
        suggestions = suggest_column_mapping(headers)

        assert len(suggestions) == len(headers)

        for header in headers:
            assert header in suggestions
            assert suggestions[header]["confidence"] == "high"
            assert suggestions[header]["suggested"] in ColName.sumstat_cols

        # Verify specific mappings
        assert suggestions["CHR"]["suggested"] == ColName.CHR
        assert suggestions["BP"]["suggested"] == ColName.BP
        assert suggestions["BETA"]["suggested"] == ColName.BETA
        assert suggestions["P"]["suggested"] == ColName.P
        assert suggestions["EA"]["suggested"] == ColName.EA
        assert suggestions["NEA"]["suggested"] == ColName.NEA

    def test_suggest_column_mapping_fuzzy_match(self):
        """Test fuzzy pattern matches produce medium confidence suggestions.

        When column names contain common patterns (chr, pos, beta, pval)
        but don't exactly match COMMON_COLNAMES, the function should use
        fuzzy matching with medium confidence.
        """
        headers = ["my_chr", "position_bp", "effect_beta", "p_val_adjusted"]
        suggestions = suggest_column_mapping(headers)

        assert len(suggestions) == len(headers)

        # Test chromosome fuzzy match
        assert suggestions["my_chr"]["confidence"] == "medium"
        assert suggestions["my_chr"]["suggested"] == ColName.CHR

        # Test position fuzzy match
        assert suggestions["position_bp"]["confidence"] == "medium"
        assert suggestions["position_bp"]["suggested"] == ColName.BP

        # Test beta fuzzy match
        assert suggestions["effect_beta"]["confidence"] == "medium"
        assert suggestions["effect_beta"]["suggested"] == ColName.BETA

        # Test p-value fuzzy match
        assert suggestions["p_val_adjusted"]["confidence"] == "medium"
        assert suggestions["p_val_adjusted"]["suggested"] == ColName.P

        # Test additional patterns
        headers2 = ["chrom_num", "pos_hg19", "effect_size", "pvalue_meta"]
        suggestions2 = suggest_column_mapping(headers2)

        assert suggestions2["chrom_num"]["confidence"] == "medium"
        assert suggestions2["chrom_num"]["suggested"] == ColName.CHR
        assert suggestions2["pos_hg19"]["confidence"] == "medium"
        assert suggestions2["pos_hg19"]["suggested"] == ColName.BP
        assert suggestions2["effect_size"]["confidence"] == "medium"
        assert suggestions2["effect_size"]["suggested"] == ColName.BETA
        assert suggestions2["pvalue_meta"]["confidence"] == "medium"
        assert suggestions2["pvalue_meta"]["suggested"] == ColName.P

    def test_suggest_column_mapping_no_match(self):
        """Test unrecognized columns produce low confidence suggestions.

        When column names don't match any known patterns, the function
        should return the original name with low confidence.
        """
        headers = ["unknown_col", "random_field", "mystery_data", "xyz123"]
        suggestions = suggest_column_mapping(headers)

        assert len(suggestions) == len(headers)

        for header in headers:
            assert header in suggestions
            assert suggestions[header]["confidence"] == "low"
            # Should suggest keeping the original name
            assert suggestions[header]["suggested"] == header

    def test_suggest_column_mapping_mixed_confidence(self):
        """Test mixed input with high, medium, and low confidence matches.

        Real-world headers often contain a mixture of exact matches,
        fuzzy matches, and unrecognized columns.
        """
        headers = [
            "CHR",  # exact match -> high
            "my_position",  # fuzzy match -> medium
            "unknown_stat",  # no match -> low
            "BETA",  # exact match -> high
            "pval_gwas",  # fuzzy match -> medium
        ]
        suggestions = suggest_column_mapping(headers)

        assert len(suggestions) == len(headers)
        assert suggestions["CHR"]["confidence"] == "high"
        assert suggestions["CHR"]["suggested"] == ColName.CHR
        assert suggestions["my_position"]["confidence"] == "medium"
        assert suggestions["my_position"]["suggested"] == ColName.BP
        assert suggestions["unknown_stat"]["confidence"] == "low"
        assert suggestions["unknown_stat"]["suggested"] == "unknown_stat"
        assert suggestions["BETA"]["confidence"] == "high"
        assert suggestions["BETA"]["suggested"] == ColName.BETA
        assert suggestions["pval_gwas"]["confidence"] == "medium"
        assert suggestions["pval_gwas"]["suggested"] == ColName.P


class TestColTypeDataTypes:
    """Tests for ColType numpy data type specifications."""

    def test_coltype_numeric_types_are_numpy(self):
        """Test that numeric column types are numpy dtypes.

        Ensures proper type specifications for efficient pandas operations.
        """
        # Integer types
        assert ColType.CHR == np.int8
        assert ColType.BP == np.int32
        assert ColType.N == np.int32

        # Float types
        assert ColType.EAF == np.float32
        assert ColType.BETA == np.float32
        assert ColType.SE == np.float32
        assert ColType.P == np.float64  # Higher precision for p-values
        assert ColType.MAF == np.float32
        assert ColType.INFO == np.float32
        assert ColType.Z == np.float32
        assert ColType.OR == np.float32
        assert ColType.OR_SE == np.float32
        assert ColType.NEGLOG10P == np.float32

        # String types
        assert ColType.SNPID == str
        assert ColType.EA == str
        assert ColType.NEA == str
        assert ColType.RSID == str


class TestColAllowNA:
    """Tests for ColAllowNA NA value permission flags."""

    def test_mandatory_cols_do_not_allow_na(self):
        """Test that mandatory columns do not allow NA values.

        Core required columns (CHR, BP, EA, NEA, BETA, SE, P) should never
        contain missing values.
        """
        assert ColAllowNA.CHR is False
        assert ColAllowNA.BP is False
        assert ColAllowNA.EA is False
        assert ColAllowNA.NEA is False
        assert ColAllowNA.BETA is False
        assert ColAllowNA.SE is False
        assert ColAllowNA.P is False

        # SNPID is also required for identification
        assert ColAllowNA.SNPID is False

    def test_optional_cols_allow_na(self):
        """Test that optional columns allow NA values.

        Optional columns (EAF, RSID, MAF, N, INFO, Z, OR, etc.) can contain
        missing values as they may not be available in all datasets.
        """
        assert ColAllowNA.EAF is True
        assert ColAllowNA.RSID is True
        assert ColAllowNA.MAF is True
        assert ColAllowNA.N is True
        assert ColAllowNA.INFO is True
        assert ColAllowNA.Z is True
        assert ColAllowNA.OR is True
        assert ColAllowNA.OR_SE is True
        assert ColAllowNA.NEGLOG10P is True
