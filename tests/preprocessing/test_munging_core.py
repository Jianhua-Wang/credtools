"""Unit tests for credtools/preprocessing/munging/core.py."""
import numpy as np
import pandas as pd
import pytest
from credtools.preprocessing.munging.core import (
    _finalize_columns,
    _munge_alleles,
    _munge_beta,
    _munge_bp,
    _munge_chr,
    _munge_eaf,
    _munge_pvalue,
    _munge_se,
    _remove_all_na_columns,
    transform_allele,
    transform_chr,
    make_SNPID_unique,
    munge,
)
from credtools.preprocessing.munging.constants import OUTPUT_COLS, ColName


class TestMunge:
    """Tests for the main munge() function."""

    def test_munge_minimal_valid_df(self, minimal_gwas_df):
        """Minimal valid df returns cleaned df with output_cols."""
        result = munge(minimal_gwas_df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert list(result.columns) == OUTPUT_COLS

    def test_munge_with_eaf_column(self, full_gwas_df):
        """With EAF column, EAF is processed correctly."""
        result = munge(full_gwas_df)
        assert ColName.EAF in result.columns
        assert result[ColName.EAF].notna().any()

    def test_munge_without_eaf_column(self, minimal_gwas_df):
        """Without EAF column, munging still succeeds."""
        result = munge(minimal_gwas_df)
        assert ColName.EAF in result.columns
        assert result[ColName.EAF].isna().all()

    def test_munge_output_columns(self, minimal_gwas_df):
        """Output columns match OUTPUT_COLS in correct order."""
        result = munge(minimal_gwas_df)
        assert list(result.columns) == OUTPUT_COLS
        assert len(result.columns) == 11

    def test_munge_sorted_by_chr_bp(self, minimal_gwas_df):
        """Output is sorted by CHR and BP."""
        df = minimal_gwas_df.copy()
        df.loc[0, "CHR"] = 10
        df.loc[1, "CHR"] = 5
        result = munge(df)
        assert result[ColName.CHR].is_monotonic_increasing or (
            result.groupby(ColName.CHR)[ColName.BP].apply(
                lambda x: x.is_monotonic_increasing
            ).all()
        )

    def test_munge_missing_mandatory_column_raises_error(self, minimal_gwas_df):
        """Missing mandatory column raises ValueError."""
        df = minimal_gwas_df.drop(columns=["BETA"])
        with pytest.raises(ValueError, match="Missing mandatory"):
            munge(df)

    def test_munge_removes_invalid_chr(self, minimal_gwas_df):
        """Invalid CHR values (chr=99) are removed."""
        df = minimal_gwas_df.copy()
        df.loc[0, "CHR"] = 99
        result = munge(df)
        assert 99 not in result[ColName.CHR].values
        assert len(result) == len(df) - 1

    def test_munge_removes_pvalue_zero(self, minimal_gwas_df):
        """P=0 rows are removed (exclude_min)."""
        df = minimal_gwas_df.copy()
        df.loc[0, "P"] = 0
        result = munge(df)
        assert 0 not in result[ColName.P].values
        assert len(result) == len(df) - 1

    def test_munge_removes_pvalue_above_one(self, minimal_gwas_df):
        """P>1.0 rows are removed."""
        df = minimal_gwas_df.copy()
        df.loc[0, "P"] = 2.0
        result = munge(df)
        assert (result[ColName.P] <= 1.0).all()
        assert len(result) == len(df) - 1

    def test_munge_removes_se_zero(self, minimal_gwas_df):
        """SE=0 rows are removed (exclude_min)."""
        df = minimal_gwas_df.copy()
        df.loc[0, "SE"] = 0
        result = munge(df)
        assert 0 not in result[ColName.SE].values
        assert len(result) == len(df) - 1

    def test_munge_removes_identical_alleles(self, minimal_gwas_df):
        """EA==NEA rows are removed."""
        df = minimal_gwas_df.copy()
        df.loc[0, "EA"] = "A"
        df.loc[0, "NEA"] = "A"
        result = munge(df)
        assert len(result) == len(df) - 1

    def test_munge_removes_invalid_alleles(self, minimal_gwas_df):
        """Invalid allele values like '123' are removed."""
        df = minimal_gwas_df.copy()
        df.loc[0, "EA"] = "123"
        result = munge(df)
        assert len(result) == len(df) - 1

    def test_munge_copy_semantics(self, minimal_gwas_df):
        """Input dataframe is not modified (copy semantics)."""
        df = minimal_gwas_df.copy()
        original_columns = df.columns.tolist()
        original_len = len(df)
        munge(df)
        assert df.columns.tolist() == original_columns
        assert len(df) == original_len


class TestMakeSNPIDUnique:
    """Tests for make_SNPID_unique() function."""

    def test_snpid_format(self, minimal_gwas_df):
        """SNPID format is chr-bp-allele1-allele2."""
        result = make_SNPID_unique(minimal_gwas_df, remove_duplicates=False)
        snpid = result[ColName.SNPID].iloc[0]
        parts = snpid.split("-")
        assert len(parts) == 4
        assert parts[0].isdigit()  # CHR
        assert parts[1].isdigit()  # BP
        assert parts[2] in ["A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT"]
        assert parts[3] in ["A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT"]

    def test_remove_duplicates_keeps_smallest_p(self, minimal_gwas_df):
        """remove_duplicates=True keeps row with smallest P-value."""
        df = minimal_gwas_df.copy()
        df.loc[1] = df.loc[0].copy()
        df.loc[0, "P"] = 0.01
        df.loc[1, "P"] = 0.001
        result = make_SNPID_unique(df, remove_duplicates=True)
        duplicated_snpid = f"{df.loc[0, 'CHR']}-{df.loc[0, 'BP']}-{min(df.loc[0, 'EA'], df.loc[0, 'NEA'])}-{max(df.loc[0, 'EA'], df.loc[0, 'NEA'])}"
        filtered = result[result[ColName.SNPID] == duplicated_snpid]
        assert len(filtered) == 1
        assert filtered[ColName.P].iloc[0] == 0.001

    def test_remove_duplicates_false_adds_suffix(self, minimal_gwas_df):
        """remove_duplicates=False adds suffix -1, -2 etc."""
        df = minimal_gwas_df.copy()
        df.loc[1] = df.loc[0].copy()
        df.loc[2] = df.loc[0].copy()
        result = make_SNPID_unique(df, remove_duplicates=False)
        base_snpid = f"{df.loc[0, 'CHR']}-{df.loc[0, 'BP']}-{min(df.loc[0, 'EA'], df.loc[0, 'NEA'])}-{max(df.loc[0, 'EA'], df.loc[0, 'NEA'])}"
        snpids = result[ColName.SNPID].tolist()
        # First occurrence has no suffix, second has -1, third has -2
        assert base_snpid in snpids
        assert f"{base_snpid}-1" in snpids
        assert f"{base_snpid}-2" in snpids

    def test_allele_sorting(self):
        """G-A and A-G produce same SNPID (alleles sorted alphabetically)."""
        df1 = pd.DataFrame({
            "CHR": [1], "BP": [1000], "EA": ["G"], "NEA": ["A"],
            "BETA": [0.1], "SE": [0.01], "P": [0.05]
        })
        df2 = pd.DataFrame({
            "CHR": [1], "BP": [1000], "EA": ["A"], "NEA": ["G"],
            "BETA": [0.1], "SE": [0.01], "P": [0.05]
        })
        result1 = make_SNPID_unique(df1, remove_duplicates=False)
        result2 = make_SNPID_unique(df2, remove_duplicates=False)
        assert result1[ColName.SNPID].iloc[0] == result2[ColName.SNPID].iloc[0]
        assert "1-1000-A-G" in result1[ColName.SNPID].iloc[0]

    def test_snpid_is_first_column(self, minimal_gwas_df):
        """SNPID is the first column."""
        result = make_SNPID_unique(minimal_gwas_df, remove_duplicates=False)
        assert result.columns[0] == ColName.SNPID

    def test_custom_column_names(self):
        """Custom column names work correctly."""
        df = pd.DataFrame({
            "chromosome": [1], "position": [1000],
            "effect_allele": ["A"], "other_allele": ["G"],
            "pvalue": [0.05]
        })
        result = make_SNPID_unique(
            df, remove_duplicates=False,
            col_chr="chromosome", col_bp="position",
            col_ea="effect_allele", col_nea="other_allele",
            col_p="pvalue"
        )
        assert ColName.SNPID in result.columns
        assert result[ColName.SNPID].iloc[0] == "1-1000-A-G"

    def test_no_duplicates_no_change(self, minimal_gwas_df):
        """No duplicates results in no change to data."""
        original_len = len(minimal_gwas_df)
        result = make_SNPID_unique(minimal_gwas_df, remove_duplicates=True)
        assert len(result) == original_len


class TestRemoveAllNAColumns:
    """Tests for _remove_all_na_columns() function."""

    def test_removes_all_na_column(self, minimal_gwas_df):
        """All-NA column is removed."""
        df = minimal_gwas_df.copy()
        df["all_na_col"] = np.nan
        result = _remove_all_na_columns(df)
        assert "all_na_col" not in result.columns

    def test_empty_string_column_removed(self, minimal_gwas_df):
        """Empty string column is converted to NaN and removed."""
        df = minimal_gwas_df.copy()
        df["empty_str_col"] = ""
        result = _remove_all_na_columns(df)
        assert "empty_str_col" not in result.columns

    def test_no_all_na_columns_unchanged(self, minimal_gwas_df):
        """No all-NA columns results in unchanged dataframe."""
        result = _remove_all_na_columns(minimal_gwas_df)
        assert set(result.columns) == set(minimal_gwas_df.columns)


class TestTransformChr:
    """Tests for transform_chr() function."""

    def test_chr_prefix_removal(self):
        """'chr1' is converted to 1."""
        series = pd.Series(["chr1", "CHR2", "Chr3"])
        result = transform_chr(series)
        assert result.tolist() == [1, 2, 3]

    def test_x_chromosome_conversion(self):
        """'X' and 'x' are converted to 23."""
        series = pd.Series(["X", "x"])
        result = transform_chr(series)
        assert result.tolist() == [23, 23]

    def test_invalid_chr_to_nan(self):
        """'abc' is converted to NaN."""
        series = pd.Series(["abc", "xyz"])
        result = transform_chr(series)
        assert result.isna().all()

    def test_numeric_string_conversion(self):
        """Numeric string '5' is converted to 5."""
        series = pd.Series(["5", "10", "22"])
        result = transform_chr(series)
        assert result.tolist() == [5, 10, 22]


class TestTransformAllele:
    """Tests for transform_allele() function."""

    def test_lowercase_to_uppercase(self):
        """Lowercase alleles converted to uppercase: 'acgt' -> 'ACGT'."""
        series = pd.Series(["a", "c", "g", "t", "acgt"])
        result = transform_allele(series)
        assert result.tolist() == ["A", "C", "G", "T", "ACGT"]

    def test_invalid_characters_to_nan(self):
        """Invalid characters converted to NaN: '123', 'N'."""
        series = pd.Series(["123", "N", "X", "abc123"])
        result = transform_allele(series)
        assert result.isna().all()

    def test_multibase_valid_preserved(self):
        """Multi-base valid alleles like 'ACG' are preserved."""
        series = pd.Series(["ACG", "ATCG", "GC"])
        result = transform_allele(series)
        assert result.tolist() == ["ACG", "ATCG", "GC"]


class TestMungeChr:
    """Tests for _munge_chr() function."""

    def test_chr_prefix_removal_and_x_conversion(self):
        """CHR prefix removal and X->23 conversion."""
        df = pd.DataFrame({
            "CHR": ["chr1", "X", "chr5"],
            "BP": [1000, 2000, 3000],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "A"],
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.01, 0.02, 0.03],
            "P": [0.05, 0.01, 0.001]
        })
        result = _munge_chr(df)
        assert result[ColName.CHR].tolist() == [1, 23, 5]


class TestMungeBp:
    """Tests for _munge_bp() function."""

    def test_bp_range_validation(self):
        """BP range validation removes -1 and values > 300000000."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "BP": [-1, 1000, 300000001, 100000],
            "EA": ["A", "C", "G", "T"],
            "NEA": ["G", "T", "A", "C"],
            "BETA": [0.1, 0.2, 0.3, 0.4],
            "SE": [0.01, 0.02, 0.03, 0.04],
            "P": [0.05, 0.01, 0.001, 0.0001]
        })
        result = _munge_bp(df)
        assert len(result) == 2
        assert -1 not in result[ColName.BP].values
        assert 300000001 not in result[ColName.BP].values


class TestMungeAlleles:
    """Tests for _munge_alleles() function."""

    def test_invalid_allele_and_identical_removed(self):
        """Invalid allele and EA==NEA rows are removed."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "BP": [1000, 2000, 3000, 4000],
            "EA": ["A", "123", "G", "A"],
            "NEA": ["G", "C", "G", "C"],
            "BETA": [0.1, 0.2, 0.3, 0.4],
            "SE": [0.01, 0.02, 0.03, 0.04],
            "P": [0.05, 0.01, 0.001, 0.0001]
        })
        result = _munge_alleles(df)
        assert len(result) == 2
        assert "123" not in result[ColName.EA].values
        # Check no identical alleles
        assert (result[ColName.EA] != result[ColName.NEA]).all()


class TestMungePvalue:
    """Tests for _munge_pvalue() function."""

    def test_pvalue_zero_removed(self):
        """P=0 is removed (exclude_min)."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "BP": [1000, 2000, 3000],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "A"],
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0.01, 0.02, 0.03],
            "P": [0, 0.05, 0.01]
        })
        result = _munge_pvalue(df)
        assert len(result) == 2
        assert 0 not in result[ColName.P].values


class TestMungeSe:
    """Tests for _munge_se() function."""

    def test_se_zero_removed(self):
        """SE=0 is removed (exclude_min)."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1],
            "BP": [1000, 2000, 3000],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "A"],
            "BETA": [0.1, 0.2, 0.3],
            "SE": [0, 0.02, 0.03],
            "P": [0.05, 0.01, 0.001]
        })
        result = _munge_se(df)
        assert len(result) == 2
        assert 0 not in result[ColName.SE].values


class TestMungeEaf:
    """Tests for _munge_eaf() function."""

    def test_eaf_range_validation(self):
        """EAF range [0,1] is validated."""
        df = pd.DataFrame({
            "CHR": [1, 1, 1, 1],
            "BP": [1000, 2000, 3000, 4000],
            "EA": ["A", "C", "G", "T"],
            "NEA": ["G", "T", "A", "C"],
            "BETA": [0.1, 0.2, 0.3, 0.4],
            "SE": [0.01, 0.02, 0.03, 0.04],
            "P": [0.05, 0.01, 0.001, 0.0001],
            "EAF": [-0.1, 0.5, 1.5, 0.3]
        })
        result = _munge_eaf(df)
        assert len(result) == 2
        assert (result[ColName.EAF] >= 0).all()
        assert (result[ColName.EAF] <= 1).all()


class TestFinalizeColumns:
    """Tests for _finalize_columns() function."""

    def test_missing_columns_filled_with_none(self, minimal_gwas_df):
        """Missing columns are filled with None."""
        df = minimal_gwas_df.copy()
        result = _finalize_columns(df)
        assert ColName.N in result.columns
        assert ColName.RSID in result.columns
        assert result[ColName.N].isna().all()
        assert result[ColName.RSID].isna().all()

    def test_correct_column_order(self, minimal_gwas_df):
        """Output has correct column order."""
        result = _finalize_columns(minimal_gwas_df)
        assert list(result.columns) == OUTPUT_COLS
