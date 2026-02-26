"""Unit tests for GWAS summary statistics validation functions."""

import numpy as np
import pandas as pd
import pytest
from credtools.preprocessing.munging.validation import (
    check_mandatory_cols,
    validate_allele_consistency,
    validate_and_clean_column,
    validate_frequency_consistency,
    validate_statistical_consistency,
)
from credtools.preprocessing.munging.constants import ColName


class TestCheckMandatoryCols:
    """Tests for check_mandatory_cols function."""

    def test_all_columns_present(self):
        """All mandatory columns present should not raise error."""
        df = pd.DataFrame(
            {
                "CHR": [1],
                "BP": [1000],
                "EA": ["A"],
                "NEA": ["G"],
                "BETA": [0.1],
                "SE": [0.01],
                "P": [0.05],
            }
        )
        check_mandatory_cols(df)  # Should not raise

    def test_missing_one_column(self):
        """Missing one mandatory column should raise ValueError with column name."""
        df = pd.DataFrame(
            {
                "CHR": [1],
                "BP": [1000],
                "EA": ["A"],
                "NEA": ["G"],
                "BETA": [0.1],
                "SE": [0.01],
                # Missing "P"
            }
        )
        with pytest.raises(ValueError, match="P"):
            check_mandatory_cols(df)

    def test_missing_multiple_columns(self):
        """Missing multiple columns should raise ValueError."""
        df = pd.DataFrame(
            {
                "CHR": [1],
                "BP": [1000],
                "EA": ["A"],
                # Missing NEA, BETA, SE, P
            }
        )
        with pytest.raises(ValueError, match="Missing mandatory columns"):
            check_mandatory_cols(df)

    def test_extra_columns_present(self):
        """Extra columns beyond mandatory ones should not raise error."""
        df = pd.DataFrame(
            {
                "CHR": [1],
                "BP": [1000],
                "EA": ["A"],
                "NEA": ["G"],
                "BETA": [0.1],
                "SE": [0.01],
                "P": [0.05],
                "EXTRA1": ["value"],
                "EXTRA2": [123],
            }
        )
        check_mandatory_cols(df)  # Should not raise


class TestValidateAndCleanColumn:
    """Tests for validate_and_clean_column function."""

    def test_column_not_in_df_allow_na_true(self):
        """Column not in df with allow_na=True should return original df unchanged."""
        df = pd.DataFrame({"A": [1, 2, 3]})
        result = validate_and_clean_column(df, "B", np.float64, allow_na=True)
        pd.testing.assert_frame_equal(result, df)

    def test_column_not_in_df_allow_na_false(self):
        """Column not in df with allow_na=False should raise ValueError."""
        df = pd.DataFrame({"A": [1, 2, 3]})
        with pytest.raises(ValueError, match="Required column 'B' not found"):
            validate_and_clean_column(df, "B", np.float64, allow_na=False)

    def test_transform_func_applied(self):
        """Transform function should be applied correctly."""
        df = pd.DataFrame({"A": [1.0, 2.0, 3.0]})
        result = validate_and_clean_column(
            df, "A", np.float64, transform_func=lambda x: x * 2
        )
        expected = pd.DataFrame({"A": [2.0, 4.0, 6.0]})
        pd.testing.assert_frame_equal(result, expected)

    def test_allow_na_false_removes_na_rows(self):
        """allow_na=False should remove NA rows."""
        df = pd.DataFrame({"A": [1.0, np.nan, 3.0]})
        result = validate_and_clean_column(df, "A", np.float64, allow_na=False)
        expected = pd.DataFrame({"A": [1.0, 3.0]}, index=[0, 2])
        pd.testing.assert_frame_equal(result, expected)

    def test_allow_na_true_preserves_na_rows(self):
        """allow_na=True should preserve NA rows."""
        df = pd.DataFrame({"A": [1.0, np.nan, 3.0]})
        result = validate_and_clean_column(df, "A", np.float64, allow_na=True)
        pd.testing.assert_frame_equal(result, df)

    def test_non_numeric_strings_coerced_to_nan(self):
        """Non-numeric strings should be coerced to NaN."""
        df = pd.DataFrame({"A": ["1.0", "invalid", "3.0"]})
        result = validate_and_clean_column(df, "A", np.float64, allow_na=True)
        assert pd.isna(result.loc[1, "A"])
        assert result.loc[0, "A"] == 1.0
        assert result.loc[2, "A"] == 3.0

    def test_min_val_inclusive_boundary_preserved(self):
        """min_val inclusive (exclude_min=False) should preserve boundary value."""
        df = pd.DataFrame({"A": [0.9, 1.0, 1.1]})
        result = validate_and_clean_column(
            df, "A", np.float64, min_val=1.0, exclude_min=False
        )
        expected = pd.DataFrame({"A": [1.0, 1.1]}, index=[1, 2])
        pd.testing.assert_frame_equal(result, expected)

    def test_min_val_exclusive_boundary_removed(self):
        """min_val exclusive (exclude_min=True) should remove boundary value."""
        df = pd.DataFrame({"A": [0.9, 1.0, 1.1]})
        result = validate_and_clean_column(
            df, "A", np.float64, min_val=1.0, exclude_min=True
        )
        expected = pd.DataFrame({"A": [1.1]}, index=[2])
        pd.testing.assert_frame_equal(result, expected)

    def test_max_val_inclusive_boundary_preserved(self):
        """max_val inclusive (exclude_max=False) should preserve boundary value."""
        df = pd.DataFrame({"A": [0.9, 1.0, 1.1]})
        result = validate_and_clean_column(
            df, "A", np.float64, max_val=1.0, exclude_max=False
        )
        expected = pd.DataFrame({"A": [0.9, 1.0]}, index=[0, 1])
        pd.testing.assert_frame_equal(result, expected)

    def test_max_val_exclusive_boundary_removed(self):
        """max_val exclusive (exclude_max=True) should remove boundary value."""
        df = pd.DataFrame({"A": [0.9, 1.0, 1.1]})
        result = validate_and_clean_column(
            df, "A", np.float64, max_val=1.0, exclude_max=True
        )
        expected = pd.DataFrame({"A": [0.9]}, index=[0])
        pd.testing.assert_frame_equal(result, expected)

    def test_range_filtering_preserves_nan_with_allow_na(self):
        """Range filtering should preserve NaN rows when allow_na=True."""
        df = pd.DataFrame({"A": [0.5, np.nan, 1.5, 2.5]})
        result = validate_and_clean_column(
            df, "A", np.float64, min_val=1.0, max_val=2.0, allow_na=True
        )
        assert len(result) == 2
        assert pd.isna(result.loc[1, "A"])
        assert result.loc[2, "A"] == 1.5

    def test_col_type_str_converts_nan_string(self):
        """col_type=str should convert 'nan' string to np.nan."""
        df = pd.DataFrame({"A": ["value1", "nan", "value2"]})
        result = validate_and_clean_column(df, "A", str)
        assert result.loc[0, "A"] == "value1"
        assert pd.isna(result.loc[1, "A"])
        assert result.loc[2, "A"] == "value2"

    @pytest.mark.parametrize(
        "col_type,expected_dtype",
        [
            (np.int8, np.int8),
            (np.float32, np.float32),
            (np.float64, np.float64),
        ],
    )
    def test_col_type_numeric_conversion(self, col_type, expected_dtype):
        """Different numeric types should be converted correctly."""
        df = pd.DataFrame({"A": [1.0, 2.0, 3.0]})
        result = validate_and_clean_column(df, "A", col_type)
        assert result["A"].dtype == expected_dtype

    def test_all_valid_data_row_count_unchanged(self):
        """All valid data should result in unchanged row count."""
        df = pd.DataFrame({"A": [1.0, 2.0, 3.0, 4.0, 5.0]})
        result = validate_and_clean_column(
            df, "A", np.float64, min_val=0.0, max_val=10.0
        )
        assert len(result) == len(df)

    def test_all_invalid_data_returns_empty_df(self):
        """All invalid data should return empty dataframe."""
        df = pd.DataFrame({"A": [1.0, 2.0, 3.0]})
        result = validate_and_clean_column(
            df, "A", np.float64, min_val=10.0, allow_na=False
        )
        assert len(result) == 0

    def test_input_df_not_modified(self):
        """Input dataframe should not be modified (copy semantics)."""
        df = pd.DataFrame({"A": [1.0, 2.0, 3.0]})
        df_original = df.copy()
        validate_and_clean_column(df, "A", np.float64, min_val=2.0)
        pd.testing.assert_frame_equal(df, df_original)


class TestValidateAlleleConsistency:
    """Tests for validate_allele_consistency function."""

    def test_ea_equals_nea_rows_removed(self):
        """Rows where EA equals NEA should be removed."""
        df = pd.DataFrame(
            {
                "EA": ["A", "A", "C"],
                "NEA": ["G", "A", "T"],
            }
        )
        result = validate_allele_consistency(df)
        assert len(result) == 2
        assert 1 not in result.index

    def test_valid_rows_preserved(self):
        """Valid rows with different alleles should be preserved."""
        df = pd.DataFrame(
            {
                "EA": ["A", "C", "G"],
                "NEA": ["G", "T", "A"],
            }
        )
        result = validate_allele_consistency(df)
        assert len(result) == 3
        pd.testing.assert_frame_equal(result, df)

    def test_na_allele_rows_removed(self):
        """Rows with NA alleles should be removed."""
        df = pd.DataFrame(
            {
                "EA": ["A", np.nan, "C", "G"],
                "NEA": ["G", "T", np.nan, "A"],
            }
        )
        result = validate_allele_consistency(df)
        assert len(result) == 2
        assert result.loc[0, "EA"] == "A"
        assert result.loc[3, "EA"] == "G"

    def test_no_ea_nea_columns_returns_original(self):
        """No EA/NEA columns should return original dataframe."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        result = validate_allele_consistency(df)
        pd.testing.assert_frame_equal(result, df)


class TestValidateStatisticalConsistency:
    """Tests for validate_statistical_consistency function."""

    def test_se_zero_or_negative_rows_removed(self):
        """Rows with SE=0 or SE<0 should be removed."""
        df = pd.DataFrame(
            {
                "BETA": [0.1, 0.2, 0.3, 0.4],
                "SE": [0.01, 0.0, -0.01, 0.02],
            }
        )
        result = validate_statistical_consistency(df)
        assert len(result) == 2
        assert 1 not in result.index
        assert 2 not in result.index

    def test_se_positive_preserved(self):
        """Rows with SE>0 should be preserved."""
        df = pd.DataFrame(
            {
                "BETA": [0.1, 0.2, 0.3],
                "SE": [0.01, 0.02, 0.03],
            }
        )
        result = validate_statistical_consistency(df)
        pd.testing.assert_frame_equal(result, df)

    def test_no_beta_se_columns_returns_original(self):
        """No BETA/SE columns should return original dataframe."""
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        result = validate_statistical_consistency(df)
        pd.testing.assert_frame_equal(result, df)


class TestValidateFrequencyConsistency:
    """Tests for validate_frequency_consistency function."""

    def test_maf_greater_than_half_flipped(self):
        """MAF>0.5 should be flipped to 1-MAF."""
        df = pd.DataFrame({"MAF": [0.2, 0.6, 0.8, 0.3]})
        result = validate_frequency_consistency(df)
        assert result.loc[0, "MAF"] == pytest.approx(0.2)
        assert result.loc[1, "MAF"] == pytest.approx(0.4)
        assert result.loc[2, "MAF"] == pytest.approx(0.2)
        assert result.loc[3, "MAF"] == pytest.approx(0.3)

    def test_eaf_maf_discrepancy_corrected(self):
        """EAF/MAF discrepancy should be corrected."""
        df = pd.DataFrame(
            {
                "EAF": [0.2, 0.8, 0.3],
                "MAF": [0.5, 0.2, 0.3],  # Wrong values
            }
        )
        result = validate_frequency_consistency(df)
        # Expected MAF from EAF: min(0.2, 0.8)=0.2, min(0.8, 0.2)=0.2, min(0.3, 0.7)=0.3
        assert abs(result.loc[0, "MAF"] - 0.2) < 0.01
        assert abs(result.loc[1, "MAF"] - 0.2) < 0.01
        assert abs(result.loc[2, "MAF"] - 0.3) < 0.01

    def test_no_maf_column_returns_original(self):
        """No MAF column should return original dataframe."""
        df = pd.DataFrame({"EAF": [0.2, 0.3, 0.4], "OTHER": [1, 2, 3]})
        result = validate_frequency_consistency(df)
        pd.testing.assert_frame_equal(result, df)


# ============================================================================
# Tests for validate_pvalue_consistency
# ============================================================================

import logging
from scipy.stats import norm
from credtools.preprocessing.munging.validation import validate_pvalue_consistency


class TestValidatePvalueConsistency:
    """Tests for validate_pvalue_consistency function."""

    def test_consistent_pvalues_no_warning(self, caplog):
        """Consistent BETA, SE, P should not trigger a warning."""
        betas = [0.5, -1.0, 0.2]
        ses = [0.1, 0.25, 0.05]
        # Compute truly consistent p-values from BETA/SE
        pvals = [
            float(2 * (1 - norm.cdf(abs(b / s)))) for b, s in zip(betas, ses)
        ]
        df = pd.DataFrame({"BETA": betas, "SE": ses, "P": pvals})

        with caplog.at_level(logging.WARNING, logger="Munging"):
            result = validate_pvalue_consistency(df)

        # No warning should be logged
        assert "P-value discrepancies" not in caplog.text
        # All rows should be preserved
        assert len(result) == len(df)

    def test_large_discrepancy_triggers_warning(self, caplog):
        """P-value far from expected should trigger a discrepancy warning."""
        # BETA=1.0, SE=0.1 => z=10 => expected_p ~ 1.5e-23
        # We set P=0.5, which is wildly different
        df = pd.DataFrame(
            {"BETA": [1.0], "SE": [0.1], "P": [0.5]}
        )

        with caplog.at_level(logging.WARNING, logger="Munging"):
            result = validate_pvalue_consistency(df)

        assert "P-value discrepancies" in caplog.text
        # Rows are preserved (function only warns, does not drop)
        assert len(result) == 1

    def test_missing_beta_returns_unchanged(self):
        """Return unchanged when BETA column is missing."""
        df = pd.DataFrame({"SE": [0.1, 0.2], "P": [0.05, 0.01]})
        result = validate_pvalue_consistency(df)
        pd.testing.assert_frame_equal(result, df)

    def test_missing_se_returns_unchanged(self):
        """Return unchanged when SE column is missing."""
        df = pd.DataFrame({"BETA": [0.1, 0.2], "P": [0.05, 0.01]})
        result = validate_pvalue_consistency(df)
        pd.testing.assert_frame_equal(result, df)
