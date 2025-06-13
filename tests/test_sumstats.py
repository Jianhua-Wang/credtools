import gzip
import io
from unittest.mock import mock_open, patch

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.sumstats import (
    check_mandatory_cols,
    get_significant_snps,
    load_sumstats,
    make_SNPID_unique,
    munge,
    rm_col_allna,
)


@pytest.fixture
def sample_df():
    return pd.DataFrame(
        {
            ColName.SNPID: ["rs1", "rs2", "rs3", "rs4", "rs5"],
            ColName.P: [1e-9, 0.05, 1e-7, 0.1, 1e-8],
        }
    )


def test_get_significant_snps_default():
    df = pd.DataFrame(
        {ColName.SNPID: ["rs1", "rs2", "rs3"], ColName.P: [1e-9, 0.05, 1e-8]}
    )
    result = get_significant_snps(df)
    assert len(result) == 2
    assert result.iloc[0][ColName.SNPID] == "rs1"
    assert result.iloc[1][ColName.SNPID] == "rs3"
    assert result.iloc[0][ColName.P] == 1e-9
    assert result.iloc[1][ColName.P] == 1e-8


def test_get_significant_snps_custom_threshold():
    df = pd.DataFrame(
        {ColName.SNPID: ["rs1", "rs2", "rs3"], ColName.P: [1e-9, 0.05, 1e-7]}
    )
    result = get_significant_snps(df, pvalue_threshold=1e-8)
    assert len(result) == 1
    assert result.iloc[0][ColName.SNPID] == "rs1"
    assert result.iloc[0][ColName.P] == 1e-9


def test_get_significant_snps_no_sig_use_most_sig():
    df = pd.DataFrame(
        {ColName.SNPID: ["rs1", "rs2", "rs3"], ColName.P: [0.1, 0.05, 0.2]}
    )
    result = get_significant_snps(df, pvalue_threshold=1e-8)
    assert len(result) == 1
    assert result.iloc[0][ColName.SNPID] == "rs2"
    assert result.iloc[0][ColName.P] == 0.05


def test_get_significant_snps_no_sig_dont_use_most_sig():
    df = pd.DataFrame(
        {ColName.SNPID: ["rs1", "rs2", "rs3"], ColName.P: [0.1, 0.05, 0.2]}
    )
    with pytest.raises(ValueError, match="No significant SNPs found."):
        get_significant_snps(df, pvalue_threshold=1e-8, use_most_sig_if_no_sig=False)


def test_get_significant_snps_missing_columns():
    df = pd.DataFrame({"SNPID": ["rs1", "rs2", "rs3"], "some_column": [1, 2, 3]})
    with pytest.raises(KeyError, match="The following required columns are missing"):
        get_significant_snps(df)


def test_get_significant_snps_empty_dataframe():
    df = pd.DataFrame(columns=[ColName.SNPID, ColName.P])
    with pytest.raises(
        ValueError, match="The DataFrame is empty. No SNPs available to select."
    ):
        get_significant_snps(df)


def test_get_significant_snps_sorting(sample_df: pd.DataFrame):
    result = get_significant_snps(sample_df, pvalue_threshold=1e-7)
    assert list(result[ColName.SNPID]) == ["rs1", "rs5", "rs3"]
    assert list(result[ColName.P]) == [1e-9, 1e-8, 1e-7]


def test_get_significant_snps_float_precision():
    df = pd.DataFrame(
        {
            ColName.SNPID: ["rs1", "rs2", "rs3"],
            ColName.P: [4.9999999e-8, 5.0000001e-8, 5e-8],
        }
    )
    result = get_significant_snps(df)
    assert len(result) == 2
    assert result.iloc[0][ColName.SNPID] == "rs1"


def test_get_significant_snps_large_dataframe():
    np.random.seed(42)
    large_df = pd.DataFrame(
        {
            ColName.SNPID: [f"rs{i}" for i in range(1000000)],
            ColName.P: np.random.random(1000000),
        }
    )
    large_df.loc[0, ColName.P] = 1e-9  # Ensure at least one significant SNP
    result = get_significant_snps(large_df)
    assert len(result) > 0
    assert result.iloc[0][ColName.P] == 1e-9


@pytest.fixture
def sample_df1():
    return pd.DataFrame(
        {
            ColName.CHR: ["1", "1", "2", "2"],
            ColName.BP: [12345, 12345, 67890, 67890],
            ColName.EA: ["A", "A", "G", "A"],
            ColName.NEA: ["G", "G", "A", "G"],
            "rsID": ["rs1", "rs2", "rs3", "rs4"],
            ColName.P: [1e-5, 1e-6, 1e-7, 1e-8],
        }
    )


def test_make_SNPID_unique_basic(sample_df1: pd.DataFrame):
    result = make_SNPID_unique(sample_df1)
    assert len(result) == 2
    assert list(result[ColName.SNPID]) == ["1-12345-A-G", "2-67890-A-G"]
    assert result.iloc[0][ColName.P] == 1e-6
    assert result.iloc[1][ColName.P] == 1e-8


def test_make_SNPID_unique_no_remove_duplicates(sample_df1: pd.DataFrame):
    result = make_SNPID_unique(sample_df1, remove_duplicates=False)
    assert len(result) == 4
    assert list(result[ColName.SNPID]) == [
        "1-12345-A-G",
        "1-12345-A-G-1",
        "2-67890-A-G",
        "2-67890-A-G-1",
    ]


def test_make_SNPID_unique_custom_column_names():
    df = pd.DataFrame(
        {
            "CHROM": ["1", "2"],
            "POS": [12345, 67890],
            "ALT": ["A", "G"],
            "REF": ["G", "A"],
            "PVAL": [1e-5, 1e-7],
        }
    )
    result = make_SNPID_unique(
        df, col_chr="CHROM", col_bp="POS", col_ea="ALT", col_nea="REF", col_p="PVAL"
    )
    assert len(result) == 2
    assert list(result[ColName.SNPID]) == ["1-12345-A-G", "2-67890-A-G"]


def test_make_SNPID_unique_missing_columns():
    df = pd.DataFrame(
        {ColName.CHR: ["1", "2"], ColName.BP: [12345, 67890], ColName.EA: ["A", "G"]}
    )
    with pytest.raises(KeyError, match="The following required columns are missing"):
        make_SNPID_unique(df)


def test_make_SNPID_unique_empty_dataframe():
    df = pd.DataFrame(columns=[ColName.CHR, ColName.BP, ColName.EA, ColName.NEA])
    with pytest.raises(ValueError, match="The input DataFrame is empty."):
        make_SNPID_unique(df)


@patch("logging.Logger.debug")
@patch("logging.Logger.warning")
def test_make_SNPID_unique_logging(mock_warning, mock_debug, sample_df1: pd.DataFrame):
    result = make_SNPID_unique(sample_df1)
    mock_debug.assert_any_call("Number of duplicated SNPs: 2")
    mock_debug.assert_any_call("Unique SNPIDs have been successfully created.")
    mock_debug.assert_any_call("Total unique SNPs: 2")

    make_SNPID_unique(sample_df1, remove_duplicates=False)
    mock_warning.assert_called_with(
        "Duplicated SNPs detected. To remove duplicates, set `remove_duplicates=True`.\n"
        "            Change the Unique SNP identifier to make it unique."
    )


def test_make_SNPID_unique_allele_sorting():
    df = pd.DataFrame(
        {
            ColName.CHR: ["1", "1"],
            ColName.BP: [12345, 12345],
            ColName.EA: ["G", "A"],
            ColName.NEA: ["A", "G"],
            ColName.P: [1e-5, 1e-6],
        }
    )
    result = make_SNPID_unique(df)
    assert len(result) == 1
    assert result.iloc[0][ColName.SNPID] == "1-12345-A-G"
    assert result.iloc[0][ColName.P] == 1e-6


def test_make_SNPID_unique_large_dataframe():
    np.random.seed(42)
    large_df = pd.DataFrame(
        {
            ColName.CHR: np.random.choice(["1", "2", "3", "4", "5"], 100000),
            ColName.BP: np.random.randint(1, 1000000, 100000),
            ColName.EA: np.random.choice(["A", "C", "G", "T"], 100000),
            ColName.NEA: np.random.choice(["A", "C", "G", "T"], 100000),
            ColName.P: np.random.random(100000),
        }
    )
    result = make_SNPID_unique(large_df)
    assert len(result) <= len(large_df)
    assert ColName.SNPID in result.columns
    assert result[ColName.SNPID].is_unique


def test_make_SNPID_unique_column_order():
    df = pd.DataFrame(
        {
            ColName.CHR: ["1", "2"],
            ColName.BP: [12345, 67890],
            ColName.EA: ["A", "G"],
            ColName.NEA: ["G", "A"],
            ColName.P: [1e-5, 1e-7],
            "OTHER": [1, 2],
        }
    )
    result = make_SNPID_unique(df)
    assert list(result.columns) == [
        ColName.SNPID,
        ColName.CHR,
        ColName.BP,
        ColName.EA,
        ColName.NEA,
        ColName.P,
        "OTHER",
    ]


def test_make_SNPID_unique_all_duplicates_removed():
    df = pd.DataFrame(
        {
            ColName.CHR: ["1", "1", "1"],
            ColName.BP: [12345, 12345, 12345],
            ColName.EA: ["A", "A", "A"],
            ColName.NEA: ["G", "G", "G"],
            ColName.P: [1e-5, 1e-6, 1e-7],
        }
    )
    result = make_SNPID_unique(df)
    assert len(result) == 1
    assert result.iloc[0][ColName.P] == 1e-7


@pytest.fixture
def sample_df3():
    return pd.DataFrame(
        {
            ColName.CHR: ["1", "2", "3", "X"],
            ColName.BP: [1000, 2000, 3000, 4000],
            ColName.EA: ["A", "C", "G", "T"],
            ColName.NEA: ["G", "T", "C", "A"],
            ColName.P: [0.01, 0.001, 0.0001, 0.00001],
            ColName.BETA: [0.1, 0.2, 0.3, 0.4],
            ColName.SE: [0.05, 0.06, 0.07, 0.08],
            ColName.EAF: [0.3, 0.4, 0.5, 0.6],
            ColName.RSID: ["rs1", "rs2", "rs3", "rs4"],
        }
    )


@patch("credtools.sumstats.check_mandatory_cols")
@patch("credtools.sumstats.rm_col_allna")
@patch("credtools.sumstats.munge_chr")
@patch("credtools.sumstats.munge_bp")
@patch("credtools.sumstats.munge_allele")
@patch("credtools.sumstats.make_SNPID_unique")
@patch("credtools.sumstats.munge_pvalue")
@patch("credtools.sumstats.munge_beta")
@patch("credtools.sumstats.munge_se")
@patch("credtools.sumstats.munge_eaf")
@patch("credtools.sumstats.munge_maf")
@patch("credtools.sumstats.munge_rsid")
@patch("credtools.sumstats.check_colnames")
def test_munge_function_calls(
    mock_check_colnames,
    mock_munge_rsid,
    mock_munge_maf,
    mock_munge_eaf,
    mock_munge_se,
    mock_munge_beta,
    mock_munge_pvalue,
    mock_make_SNPID_unique,
    mock_munge_allele,
    mock_munge_bp,
    mock_munge_chr,
    mock_rm_col_allna,
    mock_check_mandatory_cols,
    sample_df3: pd.DataFrame,
):
    # Set up all mocks to return the input dataframe
    for mock in [
        mock_check_colnames,
        mock_munge_rsid,
        mock_munge_maf,
        mock_munge_eaf,
        mock_munge_se,
        mock_munge_beta,
        mock_munge_pvalue,
        mock_make_SNPID_unique,
        mock_munge_allele,
        mock_munge_bp,
        mock_munge_chr,
        mock_rm_col_allna,
    ]:
        mock.side_effect = lambda df: df

    result = munge(sample_df3)

    # Assert that all functions were called
    mock_check_mandatory_cols.assert_called_once()
    mock_rm_col_allna.assert_called_once()
    mock_munge_chr.assert_called_once()
    mock_munge_bp.assert_called_once()
    mock_munge_allele.assert_called_once()
    mock_make_SNPID_unique.assert_called_once()
    mock_munge_pvalue.assert_called_once()
    mock_munge_beta.assert_called_once()
    mock_munge_se.assert_called_once()
    mock_munge_eaf.assert_called_once()
    mock_munge_maf.assert_called_once()
    mock_munge_rsid.assert_called_once()
    mock_check_colnames.assert_called_once()

    # Check that the result is a DataFrame
    assert isinstance(result, pd.DataFrame)


def test_munge_sorting(sample_df3: pd.DataFrame):
    result = munge(sample_df3)
    assert list(result[ColName.CHR]) == [1, 2, 3, 23]
    assert list(result[ColName.BP]) == [1000, 2000, 3000, 4000]


def test_munge_maf_column(sample_df3: pd.DataFrame):
    result = munge(sample_df3)
    assert ColName.MAF in result.columns


def test_munge_without_rsid(sample_df3: pd.DataFrame):
    sample_df3 = sample_df3.drop(columns=[ColName.RSID])
    result = munge(sample_df3)
    assert ColName.RSID in result.columns


@patch("credtools.sumstats.check_mandatory_cols")
def test_munge_missing_mandatory_cols(
    mock_check_mandatory_cols, sample_df3: pd.DataFrame
):
    mock_check_mandatory_cols.side_effect = ValueError("Missing mandatory columns")
    with pytest.raises(ValueError, match="Missing mandatory columns"):
        munge(sample_df3)


def test_munge_input_unmodified(sample_df3: pd.DataFrame):
    original_df = sample_df3.copy()
    munge(sample_df3)
    pd.testing.assert_frame_equal(sample_df3, original_df)


@patch("credtools.sumstats.rm_col_allna")
def test_munge_remove_allna_columns(mock_rm_col_allna, sample_df3: pd.DataFrame):
    mock_rm_col_allna.return_value = sample_df3.drop(columns=[ColName.RSID])
    result = munge(sample_df3)
    assert ColName.RSID in result.columns


def test_munge_large_dataframe():
    np.random.seed(42)
    large_df = pd.DataFrame(
        {
            ColName.CHR: np.random.choice(["1", "2", "3", "4", "5"], 100000),
            ColName.BP: np.random.randint(1, 1000000, 100000),
            ColName.EA: np.random.choice(["A", "C", "G", "T"], 100000),
            ColName.NEA: np.random.choice(["A", "C", "G", "T"], 100000),
            ColName.P: np.random.random(100000),
            ColName.BETA: np.random.normal(0, 1, 100000),
            ColName.SE: np.random.random(100000),
            ColName.EAF: np.random.random(100000),
        }
    )
    result = munge(large_df)
    assert len(result) <= 100000
    assert ColName.MAF in result.columns


@pytest.fixture
def sample_df4():
    return pd.DataFrame(
        {
            ColName.CHR: ["1", "2", "3", "X"],
            ColName.BP: [1000, 2000, 3000, 4000],
            ColName.EA: ["A", "C", "G", "T"],
            ColName.NEA: ["G", "T", "C", "A"],
            ColName.P: [0.01, 0.001, 0.0001, 0.00001],
            ColName.EAF: [0.3, 0.4, 0.5, 0.6],
            ColName.BETA: [0.1, 0.2, 0.3, 0.4],
            ColName.SE: [0.05, 0.06, 0.07, 0.08],
            "EXTRA": [1, 2, 3, 4],
        }
    )


def test_check_mandatory_cols_all_present(sample_df4: pd.DataFrame):
    check_mandatory_cols(sample_df4)  # Should not raise an exception


def test_check_mandatory_cols_missing(sample_df4: pd.DataFrame):
    df_missing = sample_df4.drop(columns=[ColName.CHR])
    with pytest.raises(
        ValueError, match=f"Missing mandatory columns: {{'{ColName.CHR}'}}"
    ):
        check_mandatory_cols(df_missing)


def test_check_mandatory_cols_multiple_missing(sample_df4: pd.DataFrame):
    df_missing = sample_df4.drop(columns=[ColName.CHR, ColName.BP])
    with pytest.raises(ValueError, match="Missing mandatory columns:"):
        check_mandatory_cols(df_missing)


def test_check_mandatory_cols_input_unmodified(sample_df4: pd.DataFrame):
    original_df = sample_df4.copy()
    check_mandatory_cols(sample_df4)
    pd.testing.assert_frame_equal(sample_df4, original_df)


def test_rm_col_allna_no_removal(sample_df4: pd.DataFrame):
    result = rm_col_allna(sample_df4)
    pd.testing.assert_frame_equal(result, sample_df4)


def test_rm_col_allna_remove_one():
    df = pd.DataFrame({"A": [1, 2, 3], "B": [np.nan, np.nan, np.nan], "C": [4, 5, 6]})
    result = rm_col_allna(df)
    expected = pd.DataFrame({"A": [1, 2, 3], "C": [4, 5, 6]})
    pd.testing.assert_frame_equal(result, expected)


def test_rm_col_allna_remove_multiple():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3],
            "B": [np.nan, np.nan, np.nan],
            "C": [4, 5, 6],
            "D": [np.nan, np.nan, np.nan],
        }
    )
    result = rm_col_allna(df)
    expected = pd.DataFrame({"A": [1, 2, 3], "C": [4, 5, 6]})
    pd.testing.assert_frame_equal(result, expected)


def test_rm_col_allna_empty_strings():
    df = pd.DataFrame({"A": [1, 2, 3], "B": ["", "", ""], "C": [4, 5, 6]})
    result = rm_col_allna(df)
    expected = pd.DataFrame({"A": [1, 2, 3], "C": [4, 5, 6]})
    pd.testing.assert_frame_equal(result, expected)


def test_rm_col_allna_mixed_na():
    df = pd.DataFrame({"A": [1, 2, 3], "B": [np.nan, "", None], "C": [4, 5, 6]})
    result = rm_col_allna(df)
    expected = pd.DataFrame({"A": [1, 2, 3], "C": [4, 5, 6]})
    pd.testing.assert_frame_equal(result, expected)


def test_rm_col_allna_input_unmodified():
    df = pd.DataFrame({"A": [1, 2, 3], "B": [np.nan, np.nan, np.nan], "C": [4, 5, 6]})
    original_df = df.copy()
    rm_col_allna(df)
    pd.testing.assert_frame_equal(df, original_df)


@patch("logging.Logger.debug")
def test_rm_col_allna_logging(mock_debug):
    df = pd.DataFrame({"A": [1, 2, 3], "B": [np.nan, np.nan, np.nan], "C": [4, 5, 6]})
    rm_col_allna(df)
    mock_debug.assert_called_once_with("Remove column B because it is all NA.")


@pytest.fixture
def sample_data():
    return "CHR\tBP\tEA\tNEA\tEAF\tBETA\tSE\tP\n1\t1000\tA\tG\t0.5\t0.1\t0.05\t0.01\n2\t2000\tC\tT\t0.3\t-0.2\t0.06\t0.001\n"


@pytest.fixture
def create_test_file(tmp_path, sample_data):
    def _create_file(filename, content=sample_data, gzipped=False):
        file_path = tmp_path / filename
        if gzipped:
            with gzip.open(file_path, "wt") as f:
                f.write(content)
        else:
            with open(file_path, "w") as f:
                f.write(content)
        return str(file_path)

    return _create_file


def test_load_sumstats_basic(create_test_file):
    file_path = create_test_file("test.txt")
    df = load_sumstats(file_path)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ColName.sumstat_cols
    assert len(df) == 2


def test_load_sumstats_gzipped(create_test_file):
    file_path = create_test_file("test.txt.gz", gzipped=True)
    df = load_sumstats(file_path)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2


def test_load_sumstats_custom_sep(create_test_file):
    content = "CHR,BP,EA,NEA,EAF,BETA,SE,P\n1,1000,A,G,0.5,0.1,0.05,0.01\n"
    file_path = create_test_file("test_comma.txt", content=content)
    df = load_sumstats(file_path, sep=",")
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1


def test_load_sumstats_nrows(create_test_file):
    file_path = create_test_file("test.txt")
    df = load_sumstats(file_path, nrows=1)
    assert len(df) == 1


def test_load_sumstats_skiprows(create_test_file):
    content = "# Comment line\nCHR\tBP\tEA\tNEA\tEAF\tBETA\tSE\tP\n1\t1000\tA\tG\t0.5\t0.1\t0.05\t0.01\n"
    file_path = create_test_file("test_skip.txt", content=content)
    df = load_sumstats(file_path, skiprows=1)
    assert len(df) == 1


def test_load_sumstats_comment(create_test_file):
    content = "CHR\tBP\tEA\tNEA\tEAF\tBETA\tSE\tP\n1\t1000\tA\tG\t0.5\t0.1\t0.05\t0.01\n# Comment line\n2\t2000\tC\tT\t0.3\t-0.2\t0.06\t0.001\n"
    file_path = create_test_file("test_comment.txt", content=content)
    df = load_sumstats(file_path, comment="#")
    assert len(df) == 2


def test_load_sumstats_infer_sep(create_test_file):
    content = "CHR,BP,EA,NEA,EAF,BETA,SE,P\n1,1000,A,G,0.5,0.1,0.05,0.01\n"
    file_path = create_test_file("test_infer.txt", content=content)
    df = load_sumstats(file_path)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1


def test_load_sumstats_missing_columns(create_test_file):
    content = "CHR\tBP\tEA\tNEA\tEAF\tBETA\tSE\n1\t1000\tA\tG\t0.5\t0.1\t0.05\n"
    file_path = create_test_file("test_missing.txt", content=content)
    with pytest.raises(ValueError):
        load_sumstats(file_path)


def test_load_sumstats_empty_file(create_test_file):
    file_path = create_test_file("empty.txt", content="")
    with pytest.raises(pd.errors.EmptyDataError):
        load_sumstats(file_path)


def test_load_sumstats_file_not_found():
    with pytest.raises(FileNotFoundError):
        load_sumstats("nonexistent_file.txt")
