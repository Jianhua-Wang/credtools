"""Tests for credtools.qc module."""

import os
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet

# ---------------------------------------------------------------------------
# Helper: build synthetic Locus / LocusSet
# ---------------------------------------------------------------------------


def _make_locus(
    popu: str = "EUR",
    cohort: str = "UKB",
    seed: int = 42,
    n_snps: int = 20,
    sample_size: int = 10000,
    add_af2: bool = False,
) -> Locus:
    """Create a test locus with realistic LD structure.

    Parameters
    ----------
    add_af2 : bool
        If True, add an AF2 column to the LD map (required by compare_maf).
    """
    rng = np.random.default_rng(seed)
    bps = np.arange(1000, 1000 + n_snps * 100, 100)
    snpids = [f"1-{bp}-A-G" for bp in bps]

    eaf = rng.uniform(0.1, 0.5, n_snps).astype(np.float32)
    sumstats = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: np.int8(1),
            ColName.BP: bps.astype(np.int32),
            ColName.EA: ["A"] * n_snps,
            ColName.NEA: ["G"] * n_snps,
            ColName.EAF: eaf,
            ColName.MAF: np.minimum(eaf, 1 - eaf),
            ColName.A1: ["A"] * n_snps,
            ColName.A2: ["G"] * n_snps,
            ColName.BETA: rng.normal(0, 0.1, n_snps).astype(np.float32),
            ColName.SE: rng.uniform(0.01, 0.05, n_snps).astype(np.float32),
            ColName.P: rng.uniform(1e-10, 0.05, n_snps).astype(np.float64),
        }
    )

    # Create a positive-definite LD matrix via random correlation
    A = rng.normal(size=(n_snps, n_snps))
    r = A @ A.T
    d = np.sqrt(np.diag(r))
    r = r / np.outer(d, d)
    r = r.astype(np.float32)

    ld_map = sumstats[
        [ColName.SNPID, ColName.CHR, ColName.BP, ColName.A1, ColName.A2]
    ].copy()
    if add_af2:
        ld_map["AF2"] = rng.uniform(0.1, 0.5, n_snps).astype(np.float32)

    ld = LDMatrix(ld_map, r)

    return Locus(
        popu=popu,
        cohort=cohort,
        sample_size=sample_size,
        sumstats=sumstats,
        locus_start=int(bps[0]),
        locus_end=int(bps[-1]),
        ld=ld,
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def single_locus():
    return _make_locus()


@pytest.fixture
def single_locus_with_af2():
    return _make_locus(add_af2=True)


@pytest.fixture
def single_locus_set():
    return LocusSet([_make_locus()])


@pytest.fixture
def multi_cohort_locus_set():
    locus1 = _make_locus("EUR", "UKB", seed=42)
    locus2 = _make_locus("AFR", "MVP", seed=123)
    locus3 = _make_locus("EAS", "BBJ", seed=456)
    return LocusSet([locus1, locus2, locus3])


@pytest.fixture
def single_locus_set_with_af2():
    return LocusSet([_make_locus(add_af2=True)])


@pytest.fixture
def precomputed_qc_metrics(single_locus_set):
    from credtools.qc import locus_qc

    return locus_qc(single_locus_set)


# ===================================================================
# TestGetEigen
# ===================================================================


class TestGetEigen:
    """Tests for get_eigen()."""

    def test_returns_dict_with_keys(self, single_locus):
        """Verify get_eigen returns a dict containing eigvals and eigvecs keys."""
        from credtools.qc import get_eigen

        result = get_eigen(single_locus.ld.r)
        assert "eigvals" in result
        assert "eigvecs" in result

    def test_eigvals_shape(self, single_locus):
        """Verify eigenvalues array has shape equal to the number of SNPs."""
        from credtools.qc import get_eigen

        n = single_locus.ld.r.shape[0]
        result = get_eigen(single_locus.ld.r)
        assert result["eigvals"].shape == (n,)

    def test_eigvecs_shape(self, single_locus):
        """Verify eigenvectors matrix has shape (n, n)."""
        from credtools.qc import get_eigen

        n = single_locus.ld.r.shape[0]
        result = get_eigen(single_locus.ld.r)
        assert result["eigvecs"].shape == (n, n)

    def test_dtype_conversion(self, single_locus):
        """Verify eigenvalues are cast to the requested dtype."""
        from credtools.qc import get_eigen

        result = get_eigen(single_locus.ld.r, dtype=np.float64)
        assert result["eigvals"].dtype == np.float64

    def test_positive_semidefinite_eigenvalues(self):
        """A well-formed correlation matrix should have non-negative eigenvalues."""
        from credtools.qc import get_eigen

        mat = np.eye(5, dtype=np.float32)
        result = get_eigen(mat)
        assert np.all(result["eigvals"] >= -1e-7)


# ===================================================================
# TestEstimateSRss
# ===================================================================


class TestEstimateSRss:
    """Tests for estimate_s_rss()."""

    def test_null_mle_returns_float(self, single_locus):
        """Verify null-mle method returns a float value."""
        from credtools.qc import estimate_s_rss

        s = estimate_s_rss(single_locus, method="null-mle")
        assert isinstance(s, (float, np.floating))

    def test_null_mle_s_in_range(self, single_locus):
        """Verify null-mle estimate is between 0 and 1."""
        from credtools.qc import estimate_s_rss

        s = estimate_s_rss(single_locus, method="null-mle")
        assert 0 <= s <= 1

    def test_partialmle_returns_numeric(self, single_locus):
        """Verify null-partialmle method returns a numeric value."""
        from credtools.qc import estimate_s_rss

        s = estimate_s_rss(single_locus, method="null-partialmle")
        assert isinstance(s, (int, float, np.integer, np.floating))

    def test_pseudomle_returns_float(self, single_locus):
        """Verify null-pseudomle method returns a float value."""
        from credtools.qc import estimate_s_rss

        s = estimate_s_rss(single_locus, method="null-pseudomle")
        assert isinstance(s, (float, np.floating))

    def test_pseudomle_s_in_range(self, single_locus):
        """Verify null-pseudomle estimate is between 0 and 1."""
        from credtools.qc import estimate_s_rss

        s = estimate_s_rss(single_locus, method="null-pseudomle")
        assert 0 <= s <= 1

    def test_invalid_method_raises(self, single_locus):
        """Verify an invalid method name raises ValueError."""
        from credtools.qc import estimate_s_rss

        with pytest.raises(ValueError, match="not implemented"):
            estimate_s_rss(single_locus, method="invalid-method")

    def test_n_le_1_raises(self):
        """Verify sample size <= 1 raises ValueError."""
        from credtools.qc import estimate_s_rss

        lo = _make_locus(sample_size=1)
        with pytest.raises(ValueError, match="n must be greater than 1"):
            estimate_s_rss(lo)

    def test_fullrank_partialmle_returns_zero(self):
        """When LD matrix is full rank, null-partialmle should return 0."""
        from credtools.qc import estimate_s_rss

        lo = _make_locus(n_snps=10, seed=99)
        s = estimate_s_rss(lo, method="null-partialmle")
        assert s == 0

    def test_precomputed_eigens(self, single_locus):
        """Verify estimate_s_rss works with precomputed eigendecomposition."""
        from credtools.qc import estimate_s_rss, get_eigen

        eigens = get_eigen(single_locus.ld.r)
        s = estimate_s_rss(single_locus, eigvens=eigens)
        assert isinstance(s, (float, np.floating))


# ===================================================================
# TestKrigingRss
# ===================================================================


class TestKrigingRss:
    """Tests for kriging_rss()."""

    def test_returns_dataframe(self, single_locus):
        """Verify kriging_rss returns a pandas DataFrame."""
        from credtools.qc import kriging_rss

        result = kriging_rss(single_locus)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, single_locus):
        """Verify the result contains expected column names."""
        from credtools.qc import kriging_rss

        result = kriging_rss(single_locus)
        expected_cols = {"SNPID", "z", "condmean", "condvar", "z_std_diff", "logLR"}
        assert expected_cols.issubset(set(result.columns))

    def test_condvar_positive(self, single_locus):
        """Verify conditional variance values are all positive."""
        from credtools.qc import kriging_rss

        result = kriging_rss(single_locus)
        assert (result["condvar"] > 0).all()

    def test_auto_estimate_s(self, single_locus):
        """Verify kriging_rss works when s is auto-estimated (s=None)."""
        from credtools.qc import kriging_rss

        result = kriging_rss(single_locus, s=None)
        assert len(result) > 0

    def test_n_le_1_raises(self):
        """Verify sample size <= 1 raises ValueError."""
        from credtools.qc import kriging_rss

        lo = _make_locus(sample_size=1)
        with pytest.raises(ValueError, match="n must be greater than 1"):
            kriging_rss(lo)

    def test_row_count_matches_snps(self, single_locus):
        """Verify the number of rows matches the number of intersected SNPs."""
        from credtools.locus import intersect_sumstat_ld
        from credtools.qc import kriging_rss

        lo = intersect_sumstat_ld(single_locus.copy())
        result = kriging_rss(single_locus)
        assert len(result) == len(lo.sumstats)


# ===================================================================
# TestComputeDentistS
# ===================================================================


class TestComputeDentistS:
    """Tests for compute_dentist_s()."""

    def test_returns_dataframe(self, single_locus):
        """Verify compute_dentist_s returns a pandas DataFrame."""
        from credtools.qc import compute_dentist_s

        result = compute_dentist_s(single_locus)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, single_locus):
        """Verify the result contains the expected column set."""
        from credtools.qc import compute_dentist_s

        result = compute_dentist_s(single_locus)
        expected_cols = {"SNPID", "t_dentist_s", "-log10p_dentist_s", "r2"}
        assert expected_cols == set(result.columns)

    def test_lead_snp_t_is_nan(self, single_locus):
        """Verify the lead SNP has NaN for the DENTIST-S t statistic."""
        from credtools.qc import compute_dentist_s

        result = compute_dentist_s(single_locus)
        assert result["t_dentist_s"].isna().sum() == 1

    def test_r2_range(self, single_locus):
        """Verify r-squared values are between 0 and 1."""
        from credtools.qc import compute_dentist_s

        result = compute_dentist_s(single_locus)
        assert (result["r2"] >= 0).all()
        assert (result["r2"] <= 1.0 + 1e-6).all()


# ===================================================================
# TestCompareMaf
# ===================================================================


class TestCompareMaf:
    """Tests for compare_maf()."""

    def test_with_af2_returns_result(self, single_locus_with_af2):
        """Verify compare_maf returns non-empty result when AF2 column exists."""
        from credtools.qc import compare_maf

        result = compare_maf(single_locus_with_af2)
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_without_af2_returns_empty(self, single_locus):
        """Verify compare_maf returns empty DataFrame when AF2 column is missing."""
        from credtools.qc import compare_maf

        result = compare_maf(single_locus)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_expected_columns(self, single_locus_with_af2):
        """Verify the result contains MAF_sumstats, MAF_ld, and SNPID columns."""
        from credtools.qc import compare_maf

        result = compare_maf(single_locus_with_af2)
        assert "MAF_sumstats" in result.columns
        assert "MAF_ld" in result.columns
        assert ColName.SNPID in result.columns

    def test_maf_range(self, single_locus_with_af2):
        """Verify MAF values from LD reference are between 0 and 0.5."""
        from credtools.qc import compare_maf

        result = compare_maf(single_locus_with_af2)
        assert (result["MAF_ld"] >= 0).all()
        assert (result["MAF_ld"] <= 0.5 + 1e-6).all()


# ===================================================================
# TestSnpMissingness
# ===================================================================


class TestSnpMissingness:
    """Tests for snp_missingness()."""

    def test_column_name_format(self, multi_cohort_locus_set):
        """Verify column names follow the popu_cohort format."""
        from credtools.qc import snp_missingness

        result = snp_missingness(multi_cohort_locus_set)
        for col in result.columns:
            assert "_" in col  # format is "popu_cohort"

    def test_values_are_binary(self, multi_cohort_locus_set):
        """Verify all values in the missingness matrix are 0 or 1."""
        from credtools.qc import snp_missingness

        result = snp_missingness(multi_cohort_locus_set)
        assert set(result.values.flatten()).issubset({0, 0.0, 1, 1.0})

    def test_shared_snps_all_one(self, multi_cohort_locus_set):
        """SNPs shared across all cohorts should have value 1 in every column."""
        from credtools.qc import snp_missingness

        result = snp_missingness(multi_cohort_locus_set)
        all_present = result[result.sum(axis=1) == len(result.columns)]
        assert (all_present == 1).all().all()

    def test_single_locus(self, single_locus_set):
        """Verify a single-locus set produces one column with all ones."""
        from credtools.qc import snp_missingness

        result = snp_missingness(single_locus_set)
        assert result.shape[1] == 1
        assert (result == 1).all().all()


# ===================================================================
# TestLd4thMoment
# ===================================================================


class TestLd4thMoment:
    """Tests for ld_4th_moment()."""

    def test_column_name_format(self, multi_cohort_locus_set):
        """Verify column names follow the popu_cohort format."""
        from credtools.qc import ld_4th_moment

        result = ld_4th_moment(multi_cohort_locus_set)
        for col in result.columns:
            assert "_" in col  # format is "popu_cohort"

    def test_values_non_negative(self, multi_cohort_locus_set):
        """Verify LD 4th moment values are non-negative."""
        from credtools.qc import ld_4th_moment

        result = ld_4th_moment(multi_cohort_locus_set)
        # r^4 sums minus 1 should be >= 0 since diagonal contributes 1
        assert (result >= -1e-6).all().all()

    def test_index_is_snpid(self, multi_cohort_locus_set):
        """Verify the result index is named SNPID."""
        from credtools.qc import ld_4th_moment

        result = ld_4th_moment(multi_cohort_locus_set)
        assert result.index.name == ColName.SNPID

    def test_returns_dataframe(self, multi_cohort_locus_set):
        """Verify ld_4th_moment returns a pandas DataFrame."""
        from credtools.qc import ld_4th_moment

        result = ld_4th_moment(multi_cohort_locus_set)
        assert isinstance(result, pd.DataFrame)


# ===================================================================
# TestLdDecay
# ===================================================================


class TestLdDecay:
    """Tests for ld_decay()."""

    def test_returns_dataframe(self, multi_cohort_locus_set):
        """Verify ld_decay returns a pandas DataFrame."""
        from credtools.qc import ld_decay

        result = ld_decay(multi_cohort_locus_set)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, multi_cohort_locus_set):
        """Verify the result contains expected column names."""
        from credtools.qc import ld_decay

        result = ld_decay(multi_cohort_locus_set)
        expected_cols = {"distance_kb", "r2_avg", "decay_rate", "cohort"}
        assert expected_cols == set(result.columns)

    def test_r2_range(self, multi_cohort_locus_set):
        """Verify average r-squared values are non-negative."""
        from credtools.qc import ld_decay

        result = ld_decay(multi_cohort_locus_set)
        assert (result["r2_avg"] >= 0).all()

    def test_multi_cohort_label(self, multi_cohort_locus_set):
        """Verify each cohort has a distinct label in the result."""
        from credtools.qc import ld_decay

        result = ld_decay(multi_cohort_locus_set)
        assert result["cohort"].nunique() == len(multi_cohort_locus_set.loci)


# ===================================================================
# TestCochranQ
# ===================================================================


class TestCochranQ:
    """Tests for cochran_q()."""

    def test_returns_dataframe(self, multi_cohort_locus_set):
        """Verify cochran_q returns a pandas DataFrame."""
        from credtools.qc import cochran_q

        result = cochran_q(multi_cohort_locus_set)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, multi_cohort_locus_set):
        """Verify the result contains Q, Q_pvalue, and I_squared columns."""
        from credtools.qc import cochran_q

        result = cochran_q(multi_cohort_locus_set)
        expected_cols = {"Q", "Q_pvalue", "I_squared"}
        assert expected_cols == set(result.columns)

    def test_q_pvalue_range(self, multi_cohort_locus_set):
        """Verify Q p-values are between 0 and 1."""
        from credtools.qc import cochran_q

        result = cochran_q(multi_cohort_locus_set)
        assert (result["Q_pvalue"] >= 0).all()
        assert (result["Q_pvalue"] <= 1).all()

    def test_i_squared_non_negative(self, multi_cohort_locus_set):
        """Verify I-squared values are non-negative."""
        from credtools.qc import cochran_q

        result = cochran_q(multi_cohort_locus_set)
        assert (result["I_squared"] >= 0).all()

    def test_index_is_snpid(self, multi_cohort_locus_set):
        """Verify the result index is named SNPID."""
        from credtools.qc import cochran_q

        result = cochran_q(multi_cohort_locus_set)
        assert result.index.name == ColName.SNPID


# ===================================================================
# TestLocusQc
# ===================================================================


class TestLocusQc:
    """Tests for locus_qc()."""

    def test_returns_three_keys(self, single_locus_set):
        """Verify locus_qc returns a dict with exactly three expected keys."""
        from credtools.qc import locus_qc

        result = locus_qc(single_locus_set)
        assert "expected_z" in result
        assert "dentist_s" in result
        assert "compare_maf" in result
        assert len(result) == 3

    def test_no_heterogeneity_keys(self, single_locus_set):
        """Verify single-cohort QC does not include heterogeneity metrics."""
        from credtools.qc import locus_qc

        result = locus_qc(single_locus_set)
        assert "ld_4th_moment" not in result
        assert "ld_decay" not in result
        assert "cochran_q" not in result
        assert "snp_missingness" not in result

    def test_expected_z_has_cohort_column(self, single_locus_set):
        """Verify expected_z DataFrame contains a cohort column."""
        from credtools.qc import locus_qc

        result = locus_qc(single_locus_set)
        assert "cohort" in result["expected_z"].columns

    def test_expected_z_has_lambda_s_column(self, single_locus_set):
        """Verify expected_z DataFrame contains a lambda_s column."""
        from credtools.qc import locus_qc

        result = locus_qc(single_locus_set)
        assert "lambda_s" in result["expected_z"].columns

    def test_saves_files(self, tmp_path, single_locus_set):
        """Verify locus_qc saves output files to the specified directory."""
        from credtools.qc import locus_qc

        out_dir = str(tmp_path / "qc_out")
        locus_qc(single_locus_set, out_dir=out_dir)
        assert os.path.isdir(out_dir)
        for name in ["expected_z", "dentist_s", "compare_maf"]:
            assert os.path.exists(os.path.join(out_dir, f"{name}.txt.gz"))

    def test_multi_cohort(self, multi_cohort_locus_set):
        """Verify multi-cohort QC produces results for all cohorts."""
        from credtools.qc import locus_qc

        result = locus_qc(multi_cohort_locus_set)
        cohorts = result["expected_z"]["cohort"].unique()
        assert len(cohorts) == len(multi_cohort_locus_set.loci)


# ===================================================================
# TestIdentifyOutliers
# ===================================================================


class TestIdentifyOutliers:
    """Tests for identify_outliers()."""

    def test_relaxed_threshold_no_outliers(self, precomputed_qc_metrics):
        """Verify very relaxed thresholds produce no outliers."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers(
            precomputed_qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=100,
            z_threshold=100,
            z_std_diff_threshold=100,
            dentist_s_pvalue_threshold=100,
        )
        assert len(outliers) == 0

    def test_strict_threshold_finds_outliers(self, precomputed_qc_metrics):
        """Verify very strict thresholds identify some outliers."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers(
            precomputed_qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=0,
            z_threshold=0,
            z_std_diff_threshold=0,
            r_threshold=0,
            dentist_s_pvalue_threshold=0,
            dentist_s_r2_threshold=0,
        )
        assert len(outliers) > 0

    def test_nonexistent_cohort_returns_empty(self, precomputed_qc_metrics):
        """Verify a nonexistent cohort name returns an empty list."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers(
            precomputed_qc_metrics, cohort="NONEXISTENT_COHORT"
        )
        assert outliers == []

    def test_empty_metrics_returns_empty(self):
        """Verify empty QC metrics dict returns an empty list."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers({}, cohort="EUR_UKB")
        assert outliers == []

    def test_returns_list_of_strings(self, precomputed_qc_metrics):
        """Verify outlier SNP IDs are returned as a list of strings."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers(
            precomputed_qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=0,
            z_threshold=0,
            z_std_diff_threshold=0,
            r_threshold=0,
            dentist_s_pvalue_threshold=0,
            dentist_s_r2_threshold=0,
        )
        assert all(isinstance(s, str) for s in outliers)


# ===================================================================
# TestRemoveSnpsFromLocus
# ===================================================================


class TestRemoveSnpsFromLocus:
    """Tests for remove_snps_from_locus()."""

    def test_removes_correct_snps(self, single_locus):
        """Verify specified SNPs are removed from the locus sumstats."""
        from credtools.qc import remove_snps_from_locus

        snps_to_remove = [single_locus.sumstats[ColName.SNPID].iloc[0]]
        result = remove_snps_from_locus(single_locus, snps_to_remove)
        assert snps_to_remove[0] not in result.sumstats[ColName.SNPID].values

    def test_empty_list_returns_copy(self, single_locus):
        """Verify an empty removal list returns a locus with the same SNP count."""
        from credtools.qc import remove_snps_from_locus

        result = remove_snps_from_locus(single_locus, [])
        assert len(result.sumstats) == len(single_locus.sumstats)

    def test_ld_dimensions_shrink(self, single_locus):
        """Verify LD matrix dimensions shrink after removing SNPs."""
        from credtools.qc import remove_snps_from_locus

        snps_to_remove = [single_locus.sumstats[ColName.SNPID].iloc[0]]
        result = remove_snps_from_locus(single_locus, snps_to_remove)
        assert result.ld.r.shape[0] == single_locus.ld.r.shape[0] - 1
        assert result.ld.r.shape[1] == single_locus.ld.r.shape[1] - 1

    def test_sumstats_dimensions_shrink(self, single_locus):
        """Verify sumstats row count decreases by the number of removed SNPs."""
        from credtools.qc import remove_snps_from_locus

        snps_to_remove = list(single_locus.sumstats[ColName.SNPID].iloc[:3])
        result = remove_snps_from_locus(single_locus, snps_to_remove)
        assert len(result.sumstats) == len(single_locus.sumstats) - 3

    def test_metadata_preserved(self, single_locus):
        """Verify locus metadata is preserved after SNP removal."""
        from credtools.qc import remove_snps_from_locus

        result = remove_snps_from_locus(
            single_locus, [single_locus.sumstats[ColName.SNPID].iloc[0]]
        )
        assert result.popu == single_locus.popu
        assert result.cohort == single_locus.cohort
        assert result.sample_size == single_locus.sample_size


# ===================================================================
# TestSaveCleanedLocus
# ===================================================================


class TestSaveCleanedLocus:
    """Tests for save_cleaned_locus()."""

    def test_creates_directory(self, tmp_path, single_locus):
        """Verify save_cleaned_locus creates the output directory."""
        from credtools.qc import save_cleaned_locus

        out_dir = str(tmp_path / "nested" / "dir")
        save_cleaned_locus(single_locus, out_dir, "test_prefix")
        assert os.path.isdir(out_dir)

    def test_saves_sumstats(self, tmp_path, single_locus):
        """Verify a sumstats file is saved with the correct prefix."""
        from credtools.qc import save_cleaned_locus

        out_dir = str(tmp_path / "cleaned")
        save_cleaned_locus(single_locus, out_dir, "test_prefix")
        assert os.path.exists(os.path.join(out_dir, "test_prefix.sumstats.gz"))

    def test_saves_ld_npz(self, tmp_path, single_locus):
        """Verify an LD matrix npz file is saved with the correct prefix."""
        from credtools.qc import save_cleaned_locus

        out_dir = str(tmp_path / "cleaned")
        save_cleaned_locus(single_locus, out_dir, "test_prefix")
        assert os.path.exists(os.path.join(out_dir, "test_prefix.ld.npz"))

    def test_saves_ldmap(self, tmp_path, single_locus):
        """Verify an LD map file is saved with the correct prefix."""
        from credtools.qc import save_cleaned_locus

        out_dir = str(tmp_path / "cleaned")
        save_cleaned_locus(single_locus, out_dir, "test_prefix")
        assert os.path.exists(os.path.join(out_dir, "test_prefix.ldmap.gz"))


# ===================================================================
# TestRemoveOutliersAndRerunQc
# ===================================================================


class TestRemoveOutliersAndRerunQc:
    """Tests for remove_outliers_and_rerun_qc()."""

    def test_returns_three_element_tuple(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify the function returns a tuple of (LocusSet, dict, DataFrame)."""
        from credtools.qc import remove_outliers_and_rerun_qc

        cleaned_ls, cleaned_qc, summary = remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        assert isinstance(cleaned_ls, LocusSet)
        assert isinstance(cleaned_qc, dict)
        assert isinstance(summary, pd.DataFrame)

    def test_loci_count_unchanged(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify the number of loci is unchanged after outlier removal."""
        from credtools.qc import remove_outliers_and_rerun_qc

        cleaned_ls, _, _ = remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        assert len(cleaned_ls.loci) == len(single_locus_set.loci)

    def test_outlier_summary_columns(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify the outlier summary DataFrame contains expected columns."""
        from credtools.qc import remove_outliers_and_rerun_qc

        _, _, summary = remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        expected_cols = {
            "locus_id",
            "popu",
            "cohort",
            "original_snps",
            "outliers_removed",
            "cleaned_snps",
            "retention_rate",
        }
        assert expected_cols.issubset(set(summary.columns))

    def test_saves_cleaned_files(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify cleaned locus files are saved to the output directory."""
        from credtools.qc import remove_outliers_and_rerun_qc

        remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        cleaned_dir = os.path.join(str(tmp_path), "cleaned", "test_locus")
        assert os.path.isdir(cleaned_dir)


# ===================================================================
# TestLocusQcSummary
# ===================================================================


class TestLocusQcSummary:
    """Tests for locus_qc_summary()."""

    def test_returns_dataframe(self, precomputed_qc_metrics):
        """Verify locus_qc_summary returns a pandas DataFrame."""
        from credtools.qc import locus_qc_summary

        result = locus_qc_summary(precomputed_qc_metrics)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, precomputed_qc_metrics):
        """Verify the summary contains expected QC metric columns."""
        from credtools.qc import locus_qc_summary

        result = locus_qc_summary(precomputed_qc_metrics)
        expected_cols = {
            "popu",
            "cohort",
            "n_snps",
            "n_1e-5",
            "n_5e-8",
            "maf_corr",
            "lambda_s",
            "n_lambda_s_outlier",
            "n_dentist_s_outlier",
        }
        assert expected_cols.issubset(set(result.columns))

    def test_row_count_equals_cohort_count(self, precomputed_qc_metrics):
        """Verify the number of rows equals the number of unique cohorts."""
        from credtools.qc import locus_qc_summary

        result = locus_qc_summary(precomputed_qc_metrics)
        n_cohorts = precomputed_qc_metrics["expected_z"]["cohort"].nunique()
        assert len(result) == n_cohorts

    def test_multi_cohort_summary(self, multi_cohort_locus_set):
        """Verify multi-cohort summary has one row per cohort."""
        from credtools.qc import locus_qc, locus_qc_summary

        qc = locus_qc(multi_cohort_locus_set)
        result = locus_qc_summary(qc)
        assert len(result) == len(multi_cohort_locus_set.loci)

    def test_empty_metrics_returns_empty(self):
        """Verify empty QC metrics dict returns an empty DataFrame."""
        from credtools.qc import locus_qc_summary

        result = locus_qc_summary({})
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0


# ===================================================================
# TestQcLocusCli
# ===================================================================


class TestQcLocusCli:
    """Tests for qc_locus_cli()."""

    def test_mock_returns_four_tuple(self, tmp_path, single_locus_set):
        """Verify qc_locus_cli returns a 4-element tuple with expected types."""
        from credtools.qc import qc_locus_cli

        locus_info = pd.DataFrame(
            {
                "locus_id": ["test_locus"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "prefix": ["dummy"],
            }
        )

        with patch("credtools.qc.load_locus_set", return_value=single_locus_set):
            args = (
                "test_locus",
                locus_info,
                str(tmp_path),
                False,  # remove_outlier
                2.0,  # logLR_threshold
                2.0,  # z_threshold
                3.0,  # z_std_diff_threshold
                0.8,  # r_threshold
                4.0,  # dentist_s_pvalue_threshold
                0.6,  # dentist_s_r2_threshold
            )
            result = qc_locus_cli(args)

        assert len(result) == 4
        locus_id, summary, outlier_summary, cleaned_summary = result
        assert locus_id == "test_locus"
        assert isinstance(summary, pd.DataFrame)

    def test_saves_output_files(self, tmp_path, single_locus_set):
        """Verify qc_locus_cli saves QC output files to the locus directory."""
        from credtools.qc import qc_locus_cli

        locus_info = pd.DataFrame(
            {
                "locus_id": ["test_locus"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "prefix": ["dummy"],
            }
        )

        with patch("credtools.qc.load_locus_set", return_value=single_locus_set):
            args = (
                "test_locus",
                locus_info,
                str(tmp_path),
                False,
                2.0,
                2.0,
                3.0,
                0.8,
                4.0,
                0.6,
            )
            qc_locus_cli(args)

        locus_dir = os.path.join(str(tmp_path), "test_locus")
        assert os.path.isdir(locus_dir)
        assert os.path.exists(os.path.join(locus_dir, "qc.txt.gz"))


# ===================================================================
# TestSafeQcLocusCli
# ===================================================================


class TestSafeQcLocusCli:
    """Tests for safe_qc_locus_cli()."""

    def test_exception_captured(self):
        """Verify exceptions are captured and returned as the error element."""
        from credtools.qc import safe_qc_locus_cli

        locus_info = pd.DataFrame({"locus_id": ["test"]})
        args = (
            "test_locus",
            locus_info,
            "/nonexistent/path",
            False,
            2.0,
            2.0,
            3.0,
            0.8,
            4.0,
            0.6,
        )

        result = safe_qc_locus_cli(args)
        assert len(result) == 5
        locus_id, summary, outlier_summary, cleaned_summary, error = result
        assert locus_id == "test_locus"
        assert error is not None  # Should have captured the error

    def test_success_has_no_error(self, tmp_path, single_locus_set):
        """Verify successful execution returns None as the error element."""
        from credtools.qc import safe_qc_locus_cli

        locus_info = pd.DataFrame(
            {
                "locus_id": ["test_locus"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "prefix": ["dummy"],
            }
        )

        with patch("credtools.qc.load_locus_set", return_value=single_locus_set):
            args = (
                "test_locus",
                locus_info,
                str(tmp_path),
                False,
                2.0,
                2.0,
                3.0,
                0.8,
                4.0,
                0.6,
            )
            result = safe_qc_locus_cli(args)

        assert len(result) == 5
        assert result[4] is None  # error should be None


# ===================================================================
# TestLociQc
# ===================================================================


class TestLociQc:
    """Tests for loci_qc()."""

    def test_threads_less_than_1_raises(self):
        """Verify threads=0 raises ValueError."""
        from credtools.qc import loci_qc

        with pytest.raises(ValueError, match="threads must be a positive integer"):
            loci_qc("dummy.txt", "dummy_out", threads=0)

    def test_creates_output_dir(self, tmp_path):
        """Verify loci_qc creates the output directory."""
        from credtools.qc import loci_qc

        out_dir = str(tmp_path / "new_qc_output")

        # Create a minimal valid input file
        input_df = pd.DataFrame(
            {
                "locus_id": ["locus1"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "prefix": ["dummy"],
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "n_snps": [20],
            }
        )
        input_file = str(tmp_path / "input.txt")
        input_df.to_csv(input_file, sep="\t", index=False)

        # Mock the pool to avoid real multiprocessing
        mock_result = ("locus1", pd.DataFrame(), None, None, None)

        with patch("credtools.qc.Pool") as mock_pool_cls:
            mock_pool = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_pool)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
            mock_pool.imap_unordered.return_value = [mock_result]
            loci_qc(input_file, out_dir, threads=1)

        assert os.path.isdir(out_dir)
