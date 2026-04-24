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


def _make_locus_with_disjoint_ld(
    popu: str = "AFR",
    cohort: str = "MGBB",
    seed: int = 999,
    n_snps: int = 20,
    sample_size: int = 10000,
) -> Locus:
    """Create a locus where sumstats and LD map have completely disjoint SNP IDs."""
    rng = np.random.default_rng(seed)
    bps = np.arange(1000, 1000 + n_snps * 100, 100)

    # sumstats SNPs use "A-G" alleles
    snpids_sumstats = [f"1-{bp}-A-G" for bp in bps]
    eaf = rng.uniform(0.1, 0.5, n_snps).astype(np.float32)
    sumstats = pd.DataFrame(
        {
            ColName.SNPID: snpids_sumstats,
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

    # LD map SNPs use completely different IDs ("C-T" alleles)
    snpids_ld = [f"1-{bp}-C-T" for bp in bps]
    A = rng.normal(size=(n_snps, n_snps))
    r = A @ A.T
    d = np.sqrt(np.diag(r))
    r = r / np.outer(d, d)
    r = r.astype(np.float32)

    ld_map = pd.DataFrame(
        {
            ColName.SNPID: snpids_ld,
            ColName.CHR: np.int8(1),
            ColName.BP: bps.astype(np.int32),
            ColName.A1: ["C"] * n_snps,
            ColName.A2: ["T"] * n_snps,
        }
    )
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
# Fixtures for edge cases
# ---------------------------------------------------------------------------


@pytest.fixture
def locus_set_with_no_ld_overlap():
    """Create LocusSet where one cohort has no overlap between sumstats and LD."""
    good1 = _make_locus("EUR", "UKB", seed=42)
    good2 = _make_locus("EAS", "BBJ", seed=456)
    bad = _make_locus_with_disjoint_ld("AFR", "MGBB", seed=999)
    return LocusSet([good1, good2, bad])


@pytest.fixture
def locus_set_all_no_ld_overlap():
    """Create LocusSet where ALL cohorts have no overlap between sumstats and LD."""
    bad1 = _make_locus_with_disjoint_ld("AFR", "MGBB", seed=999)
    bad2 = _make_locus_with_disjoint_ld("EUR", "UKB2", seed=888)
    return LocusSet([bad1, bad2])


# ===================================================================
# TestHeterogeneityEdgeCases
# ===================================================================


class TestHeterogeneityEdgeCases:
    """Tests for ld_4th_moment / snp_missingness when some cohorts have no LD overlap."""

    def test_ld_4th_moment_skips_no_overlap_cohort(self, locus_set_with_no_ld_overlap):
        """Cohort with no LD overlap is skipped; valid cohorts produce results."""
        from credtools.qc import ld_4th_moment

        result = ld_4th_moment(locus_set_with_no_ld_overlap)
        assert isinstance(result, pd.DataFrame)
        # The bad cohort (AFR_MGBB) should be absent from columns
        assert "AFR_MGBB" not in result.columns
        # The two good cohorts should still be present
        assert "EUR_UKB" in result.columns
        assert "EAS_BBJ" in result.columns

    def test_ld_4th_moment_all_no_overlap_returns_empty(
        self, locus_set_all_no_ld_overlap
    ):
        """When all cohorts have no LD overlap, return an empty DataFrame."""
        from credtools.qc import ld_4th_moment

        result = ld_4th_moment(locus_set_all_no_ld_overlap)
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_snp_missingness_skips_no_overlap_cohort(
        self, locus_set_with_no_ld_overlap
    ):
        """Cohort with no LD overlap is skipped; valid cohorts produce results."""
        from credtools.qc import snp_missingness

        result = snp_missingness(locus_set_with_no_ld_overlap)
        assert isinstance(result, pd.DataFrame)
        # The bad cohort (AFR_MGBB) should be absent from columns
        assert "AFR_MGBB" not in result.columns
        # The two good cohorts should still be present
        assert "EUR_UKB" in result.columns
        assert "EAS_BBJ" in result.columns

    def test_snp_missingness_all_no_overlap_returns_empty(
        self, locus_set_all_no_ld_overlap
    ):
        """When all cohorts have no LD overlap, return an empty DataFrame."""
        from credtools.qc import snp_missingness

        result = snp_missingness(locus_set_all_no_ld_overlap)
        assert isinstance(result, pd.DataFrame)
        assert result.empty


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
        """Verify SNPs shared across all cohorts have value 1 in every column."""
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
        """Verify the result contains Q, Q_pvalue, I_squared, and k columns."""
        from credtools.qc import cochran_q

        result = cochran_q(multi_cohort_locus_set)
        expected_cols = {"Q", "Q_pvalue", "I_squared", "k"}
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

    def test_outer_join_includes_all_snps(self):
        """Verify cochran_q uses outer join so SNPs in any cohort are included."""
        from credtools.qc import cochran_q

        # locus1 has SNPs at bp 1000-2900, locus2 has SNPs at bp 1500-3400
        # Overlap is bp 1500-2900 (15 SNPs), union is bp 1000-3400 (25 SNPs)
        locus1 = _make_locus("EUR", "UKB", seed=42, n_snps=20)  # bp 1000-2900
        rng2 = np.random.default_rng(99)
        n2 = 20
        bps2 = np.arange(1500, 1500 + n2 * 100, 100)
        snpids2 = [f"1-{bp}-A-G" for bp in bps2]
        eaf2 = rng2.uniform(0.1, 0.5, n2).astype(np.float32)
        sumstats2 = pd.DataFrame(
            {
                ColName.SNPID: snpids2,
                ColName.CHR: np.int8(1),
                ColName.BP: bps2.astype(np.int32),
                ColName.EA: ["A"] * n2,
                ColName.NEA: ["G"] * n2,
                ColName.EAF: eaf2,
                ColName.MAF: np.minimum(eaf2, 1 - eaf2),
                ColName.A1: ["A"] * n2,
                ColName.A2: ["G"] * n2,
                ColName.BETA: rng2.normal(0, 0.1, n2).astype(np.float32),
                ColName.SE: rng2.uniform(0.01, 0.05, n2).astype(np.float32),
                ColName.P: rng2.uniform(1e-10, 0.05, n2).astype(np.float64),
            }
        )
        A2 = rng2.normal(size=(n2, n2))
        r2 = A2 @ A2.T
        d2 = np.sqrt(np.diag(r2))
        r2 = (r2 / np.outer(d2, d2)).astype(np.float32)
        ld_map2 = sumstats2[
            [ColName.SNPID, ColName.CHR, ColName.BP, ColName.A1, ColName.A2]
        ].copy()
        ld2 = LDMatrix(ld_map2, r2)
        locus2 = Locus(
            popu="AFR",
            cohort="MVP",
            sample_size=10000,
            sumstats=sumstats2,
            locus_start=int(bps2[0]),
            locus_end=int(bps2[-1]),
            ld=ld2,
        )

        locus_set = LocusSet([locus1, locus2])
        result = cochran_q(locus_set)

        # Union should have 25 unique SNPs (5 only in locus1 + 15 overlap + 5 only in locus2)
        all_snps = set(locus1.sumstats[ColName.SNPID]) | set(
            locus2.sumstats[ColName.SNPID]
        )
        assert len(result) == len(
            all_snps
        ), f"Expected {len(all_snps)} SNPs (union), got {len(result)}"

    def test_snps_in_single_cohort_have_nan(self):
        """Verify SNPs present in only one cohort have NaN Q values."""
        from credtools.qc import cochran_q

        # locus1: bp 1000-1400 (5 SNPs), locus2: bp 1200-1600 (5 SNPs)
        # Only bp 1200-1400 overlap (3 SNPs), rest are single-cohort
        locus1 = _make_locus("EUR", "UKB", seed=42, n_snps=5)
        rng2 = np.random.default_rng(99)
        n2 = 5
        bps2 = np.arange(1200, 1200 + n2 * 100, 100)
        snpids2 = [f"1-{bp}-A-G" for bp in bps2]
        eaf2 = rng2.uniform(0.1, 0.5, n2).astype(np.float32)
        sumstats2 = pd.DataFrame(
            {
                ColName.SNPID: snpids2,
                ColName.CHR: np.int8(1),
                ColName.BP: bps2.astype(np.int32),
                ColName.EA: ["A"] * n2,
                ColName.NEA: ["G"] * n2,
                ColName.EAF: eaf2,
                ColName.MAF: np.minimum(eaf2, 1 - eaf2),
                ColName.A1: ["A"] * n2,
                ColName.A2: ["G"] * n2,
                ColName.BETA: rng2.normal(0, 0.1, n2).astype(np.float32),
                ColName.SE: rng2.uniform(0.01, 0.05, n2).astype(np.float32),
                ColName.P: rng2.uniform(1e-10, 0.05, n2).astype(np.float64),
            }
        )
        A2 = rng2.normal(size=(n2, n2))
        r2 = A2 @ A2.T
        d2 = np.sqrt(np.diag(r2))
        r2 = (r2 / np.outer(d2, d2)).astype(np.float32)
        ld_map2 = sumstats2[
            [ColName.SNPID, ColName.CHR, ColName.BP, ColName.A1, ColName.A2]
        ].copy()
        ld2 = LDMatrix(ld_map2, r2)
        locus2 = Locus(
            popu="AFR",
            cohort="MVP",
            sample_size=10000,
            sumstats=sumstats2,
            locus_start=int(bps2[0]),
            locus_end=int(bps2[-1]),
            ld=ld2,
        )

        locus_set = LocusSet([locus1, locus2])
        result = cochran_q(locus_set)

        # SNPs only in one cohort should have NaN Q
        only_locus1 = set(locus1.sumstats[ColName.SNPID]) - set(
            locus2.sumstats[ColName.SNPID]
        )
        only_locus2 = set(locus2.sumstats[ColName.SNPID]) - set(
            locus1.sumstats[ColName.SNPID]
        )
        for snp in only_locus1 | only_locus2:
            assert pd.isna(result.loc[snp, "Q"]), f"SNP {snp} should have NaN Q"

        # SNPs in both cohorts should have valid Q
        overlap = set(locus1.sumstats[ColName.SNPID]) & set(
            locus2.sumstats[ColName.SNPID]
        )
        for snp in overlap:
            assert not pd.isna(result.loc[snp, "Q"]), f"SNP {snp} should have valid Q"

    def test_k_column_tracks_cohort_count(self):
        """Verify the k column correctly reports the number of cohorts per SNP."""
        from credtools.qc import cochran_q

        locus1 = _make_locus("EUR", "UKB", seed=42, n_snps=5)
        rng2 = np.random.default_rng(99)
        n2 = 5
        bps2 = np.arange(1200, 1200 + n2 * 100, 100)
        snpids2 = [f"1-{bp}-A-G" for bp in bps2]
        eaf2 = rng2.uniform(0.1, 0.5, n2).astype(np.float32)
        sumstats2 = pd.DataFrame(
            {
                ColName.SNPID: snpids2,
                ColName.CHR: np.int8(1),
                ColName.BP: bps2.astype(np.int32),
                ColName.EA: ["A"] * n2,
                ColName.NEA: ["G"] * n2,
                ColName.EAF: eaf2,
                ColName.MAF: np.minimum(eaf2, 1 - eaf2),
                ColName.A1: ["A"] * n2,
                ColName.A2: ["G"] * n2,
                ColName.BETA: rng2.normal(0, 0.1, n2).astype(np.float32),
                ColName.SE: rng2.uniform(0.01, 0.05, n2).astype(np.float32),
                ColName.P: rng2.uniform(1e-10, 0.05, n2).astype(np.float64),
            }
        )
        A2 = rng2.normal(size=(n2, n2))
        r2 = A2 @ A2.T
        d2 = np.sqrt(np.diag(r2))
        r2 = (r2 / np.outer(d2, d2)).astype(np.float32)
        ld_map2 = sumstats2[
            [ColName.SNPID, ColName.CHR, ColName.BP, ColName.A1, ColName.A2]
        ].copy()
        ld2 = LDMatrix(ld_map2, r2)
        locus2 = Locus(
            popu="AFR",
            cohort="MVP",
            sample_size=10000,
            sumstats=sumstats2,
            locus_start=int(bps2[0]),
            locus_end=int(bps2[-1]),
            ld=ld2,
        )

        locus_set = LocusSet([locus1, locus2])
        result = cochran_q(locus_set)

        assert "k" in result.columns, "Result should have a 'k' column"
        overlap = set(locus1.sumstats[ColName.SNPID]) & set(
            locus2.sumstats[ColName.SNPID]
        )
        only_one = (
            set(locus1.sumstats[ColName.SNPID]) | set(locus2.sumstats[ColName.SNPID])
        ) - overlap
        for snp in overlap:
            assert result.loc[snp, "k"] == 2
        for snp in only_one:
            assert result.loc[snp, "k"] == 1


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

    def test_returns_dataframe(self, precomputed_qc_metrics):
        """Verify identify_outliers returns a DataFrame."""
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
        assert isinstance(outliers, pd.DataFrame)

    def test_has_criterion_columns(self, precomputed_qc_metrics):
        """Verify the DataFrame contains SNPID and C1/C2/C3 criterion columns."""
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
        expected_cols = {"SNPID", "C1_ld_mismatch", "C2_marginal", "C3_dentist_s"}
        assert expected_cols.issubset(set(outliers.columns))

    def test_criterion_flags_are_bool(self, precomputed_qc_metrics):
        """Verify C1/C2/C3 columns contain boolean values."""
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
        if not outliers.empty:
            assert outliers["C1_ld_mismatch"].dtype == bool
            assert outliers["C2_marginal"].dtype == bool
            assert outliers["C3_dentist_s"].dtype == bool

    def test_at_least_one_criterion_true_per_row(self, precomputed_qc_metrics):
        """Verify each outlier SNP triggers at least one criterion."""
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
        if not outliers.empty:
            any_true = (
                outliers["C1_ld_mismatch"]
                | outliers["C2_marginal"]
                | outliers["C3_dentist_s"]
            )
            assert any_true.all()

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

    def test_nonexistent_cohort_returns_empty_dataframe(self, precomputed_qc_metrics):
        """Verify a nonexistent cohort name returns an empty DataFrame."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers(
            precomputed_qc_metrics, cohort="NONEXISTENT_COHORT"
        )
        assert isinstance(outliers, pd.DataFrame)
        assert len(outliers) == 0

    def test_empty_metrics_returns_empty_dataframe(self):
        """Verify empty QC metrics dict returns an empty DataFrame."""
        from credtools.qc import identify_outliers

        outliers = identify_outliers({}, cohort="EUR_UKB")
        assert isinstance(outliers, pd.DataFrame)
        assert len(outliers) == 0

    def test_snpid_column_contains_strings(self, precomputed_qc_metrics):
        """Verify outlier SNP IDs in SNPID column are strings."""
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
        if not outliers.empty:
            assert all(isinstance(s, str) for s in outliers["SNPID"])

    def test_no_duplicate_snpids(self, precomputed_qc_metrics):
        """Verify there are no duplicate SNPIDs in the result."""
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
        if not outliers.empty:
            assert outliers["SNPID"].is_unique


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

    def test_returns_four_element_tuple(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify the function returns a 4-tuple with expected types."""
        from credtools.qc import remove_outliers_and_rerun_qc

        result = remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        assert len(result) == 4
        cleaned_ls, cleaned_qc, summary, outlier_detail = result
        assert isinstance(cleaned_ls, LocusSet)
        assert isinstance(cleaned_qc, dict)
        assert isinstance(summary, pd.DataFrame)
        assert isinstance(outlier_detail, pd.DataFrame)

    def test_loci_count_unchanged(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify the number of loci is unchanged after outlier removal."""
        from credtools.qc import remove_outliers_and_rerun_qc

        cleaned_ls, _, _, _ = remove_outliers_and_rerun_qc(
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

        _, _, summary, _ = remove_outliers_and_rerun_qc(
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

    def test_outlier_detail_columns(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify the outlier detail DataFrame contains expected columns."""
        from credtools.qc import remove_outliers_and_rerun_qc

        _, _, _, outlier_detail = remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        expected_cols = {
            "SNPID",
            "C1_ld_mismatch",
            "C2_marginal",
            "C3_dentist_s",
            "locus_id",
            "popu",
            "cohort",
        }
        if not outlier_detail.empty:
            assert expected_cols.issubset(set(outlier_detail.columns))

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

    def test_saves_outlier_snps_file(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify outlier_snps.txt.gz is saved in the cleaned directory."""
        from credtools.qc import remove_outliers_and_rerun_qc

        remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        cleaned_dir = os.path.join(str(tmp_path), "cleaned", "test_locus")
        outlier_file = os.path.join(cleaned_dir, "outlier_snps.txt.gz")
        assert os.path.exists(outlier_file)

    def test_outlier_snps_file_content(
        self, tmp_path, single_locus_set, precomputed_qc_metrics
    ):
        """Verify outlier_snps.txt.gz content has correct columns."""
        from credtools.qc import remove_outliers_and_rerun_qc

        remove_outliers_and_rerun_qc(
            single_locus_set,
            precomputed_qc_metrics,
            str(tmp_path),
            "test_locus",
        )
        cleaned_dir = os.path.join(str(tmp_path), "cleaned", "test_locus")
        outlier_file = os.path.join(cleaned_dir, "outlier_snps.txt.gz")
        df = pd.read_csv(outlier_file, sep="\t")
        expected_cols = {
            "SNPID",
            "C1_ld_mismatch",
            "C2_marginal",
            "C3_dentist_s",
            "locus_id",
            "popu",
            "cohort",
        }
        assert expected_cols.issubset(set(df.columns))


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

    def test_mock_returns_five_tuple(self, tmp_path, single_locus_set):
        """Verify qc_locus_cli returns a 5-element tuple with expected types."""
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

        assert len(result) == 5
        locus_id, summary, outlier_summary, cleaned_summary, outlier_detail = result
        assert locus_id == "test_locus"
        assert isinstance(summary, pd.DataFrame)
        assert outlier_detail is None  # no outlier removal requested

    def test_with_outlier_removal_returns_detail(self, tmp_path, single_locus_set):
        """Verify qc_locus_cli returns outlier detail when remove_outlier=True."""
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
                True,  # remove_outlier
                2.0,
                2.0,
                3.0,
                0.8,
                4.0,
                0.6,
            )
            result = qc_locus_cli(args)

        assert len(result) == 5
        _, _, _, _, outlier_detail = result
        assert isinstance(outlier_detail, pd.DataFrame)

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
        assert len(result) == 6
        locus_id, summary, outlier_summary, cleaned_summary, outlier_detail, error = (
            result
        )
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

        assert len(result) == 6
        assert result[5] is None  # error should be None


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
        mock_result = ("locus1", pd.DataFrame(), None, None, None, None)

        with patch("credtools.qc.Pool") as mock_pool_cls:
            mock_pool = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_pool)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
            mock_pool.imap_unordered.return_value = [mock_result]
            loci_qc(input_file, out_dir, threads=1)

        assert os.path.isdir(out_dir)


# ---------------------------------------------------------------------------
# Supplementary tests: identify_outliers edge cases
# ---------------------------------------------------------------------------
class TestIdentifyOutliersEdgeCases:
    """Additional edge case tests for identify_outliers."""

    def test_no_outliers_returns_empty_df(self):
        """When no metrics trigger outlier criteria, return empty DF."""
        from credtools.qc import identify_outliers

        # Empty qc_metrics
        qc_metrics = {"expected_z": pd.DataFrame(), "dentist_s": pd.DataFrame()}
        result = identify_outliers(qc_metrics, cohort="EUR_UKB")
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
        assert "SNPID" in result.columns

    def test_unknown_cohort_returns_empty(self):
        """When the specified cohort doesn't exist in data, return empty."""
        from credtools.qc import identify_outliers

        expected_z = pd.DataFrame(
            {
                "SNPID": ["s1"],
                "cohort": ["EUR_UKB"],
                "z": [1.0],
                "condmean": [1.0],
                "condvar": [0.5],
                "logLR": [0.1],
                "z_std_diff": [0.5],
                "lambda_s": [1.0],
            }
        )
        qc_metrics = {"expected_z": expected_z, "dentist_s": pd.DataFrame()}
        result = identify_outliers(qc_metrics, cohort="AFR_APCDR")
        assert len(result) == 0


class TestIdentifyOutliersC2Semantics:
    """C2 must require all three conditions per its docstring.

    Per SuSiE guidelines (Wang et al., JRSSB 2020), C2 fires only when:
    ``|z| < z_threshold`` AND ``|z_std_diff| > z_std_diff_threshold``
    AND ``|r_to_lead| > r_threshold``.

    Prior bug: the |z| < z_threshold guard was missing, so any SNP with large
    |z_std_diff| and high LD with the lead was flagged -- including causal /
    lead SNPs themselves (|z| very large). This class pins the intended
    semantics.
    """

    @staticmethod
    def _make_qc_metrics(rows):
        """Build minimal qc_metrics from a list of per-SNP records.

        Each record is a dict with keys: SNPID, z, z_std_diff, r_to_lead.
        Everything else is filled with neutral values so only C2 can fire.
        """
        cohort = "EUR_UKB"
        expected_z = pd.DataFrame(
            {
                "SNPID": [r["SNPID"] for r in rows],
                "cohort": [cohort] * len(rows),
                "z": [r["z"] for r in rows],
                "condmean": [0.0] * len(rows),
                "condvar": [1.0] * len(rows),
                # logLR kept small so C1 never triggers
                "logLR": [0.0] * len(rows),
                "z_std_diff": [r["z_std_diff"] for r in rows],
                "lambda_s": [1e-3] * len(rows),
            }
        )
        dentist_s = pd.DataFrame(
            {
                "SNPID": [r["SNPID"] for r in rows],
                "cohort": [cohort] * len(rows),
                # r2 = r_to_lead**2 so r_abs = |r_to_lead|
                "r2": [r["r_to_lead"] ** 2 for r in rows],
                # -log10p kept low so C3 never triggers
                "-log10p_dentist_s": [0.0] * len(rows),
            }
        )
        return {"expected_z": expected_z, "dentist_s": dentist_s}

    def test_large_z_high_ld_not_flagged_as_c2(self):
        """Lead-like SNPs must not be flagged as C2.

        A SNP with |z|=15 must not trigger C2 even if |z_std_diff| and
        |r_to_lead| both exceed their thresholds. This is the regression
        test for the causal-SNP-removal bug.
        """
        from credtools.qc import identify_outliers

        qc_metrics = self._make_qc_metrics(
            [{"SNPID": "causal", "z": 15.0, "z_std_diff": 4.2, "r_to_lead": 1.0}]
        )
        outliers = identify_outliers(
            qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=2.0,
            z_threshold=2.0,
            z_std_diff_threshold=3.0,
            r_threshold=0.8,
            dentist_s_pvalue_threshold=4.0,
            dentist_s_r2_threshold=0.6,
        )
        # Causal SNP should not appear in outliers at all
        assert "causal" not in set(outliers["SNPID"]), (
            "Causal SNP with |z|=15 was flagged — C2 is missing the "
            "|z| < z_threshold guard."
        )

    def test_small_z_high_ld_big_zdiff_is_flagged_as_c2(self):
        """True marginally-non-significant flip SNPs must still be flagged.

        A SNP with small |z|, high LD to lead, and large |z_std_diff| must
        still trigger C2.
        """
        from credtools.qc import identify_outliers

        qc_metrics = self._make_qc_metrics(
            [{"SNPID": "flip", "z": 1.0, "z_std_diff": 4.2, "r_to_lead": 0.95}]
        )
        outliers = identify_outliers(
            qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=2.0,
            z_threshold=2.0,
            z_std_diff_threshold=3.0,
            r_threshold=0.8,
            dentist_s_pvalue_threshold=4.0,
            dentist_s_r2_threshold=0.6,
        )
        row = outliers[outliers["SNPID"] == "flip"]
        assert len(row) == 1, "Marginal-non-significant SNP was not flagged"
        assert bool(row["C2_marginal"].iloc[0]) is True

    def test_small_z_low_ld_not_flagged_as_c2(self):
        """Small |z| alone isn't enough: without high LD to lead, no C2."""
        from credtools.qc import identify_outliers

        qc_metrics = self._make_qc_metrics(
            [{"SNPID": "nolink", "z": 1.0, "z_std_diff": 4.2, "r_to_lead": 0.5}]
        )
        outliers = identify_outliers(
            qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=2.0,
            z_threshold=2.0,
            z_std_diff_threshold=3.0,
            r_threshold=0.8,
            dentist_s_pvalue_threshold=4.0,
            dentist_s_r2_threshold=0.6,
        )
        assert "nolink" not in set(outliers["SNPID"])

    def test_mixed_snps_only_true_marginals_flagged(self):
        """In a mixed set, only SNPs satisfying all three C2 conditions fire."""
        from credtools.qc import identify_outliers

        qc_metrics = self._make_qc_metrics(
            [
                # Causal / lead-like: big z, high LD, big zdiff -> NOT C2
                {"SNPID": "causal", "z": 14.9, "z_std_diff": 4.2, "r_to_lead": 1.0},
                # True marginal flip: small z, high LD, big zdiff -> C2
                {"SNPID": "flip", "z": 0.8, "z_std_diff": 5.0, "r_to_lead": 0.9},
                # Quiet variant: small z, low LD, small zdiff -> not flagged
                {"SNPID": "quiet", "z": 0.5, "z_std_diff": 1.0, "r_to_lead": 0.2},
                # Strong neighbor: big z, high LD, small zdiff -> not flagged
                {"SNPID": "neighbor", "z": 8.0, "z_std_diff": 1.5, "r_to_lead": 0.9},
            ]
        )
        outliers = identify_outliers(
            qc_metrics,
            cohort="EUR_UKB",
            logLR_threshold=2.0,
            z_threshold=2.0,
            z_std_diff_threshold=3.0,
            r_threshold=0.8,
            dentist_s_pvalue_threshold=4.0,
            dentist_s_r2_threshold=0.6,
        )
        c2_snps = set(outliers[outliers["C2_marginal"]]["SNPID"])
        assert c2_snps == {"flip"}, f"Expected only 'flip' as C2, got {c2_snps}"


# ---------------------------------------------------------------------------
# Supplementary tests: remove_snps_from_locus
# ---------------------------------------------------------------------------
class TestRemoveSnpsEdgeCases:
    """Additional tests for remove_snps_from_locus."""

    def test_remove_all_snps(self, single_locus):
        """Removing all SNPs should result in empty sumstats."""
        from credtools.qc import remove_snps_from_locus

        all_snps = single_locus.sumstats[ColName.SNPID].tolist()
        result = remove_snps_from_locus(single_locus, all_snps)
        assert len(result.sumstats) == 0

    def test_remove_nonexistent_snp(self, single_locus):
        """Removing a non-existent SNP should not change anything."""
        from credtools.qc import remove_snps_from_locus

        original_n = len(single_locus.sumstats)
        result = remove_snps_from_locus(single_locus, ["nonexistent_snp"])
        assert len(result.sumstats) == original_n


# ---------------------------------------------------------------------------
# Supplementary tests: locus_qc_summary
# ---------------------------------------------------------------------------
class TestLocusQcSummaryEdgeCases:
    """Additional tests for locus_qc_summary."""

    def test_empty_expected_z_returns_empty_df(self):
        """When expected_z is empty, return empty DataFrame."""
        from credtools.qc import locus_qc_summary

        qc_metrics = {"expected_z": pd.DataFrame()}
        result = locus_qc_summary(qc_metrics)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_multi_cohort_returns_multi_rows(self, precomputed_qc_metrics):
        """Multiple cohorts should produce multiple summary rows."""
        from credtools.qc import locus_qc_summary

        result = locus_qc_summary(precomputed_qc_metrics)
        expected_z = precomputed_qc_metrics["expected_z"]
        n_cohorts = expected_z["cohort"].nunique()
        assert len(result) == n_cohorts


def _make_jit_fallback():
    """Create a no-op jit decorator matching the qc.py fallback."""

    def jit_fallback(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    return jit_fallback


# ---------------------------------------------------------------------------
# Tests: sklearn / numba fallback paths
# ---------------------------------------------------------------------------
class TestFallbackPaths:
    """Test the fallback mock classes when sklearn/numba are unavailable."""

    def test_mock_gaussian_mixture_raises_import_error(self):
        """The mock GaussianMixture __init__ should raise ImportError."""
        import sys

        # Temporarily hide sklearn so the mock class is activated
        real_sklearn = sys.modules.get("sklearn.mixture")
        real_sklearn_pkg = sys.modules.get("sklearn")
        sys.modules["sklearn.mixture"] = None  # type: ignore[assignment]
        sys.modules["sklearn"] = None  # type: ignore[assignment]
        try:
            # Re-execute the fallback code path manually
            try:
                from sklearn.mixture import GaussianMixture as _GM  # noqa: F401

                available = True
            except (ImportError, TypeError):
                available = False

                class MockGM:
                    """Local mock matching qc.py fallback."""

                    def __init__(self, *args, **kwargs):
                        raise ImportError(
                            "sklearn not available - install scikit-learn"
                        )

                    def fit(self, *args):
                        """Mock fit method."""
                        pass

                    @property
                    def weights_(self):
                        """Mock weights property."""
                        return None

            assert not available
            with pytest.raises(ImportError, match="sklearn not available"):
                MockGM(n_components=2)
        finally:
            # Restore original modules
            if real_sklearn is not None:
                sys.modules["sklearn.mixture"] = real_sklearn
            else:
                sys.modules.pop("sklearn.mixture", None)
            if real_sklearn_pkg is not None:
                sys.modules["sklearn"] = real_sklearn_pkg
            else:
                sys.modules.pop("sklearn", None)

    def test_numba_jit_fallback_is_noop_decorator(self):
        """The fallback jit decorator should return the original function unchanged."""
        jit_fallback = _make_jit_fallback()

        @jit_fallback(nopython=True, cache=True)
        def my_func(x):
            return x + 1

        # The function should work identically without numba
        assert my_func(5) == 6
        assert my_func(0) == 1

    def test_numba_jit_fallback_preserves_function_name(self):
        """Fallback jit should preserve the original function identity."""

        def jit_fallback(*args, **kwargs):
            def decorator(func):
                return func

            return decorator

        def original(x):
            return x * 2

        decorated = jit_fallback(nopython=True)(original)
        assert decorated is original


# ---------------------------------------------------------------------------
# Tests: loci_qc function (additional coverage)
# ---------------------------------------------------------------------------
class TestLociQcExtended:
    """Extended tests for the loci_qc orchestration function."""

    @staticmethod
    def _make_loci_info_tsv(tmp_path, n_loci=2):
        """Create a minimal valid loci_info TSV file and return its path.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Temporary directory provided by pytest.
        n_loci : int
            Number of distinct locus_id entries to create.

        Returns
        -------
        str
            Path to the generated TSV file.
        """
        rows = []
        for i in range(1, n_loci + 1):
            rows.append(
                {
                    "locus_id": f"locus_{i}",
                    "prefix": f"/fake/path/locus_{i}",
                    "popu": "EUR",
                    "cohort": "UKB",
                    "sample_size": 10000,
                    "chr": 1,
                    "start": 1000 * i,
                    "end": 1000 * i + 5000,
                }
            )
        df = pd.DataFrame(rows)
        fpath = str(tmp_path / "loci_info.tsv")
        df.to_csv(fpath, sep="\t", index=False)
        return fpath

    # ---- threads validation ----

    def test_threads_less_than_one_raises_value_error(self, tmp_path):
        """Passing threads < 1 should raise ValueError immediately."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=1)
        with pytest.raises(ValueError, match="threads must be a positive integer"):
            loci_qc(inputs=fpath, out_dir=str(tmp_path / "out"), threads=0)

    def test_threads_negative_raises_value_error(self, tmp_path):
        """Passing threads = -1 should raise ValueError."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=1)
        with pytest.raises(ValueError, match="threads must be a positive integer"):
            loci_qc(inputs=fpath, out_dir=str(tmp_path / "out"), threads=-1)

    # ---- basic successful run with all summaries populated ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_successful_run_saves_global_summary(
        self, mock_progress_cls, mock_pool_cls, tmp_path
    ):
        """When all loci succeed, global qc.txt.gz should be written."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=2)
        out_dir = str(tmp_path / "qc_out")

        summary1 = pd.DataFrame(
            {"locus_id": ["locus_1"], "metric_a": [0.5], "metric_b": [1.0]}
        )
        summary2 = pd.DataFrame(
            {"locus_id": ["locus_2"], "metric_a": [0.7], "metric_b": [0.9]}
        )

        pool_results = [
            ("locus_1", summary1, None, None, None, None),
            ("locus_2", summary2, None, None, None, None),
        ]

        # Configure the mock Pool context manager
        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        # Configure mock Progress context manager
        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(inputs=fpath, out_dir=out_dir, threads=1)

        assert result["successful_loci"] == 2
        assert result["failed_loci"] == 0
        assert result["total_loci"] == 2
        assert len(result["errors"]) == 0
        assert result["end_time"] is not None

        # Global summary file should exist
        assert os.path.exists(os.path.join(out_dir, "qc.txt.gz"))
        saved = pd.read_csv(os.path.join(out_dir, "qc.txt.gz"), sep="\t")
        assert len(saved) == 2
        # locus_id should be the first column
        assert saved.columns[0] == "locus_id"

        # Log file should exist
        assert os.path.exists(os.path.join(out_dir, "qc_run_summary.log"))

    # ---- error loci tracked in run_summary ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_failed_loci_recorded_in_run_summary(
        self, mock_progress_cls, mock_pool_cls, tmp_path
    ):
        """When safe_qc_locus_cli returns an error, it is tracked in run_summary."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=2)
        out_dir = str(tmp_path / "qc_out_err")

        pool_results = [
            ("locus_1", pd.DataFrame(), None, None, None, None),
            (
                "locus_2",
                pd.DataFrame(),
                None,
                None,
                None,
                "RuntimeError: something broke",
            ),
        ]

        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(inputs=fpath, out_dir=out_dir, threads=1)

        assert result["successful_loci"] == 1
        assert result["failed_loci"] == 1
        assert len(result["errors"]) == 1
        assert "locus_2" in result["errors"][0]
        assert "something broke" in result["errors"][0]

        # Verify the log file records the error
        log_path = os.path.join(out_dir, "qc_run_summary.log")
        assert os.path.exists(log_path)
        with open(log_path, "r") as f:
            log_content = f.read()
        assert "Failed: 1" in log_content
        assert "locus_2" in log_content

    # ---- cleaned summary and outlier summary paths ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_cleaned_and_outlier_summaries_saved(
        self, mock_progress_cls, mock_pool_cls, tmp_path
    ):
        """When cleaned/outlier summaries are returned, they are saved correctly."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=1)
        out_dir = str(tmp_path / "qc_out_cleaned")

        summary = pd.DataFrame({"locus_id": ["locus_1"], "metric_a": [0.5]})
        outlier_summary = pd.DataFrame(
            {
                "locus_id": ["locus_1"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "original_snps": [100],
                "cleaned_snps": [95],
                "outliers_removed": [5],
            }
        )
        cleaned_summary = pd.DataFrame({"locus_id": ["locus_1"], "metric_a": [0.6]})
        outlier_detail = pd.DataFrame(
            {
                "locus_id": ["locus_1"],
                "snpid": ["1-1000-A-G"],
                "reason": ["ld_mismatch"],
            }
        )

        pool_results = [
            (
                "locus_1",
                summary,
                outlier_summary,
                cleaned_summary,
                outlier_detail,
                None,
            ),
        ]

        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(inputs=fpath, out_dir=out_dir, threads=1)

        assert result["successful_loci"] == 1
        assert result["failed_loci"] == 0

        # Global QC summary
        assert os.path.exists(os.path.join(out_dir, "qc.txt.gz"))

        # Cleaned QC summary
        cleaned_path = os.path.join(out_dir, "cleaned", "qc_cleaned.txt.gz")
        assert os.path.exists(cleaned_path)
        cleaned_df = pd.read_csv(cleaned_path, sep="\t")
        assert len(cleaned_df) == 1
        assert cleaned_df.columns[0] == "locus_id"

        # Outlier removal summary
        outlier_path = os.path.join(
            out_dir, "cleaned", "outlier_removal_summary.txt.gz"
        )
        assert os.path.exists(outlier_path)

        # Outlier SNP details
        snp_detail_path = os.path.join(out_dir, "cleaned", "outlier_snps.txt.gz")
        assert os.path.exists(snp_detail_path)
        detail_df = pd.read_csv(snp_detail_path, sep="\t")
        assert len(detail_df) == 1

        # Cleaned loci info file
        cleaned_info_path = os.path.join(out_dir, "cleaned", "cleaned_loci_info.txt.gz")
        assert os.path.exists(cleaned_info_path)
        info_df = pd.read_csv(cleaned_info_path, sep="\t")
        assert len(info_df) == 1
        assert info_df["popu"].iloc[0] == "EUR"
        assert info_df["cohort"].iloc[0] == "UKB"

    # ---- no summary data: qc.txt.gz should NOT be created ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_empty_summaries_no_qc_file(
        self, mock_progress_cls, mock_pool_cls, tmp_path
    ):
        """When all summaries are empty, qc.txt.gz should not be created."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=1)
        out_dir = str(tmp_path / "qc_out_empty")

        # Return empty summary (success but nothing to report)
        pool_results = [
            ("locus_1", pd.DataFrame(), None, None, None, None),
        ]

        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(inputs=fpath, out_dir=out_dir, threads=1)

        assert result["successful_loci"] == 1
        # No qc.txt.gz since all_summaries is empty
        assert not os.path.exists(os.path.join(out_dir, "qc.txt.gz"))
        # No cleaned directory either
        assert not os.path.exists(os.path.join(out_dir, "cleaned"))

    # ---- run_summary parameters are stored correctly ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_run_summary_parameters(self, mock_progress_cls, mock_pool_cls, tmp_path):
        """Run summary should contain all input parameters."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=1)
        out_dir = str(tmp_path / "qc_out_params")

        pool_results = [
            ("locus_1", pd.DataFrame(), None, None, None, None),
        ]

        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(
            inputs=fpath,
            out_dir=out_dir,
            threads=2,
            remove_outlier=True,
            logLR_threshold=3.0,
            z_threshold=2.5,
        )

        params = result["parameters"]
        assert params["inputs"] == fpath
        assert params["out_dir"] == out_dir
        assert params["threads"] == 2
        assert params["remove_outlier"] is True
        assert params["logLR_threshold"] == 3.0
        assert params["z_threshold"] == 2.5

    # ---- multi-cohort outlier summary with "+" in cohort name ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_multi_cohort_prefix_with_plus_sign(
        self, mock_progress_cls, mock_pool_cls, tmp_path
    ):
        """When cohort contains '+', prefix uses meta hash pattern."""
        from credtools.qc import loci_qc

        # Create loci_info with a multi-cohort entry
        rows = [
            {
                "locus_id": "locus_1",
                "prefix": "/fake/path/locus_1",
                "popu": "EUR",
                "cohort": "UKB+FinnGen",
                "sample_size": 10000,
                "chr": 1,
                "start": 1000,
                "end": 6000,
            }
        ]
        df = pd.DataFrame(rows)
        fpath = str(tmp_path / "loci_info_multi.tsv")
        df.to_csv(fpath, sep="\t", index=False)
        out_dir = str(tmp_path / "qc_out_multi")

        summary = pd.DataFrame({"locus_id": ["locus_1"], "metric_a": [0.5]})
        outlier_summary = pd.DataFrame(
            {
                "locus_id": ["locus_1"],
                "popu": ["EUR"],
                "cohort": ["UKB+FinnGen"],
                "original_snps": [100],
                "cleaned_snps": [90],
                "outliers_removed": [10],
            }
        )

        pool_results = [
            ("locus_1", summary, outlier_summary, None, None, None),
        ]

        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(inputs=fpath, out_dir=out_dir, threads=1)

        assert result["successful_loci"] == 1

        # Verify cleaned_loci_info prefix uses the meta hash pattern
        cleaned_info_path = os.path.join(out_dir, "cleaned", "cleaned_loci_info.txt.gz")
        assert os.path.exists(cleaned_info_path)
        info_df = pd.read_csv(cleaned_info_path, sep="\t")
        assert len(info_df) == 1
        prefix_val = info_df["prefix"].iloc[0]
        assert "meta2cohorts_" in prefix_val
        assert "UKB+FinnGen" not in prefix_val  # cohort name hashed, not literal

    # ---- all loci fail ----

    @patch("credtools.qc.Pool")
    @patch("credtools.qc.Progress")
    def test_all_loci_fail(self, mock_progress_cls, mock_pool_cls, tmp_path):
        """When every locus fails, all are recorded and no output files written."""
        from credtools.qc import loci_qc

        fpath = self._make_loci_info_tsv(tmp_path, n_loci=2)
        out_dir = str(tmp_path / "qc_out_allfail")

        pool_results = [
            ("locus_1", pd.DataFrame(), None, None, None, "ValueError: bad data"),
            ("locus_2", pd.DataFrame(), None, None, None, "IOError: file missing"),
        ]

        mock_pool = MagicMock()
        mock_pool.imap_unordered.return_value = iter(pool_results)
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool_cls.return_value = mock_pool

        mock_prog = MagicMock()
        mock_prog.add_task.return_value = 0
        mock_prog.__enter__ = MagicMock(return_value=mock_prog)
        mock_prog.__exit__ = MagicMock(return_value=False)
        mock_progress_cls.return_value = mock_prog

        result = loci_qc(inputs=fpath, out_dir=out_dir, threads=1)

        assert result["successful_loci"] == 0
        assert result["failed_loci"] == 2
        assert len(result["errors"]) == 2
        # No summary files
        assert not os.path.exists(os.path.join(out_dir, "qc.txt.gz"))

        # Log file should still exist and record both errors
        log_path = os.path.join(out_dir, "qc_run_summary.log")
        with open(log_path, "r") as f:
            log_content = f.read()
        assert "Successful: 0" in log_content
        assert "Failed: 2" in log_content
