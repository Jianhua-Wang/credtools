"""Tests for meta-analysis functions in credtools.meta module."""

import os
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet
from credtools.meta import (
    heterogeneity_summary,
    meta,
    meta_all,
    meta_by_population,
    meta_lds,
    meta_loci,
    meta_locus,
    meta_sumstats,
    recover_completed_locus,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_locus(
    popu: str,
    cohort: str,
    seed: int = 42,
    sample_size: int = 10000,
    locus_start: int = 1000,
    locus_end: int = 3000,
    n_snps: int = 20,
    add_af2: bool = False,
) -> Locus:
    """Create a test locus with configurable parameters."""
    rng = np.random.default_rng(seed)
    bps = np.arange(1000, 1000 + n_snps * 100, 100)
    snpids = [f"1-{bp}-A-G" for bp in bps]

    sumstats = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: np.int8(1),
            ColName.BP: bps.astype(np.int32),
            ColName.EA: "A",
            ColName.NEA: "G",
            ColName.EAF: rng.uniform(0.1, 0.5, n_snps).astype(np.float32),
            ColName.MAF: rng.uniform(0.1, 0.5, n_snps).astype(np.float32),
            ColName.A1: "A",
            ColName.A2: "G",
            ColName.BETA: rng.normal(0, 0.1, n_snps).astype(np.float32),
            ColName.SE: rng.uniform(0.01, 0.05, n_snps).astype(np.float32),
            ColName.P: rng.uniform(1e-10, 0.05, n_snps),
        }
    )

    # Positive-definite LD matrix via random correlation
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
        locus_start=locus_start,
        locus_end=locus_end,
        ld=ld,
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def two_population_locus_set():
    """EUR/UKB + AFR/MVP, basic meta test."""
    return LocusSet(
        [
            _make_locus("EUR", "UKB", seed=42, sample_size=10000),
            _make_locus("AFR", "MVP", seed=123, sample_size=8000),
        ]
    )


@pytest.fixture
def same_population_locus_set():
    """EUR/UKB + EUR/GWAS2, same population multiple cohorts."""
    return LocusSet(
        [
            _make_locus("EUR", "UKB", seed=42, sample_size=10000),
            _make_locus("EUR", "GWAS2", seed=99, sample_size=6000),
        ]
    )


@pytest.fixture
def mixed_population_locus_set():
    """EUR/UKB + EUR/GWAS2 + AFR/MVP, for meta_by_population."""
    return LocusSet(
        [
            _make_locus("EUR", "UKB", seed=42, sample_size=10000),
            _make_locus("EUR", "GWAS2", seed=99, sample_size=6000),
            _make_locus("AFR", "MVP", seed=123, sample_size=8000),
        ]
    )


@pytest.fixture
def single_locus_set():
    """Single EUR/UKB locus."""
    return LocusSet([_make_locus("EUR", "UKB", seed=42)])


@pytest.fixture
def locus_set_with_af2():
    """Two loci with AF2 in LD map."""
    return LocusSet(
        [
            _make_locus("EUR", "UKB", seed=42, sample_size=10000, add_af2=True),
            _make_locus("AFR", "MVP", seed=123, sample_size=8000, add_af2=True),
        ]
    )


@pytest.fixture
def mismatched_start_locus_set():
    """Loci with different locus_start."""
    return LocusSet(
        [
            _make_locus("EUR", "UKB", seed=42, locus_start=1000),
            _make_locus("AFR", "MVP", seed=123, locus_start=2000),
        ]
    )


@pytest.fixture
def mismatched_end_locus_set():
    """Loci with different locus_end."""
    return LocusSet(
        [
            _make_locus("EUR", "UKB", seed=42, locus_end=3000),
            _make_locus("AFR", "MVP", seed=123, locus_end=4000),
        ]
    )


# ===========================================================================
# TestMetaSumstats
# ===========================================================================


class TestMetaSumstats:
    """Tests for meta_sumstats (lines 34-122)."""

    def test_output_columns(self, two_population_locus_set):
        """Result DataFrame should have all required columns."""
        result = meta_sumstats(two_population_locus_set)
        for col in [
            ColName.SNPID,
            ColName.BETA,
            ColName.SE,
            ColName.P,
            ColName.EAF,
            ColName.CHR,
            ColName.BP,
            ColName.EA,
            ColName.NEA,
        ]:
            assert col in result.columns, f"Missing column {col}"

    def test_snpid_union(self, two_population_locus_set):
        """Outer join should produce the union of SNPIDs."""
        loci = two_population_locus_set.loci
        expected_snps = set(loci[0].original_sumstats[ColName.SNPID]) | set(
            loci[1].sumstats[ColName.SNPID]
        )
        result = meta_sumstats(two_population_locus_set)
        assert set(result[ColName.SNPID]) == expected_snps

    def test_ivw_formula_manual(self, two_population_locus_set):
        """Verify IVW formula by hand for a shared SNP."""
        loci = two_population_locus_set.loci
        # Pick the first SNP in common
        common = set(loci[0].sumstats[ColName.SNPID]) & set(
            loci[1].sumstats[ColName.SNPID]
        )
        snp = sorted(common)[0]

        row0 = loci[0].sumstats[loci[0].sumstats[ColName.SNPID] == snp].iloc[0]
        row1 = loci[1].sumstats[loci[1].sumstats[ColName.SNPID] == snp].iloc[0]

        w0 = 1 / (row0[ColName.SE] ** 2)
        w1 = 1 / (row1[ColName.SE] ** 2)
        expected_beta = (row0[ColName.BETA] * w0 + row1[ColName.BETA] * w1) / (w0 + w1)
        expected_se = np.sqrt(1 / (w0 + w1))

        result = meta_sumstats(two_population_locus_set)
        res_row = result[result[ColName.SNPID] == snp].iloc[0]
        np.testing.assert_allclose(res_row[ColName.BETA], expected_beta, rtol=1e-4)
        np.testing.assert_allclose(res_row[ColName.SE], expected_se, rtol=1e-4)

    def test_p_values_valid(self, two_population_locus_set):
        """All P values should be in (0, 1]."""
        result = meta_sumstats(two_population_locus_set)
        assert (result[ColName.P] > 0).all()
        assert (result[ColName.P] <= 1).all()

    def test_eaf_sample_size_weighted(self, two_population_locus_set):
        """EAF should be weighted by sample size."""
        loci = two_population_locus_set.loci
        common = set(loci[0].sumstats[ColName.SNPID]) & set(
            loci[1].sumstats[ColName.SNPID]
        )
        snp = sorted(common)[0]

        n_sum = sum(loc.sample_size for loc in loci)
        w0 = loci[0].sample_size / n_sum
        w1 = loci[1].sample_size / n_sum

        eaf0 = (
            loci[0]
            .sumstats.loc[loci[0].sumstats[ColName.SNPID] == snp, ColName.EAF]
            .iloc[0]
        )
        eaf1 = (
            loci[1]
            .sumstats.loc[loci[1].sumstats[ColName.SNPID] == snp, ColName.EAF]
            .iloc[0]
        )
        expected_eaf = eaf0 * w0 + eaf1 * w1

        result = meta_sumstats(two_population_locus_set)
        res_eaf = result.loc[result[ColName.SNPID] == snp, ColName.EAF].iloc[0]
        np.testing.assert_allclose(res_eaf, expected_eaf, rtol=1e-4)

    def test_eaf_missing_cohort_not_diluted(self):
        """EAF for SNPs only in one cohort should equal that cohort's EAF, not be diluted."""
        # Cohort A has SNPs 1-3, Cohort B has SNPs 2-4
        # SNP 1 only in A, SNP 4 only in B => meta EAF should equal original EAF
        snps_a = ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G"]
        snps_b = ["1-1100-A-G", "1-1200-A-G", "1-1300-A-G"]

        def _make_ss(snpids, eafs, betas, ses):
            bps = [int(s.split("-")[1]) for s in snpids]
            return pd.DataFrame(
                {
                    ColName.SNPID: snpids,
                    ColName.CHR: np.int8(1),
                    ColName.BP: np.array(bps, dtype=np.int32),
                    ColName.EA: "A",
                    ColName.NEA: "G",
                    ColName.EAF: np.array(eafs, dtype=np.float32),
                    ColName.MAF: np.array(eafs, dtype=np.float32),
                    ColName.A1: "A",
                    ColName.A2: "G",
                    ColName.BETA: np.array(betas, dtype=np.float32),
                    ColName.SE: np.array(ses, dtype=np.float32),
                    ColName.P: [1e-5] * len(snpids),
                }
            )

        ss_a = _make_ss(snps_a, [0.3, 0.4, 0.5], [0.1, 0.2, 0.3], [0.01, 0.02, 0.03])
        ss_b = _make_ss(snps_b, [0.35, 0.45, 0.6], [0.15, 0.25, 0.35], [0.015, 0.025, 0.035])

        def _make_ld(snpids):
            n = len(snpids)
            bps = [int(s.split("-")[1]) for s in snpids]
            ld_map = pd.DataFrame(
                {
                    ColName.SNPID: snpids,
                    ColName.CHR: np.int8(1),
                    ColName.BP: np.array(bps, dtype=np.int32),
                    ColName.A1: "A",
                    ColName.A2: "G",
                }
            )
            return LDMatrix(ld_map, np.eye(n, dtype=np.float32))

        locus_a = Locus("EUR", "UKB", 10000, ss_a, 900, 3500, ld=_make_ld(snps_a))
        locus_b = Locus("AFR", "MVP", 5000, ss_b, 900, 3500, ld=_make_ld(snps_b))
        ls = LocusSet([locus_a, locus_b])

        result = meta_sumstats(ls)

        # SNP "1-1000-A-G" only in cohort A => meta EAF should be 0.3, not 0.3 * 10000/15000
        snp_only_a = result[result[ColName.SNPID] == "1-1000-A-G"].iloc[0]
        np.testing.assert_allclose(snp_only_a[ColName.EAF], 0.3, atol=1e-4)

        # SNP "1-1300-A-G" only in cohort B => meta EAF should be 0.6, not 0.6 * 5000/15000
        snp_only_b = result[result[ColName.SNPID] == "1-1300-A-G"].iloc[0]
        np.testing.assert_allclose(snp_only_b[ColName.EAF], 0.6, atol=1e-4)

        # Shared SNP "1-1100-A-G" => weighted average
        expected = 0.4 * (10000 / 15000) + 0.35 * (5000 / 15000)
        snp_shared = result[result[ColName.SNPID] == "1-1100-A-G"].iloc[0]
        np.testing.assert_allclose(snp_shared[ColName.EAF], expected, atol=1e-4)

    def test_meta_se_smaller(self, same_population_locus_set):
        """Meta SE should be <= the smallest individual SE for shared SNPs."""
        loci = same_population_locus_set.loci
        common = set(loci[0].sumstats[ColName.SNPID]) & set(
            loci[1].sumstats[ColName.SNPID]
        )
        result = meta_sumstats(same_population_locus_set)

        for snp in sorted(common)[:5]:
            se0 = (
                loci[0]
                .sumstats.loc[loci[0].sumstats[ColName.SNPID] == snp, ColName.SE]
                .iloc[0]
            )
            se1 = (
                loci[1]
                .sumstats.loc[loci[1].sumstats[ColName.SNPID] == snp, ColName.SE]
                .iloc[0]
            )
            meta_se = result.loc[result[ColName.SNPID] == snp, ColName.SE].iloc[0]
            assert meta_se <= min(float(se0), float(se1)) + 1e-6


# ===========================================================================
# TestMetaLds
# ===========================================================================


class TestMetaLds:
    """Tests for meta_lds (lines 125-206)."""

    def test_returns_ldmatrix(self, two_population_locus_set):
        """Should return an LDMatrix."""
        result = meta_lds(two_population_locus_set)
        assert isinstance(result, LDMatrix)

    def test_dimension_matches_snp_union(self, two_population_locus_set):
        """Matrix dimension should equal the union of SNPs."""
        loci = two_population_locus_set.loci
        all_snps = set()
        for loc in loci:
            all_snps.update(loc.ld.map[ColName.SNPID].values)
        result = meta_lds(two_population_locus_set)
        assert result.r.shape[0] == len(all_snps)
        assert result.r.shape[1] == len(all_snps)

    def test_diagonal_close_to_one(self, two_population_locus_set):
        """Diagonal should be close to 1."""
        result = meta_lds(two_population_locus_set)
        np.testing.assert_allclose(np.diag(result.r), 1.0, atol=1e-4)

    def test_symmetric(self, two_population_locus_set):
        """Matrix should be symmetric."""
        result = meta_lds(two_population_locus_set)
        np.testing.assert_allclose(result.r, result.r.T, atol=1e-6)

    def test_equal_sample_size_simple_average(self):
        """When sample sizes are equal, result is simple average."""
        ls = LocusSet(
            [
                _make_locus("EUR", "UKB", seed=42, sample_size=5000),
                _make_locus("AFR", "MVP", seed=42, sample_size=5000),
            ]
        )
        # Same seed => identical data; average of identical = same
        result = meta_lds(ls)
        expected = ls.loci[0].ld.r
        np.testing.assert_allclose(result.r, expected, atol=1e-4)

    def test_af2_branch(self, locus_set_with_af2):
        """When LD maps have AF2, the result LD map should also have AF2."""
        result = meta_lds(locus_set_with_af2)
        assert "AF2" in result.map.columns

    def test_af2_missing_cohort_not_diluted(self):
        """AF2 for SNPs only in one cohort should equal that cohort's AF2."""
        snps_a = ["1-1000-A-G", "1-1100-A-G"]
        snps_b = ["1-1100-A-G", "1-1200-A-G"]

        def _make_ld_af2(snpids, af2_vals):
            n = len(snpids)
            bps = [int(s.split("-")[1]) for s in snpids]
            ld_map = pd.DataFrame(
                {
                    ColName.SNPID: snpids,
                    ColName.CHR: np.int8(1),
                    ColName.BP: np.array(bps, dtype=np.int32),
                    ColName.A1: "A",
                    ColName.A2: "G",
                    "AF2": np.array(af2_vals, dtype=np.float32),
                }
            )
            return LDMatrix(ld_map, np.eye(n, dtype=np.float32))

        def _make_ss(snpids):
            bps = [int(s.split("-")[1]) for s in snpids]
            return pd.DataFrame(
                {
                    ColName.SNPID: snpids,
                    ColName.CHR: np.int8(1),
                    ColName.BP: np.array(bps, dtype=np.int32),
                    ColName.EA: "A",
                    ColName.NEA: "G",
                    ColName.EAF: np.float32(0.3),
                    ColName.MAF: np.float32(0.3),
                    ColName.A1: "A",
                    ColName.A2: "G",
                    ColName.BETA: np.float32(0.1),
                    ColName.SE: np.float32(0.01),
                    ColName.P: 1e-5,
                }
            )

        locus_a = Locus(
            "EUR", "UKB", 10000, _make_ss(snps_a), 900, 3500,
            ld=_make_ld_af2(snps_a, [0.4, 0.5]),
        )
        locus_b = Locus(
            "AFR", "MVP", 5000, _make_ss(snps_b), 900, 3500,
            ld=_make_ld_af2(snps_b, [0.45, 0.7]),
        )
        ls = LocusSet([locus_a, locus_b])
        result = meta_lds(ls)
        af2 = result.map.set_index(ColName.SNPID)["AF2"]

        # SNP only in cohort A => AF2 should be 0.4, not diluted
        np.testing.assert_allclose(af2["1-1000-A-G"], 0.4, atol=1e-4)
        # SNP only in cohort B => AF2 should be 0.7, not diluted
        np.testing.assert_allclose(af2["1-1200-A-G"], 0.7, atol=1e-4)
        # Shared SNP => weighted average
        expected = 0.5 * (10000 / 15000) + 0.45 * (5000 / 15000)
        np.testing.assert_allclose(af2["1-1100-A-G"], expected, atol=1e-4)

    def test_output_dtype_float32(self, two_population_locus_set):
        """Output LD matrix should be float32."""
        result = meta_lds(two_population_locus_set)
        assert result.r.dtype == np.float32


# ===========================================================================
# TestMetaAll
# ===========================================================================


class TestMetaAll:
    """Tests for meta_all (lines 209-267)."""

    def test_returns_locus(self, two_population_locus_set):
        """Should return a Locus."""
        result = meta_all(two_population_locus_set)
        assert isinstance(result, Locus)

    def test_popu_sorted_joined(self, two_population_locus_set):
        """Popu should be sorted and joined with +."""
        result = meta_all(two_population_locus_set)
        assert result.popu == "AFR+EUR"

    def test_cohort_sorted_joined(self, two_population_locus_set):
        """Cohort should be sorted and joined with +."""
        result = meta_all(two_population_locus_set)
        assert result.cohort == "MVP+UKB"

    def test_sample_size_sum(self, two_population_locus_set):
        """Sample size should be summed."""
        result = meta_all(two_population_locus_set)
        expected = sum(loc.sample_size for loc in two_population_locus_set.loci)
        assert result.sample_size == expected

    def test_valueerror_start_mismatch(self, mismatched_start_locus_set):
        """Should raise ValueError when start positions differ."""
        with pytest.raises(ValueError, match="same start position"):
            meta_all(mismatched_start_locus_set)

    def test_valueerror_end_mismatch(self, mismatched_end_locus_set):
        """Should raise ValueError when end positions differ."""
        with pytest.raises(ValueError, match="same end position"):
            meta_all(mismatched_end_locus_set)

    def test_result_is_matched(self, two_population_locus_set):
        """Result should be matched (sumstats and LD aligned)."""
        result = meta_all(two_population_locus_set)
        assert result.is_matched


# ===========================================================================
# TestMetaByPopulation
# ===========================================================================


class TestMetaByPopulation:
    """Tests for meta_by_population (lines 270-310)."""

    def test_groups_by_population(self, mixed_population_locus_set):
        """Should return one entry per population."""
        result = meta_by_population(mixed_population_locus_set)
        assert set(result.keys()) == {"EUR", "AFR"}

    def test_multi_cohort_does_meta(self, mixed_population_locus_set):
        """Multi-cohort population should do meta-analysis (cohort has +)."""
        result = meta_by_population(mixed_population_locus_set)
        eur = result["EUR"]
        assert "+" in eur.cohort  # GWAS2+UKB

    def test_single_cohort_does_intersect(self, mixed_population_locus_set):
        """Single-cohort population should just intersect."""
        result = meta_by_population(mixed_population_locus_set)
        afr = result["AFR"]
        assert afr.cohort == "MVP"

    def test_all_results_are_locus(self, mixed_population_locus_set):
        """All results should be Locus objects."""
        result = meta_by_population(mixed_population_locus_set)
        for popu, locus in result.items():
            assert isinstance(locus, Locus), f"{popu} result is not Locus"

    def test_single_population_only_intersect(self, same_population_locus_set):
        """When all loci are same population, meta-analysis is performed."""
        result = meta_by_population(same_population_locus_set)
        assert "EUR" in result
        assert len(result) == 1
        # Multiple cohorts in same popu => meta was done
        assert "+" in result["EUR"].cohort


# ===========================================================================
# TestMeta
# ===========================================================================


class TestMeta:
    """Tests for meta() function (lines 313-357)."""

    def test_meta_all_returns_one(self, two_population_locus_set):
        """meta_all should return a LocusSet with 1 locus."""
        result = meta(two_population_locus_set, "meta_all")
        assert isinstance(result, LocusSet)
        assert result.n_loci == 1

    def test_meta_by_population_returns_per_popu(self, mixed_population_locus_set):
        """meta_by_population should return per-population loci."""
        result = meta(mixed_population_locus_set, "meta_by_population")
        assert isinstance(result, LocusSet)
        assert result.n_loci == 2  # EUR and AFR

    def test_no_meta_returns_all(self, two_population_locus_set):
        """no_meta should return all original loci (intersected)."""
        result = meta(two_population_locus_set, "no_meta")
        assert isinstance(result, LocusSet)
        assert result.n_loci == 2

    def test_invalid_method_raises(self, two_population_locus_set):
        """Invalid method should raise ValueError."""
        with pytest.raises(ValueError, match="Unsupported meta-analysis method"):
            meta(two_population_locus_set, "invalid_method")

    def test_default_method_is_meta_all(self, two_population_locus_set):
        """Default method should be meta_all."""
        result = meta(two_population_locus_set)
        assert result.n_loci == 1


# ===========================================================================
# TestMetaLocus
# ===========================================================================


class TestMetaLocus:
    """Tests for meta_locus (lines 489-566)."""

    def _make_locus_info_df(self):
        """Create a minimal locus_info DataFrame for mocking."""
        return pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB", "/fake/AFR_MVP"],
                "popu": ["EUR", "AFR"],
                "cohort": ["UKB", "MVP"],
                "sample_size": [10000, 8000],
                "chr": [1, 1],
                "start": [1000, 1000],
                "end": [3000, 3000],
                "locus_id": ["chr1_1000_3000", "chr1_1000_3000"],
            }
        )

    @patch("credtools.meta.load_locus_set")
    def test_returns_tuple(self, mock_load, tmp_path, two_population_locus_set):
        """Should return (results, het_summary) tuple."""
        mock_load.return_value = two_population_locus_set
        locus_info = self._make_locus_info_df()
        args = ("chr1_1000_3000", locus_info, str(tmp_path), "meta_all", False)
        results, het_summary = meta_locus(args)
        assert isinstance(results, list)
        assert isinstance(het_summary, pd.DataFrame)

    @patch("credtools.meta.load_locus_set")
    def test_creates_output_files(self, mock_load, tmp_path, two_population_locus_set):
        """Should create sumstats.gz, ld.npz, ldmap.gz, and heterogeneity files."""
        mock_load.return_value = two_population_locus_set
        locus_info = self._make_locus_info_df()
        args = ("chr1_1000_3000", locus_info, str(tmp_path), "meta_all", False)
        results, _ = meta_locus(args)

        out_dir = tmp_path / "chr1_1000_3000"
        assert out_dir.exists()
        # Check heterogeneity files exist
        assert (out_dir / "heterogeneity.txt.gz").exists()
        # Check at least one result file set
        assert len(results) >= 1
        prefix = results[0][6]  # out_prefix is 7th element
        assert os.path.exists(f"{prefix}.sumstats.gz")
        assert os.path.exists(f"{prefix}.ld.npz")
        assert os.path.exists(f"{prefix}.ldmap.gz")

    @patch("credtools.meta.load_locus_set")
    def test_result_metadata(self, mock_load, tmp_path, two_population_locus_set):
        """Result records should have correct chrom, popu, sample_size."""
        mock_load.return_value = two_population_locus_set
        locus_info = self._make_locus_info_df()
        args = ("chr1_1000_3000", locus_info, str(tmp_path), "meta_all", False)
        results, _ = meta_locus(args)

        # meta_all should produce 1 result
        assert len(results) == 1
        res = results[0]
        assert res[0] == 1  # chrom
        assert res[4] == 18000  # sample_size = 10000 + 8000

    @patch("credtools.meta.load_locus_set")
    def test_no_meta_creates_multiple(
        self, mock_load, tmp_path, two_population_locus_set
    ):
        """no_meta should create multiple result entries."""
        mock_load.return_value = two_population_locus_set
        locus_info = self._make_locus_info_df()
        args = ("chr1_1000_3000", locus_info, str(tmp_path), "no_meta", False)
        results, _ = meta_locus(args)

        assert len(results) == 2


# ===========================================================================
# TestMetaLoci
# ===========================================================================


class TestMetaLoci:
    """Tests for meta_loci (lines 569-665)."""

    @patch("credtools.meta.Pool")
    @patch("credtools.meta.check_loci_info")
    def test_creates_loci_info(self, mock_check, mock_pool_cls, tmp_path):
        """Should create loci_info.txt in outdir."""
        # Build a minimal loci_info
        loci_info = pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "locus_id": ["chr1_1000_3000"],
            }
        )
        mock_check.return_value = loci_info

        # Mock Pool to return fake results
        fake_result = (
            [[1, 1000, 3000, "EUR", 10000, "UKB", "/fake/prefix", "chr1_1000_3000"]],
            pd.DataFrame(
                {"popu": ["EUR"], "cohort": ["UKB"], "locus_id": ["chr1_1000_3000"]}
            ),
        )
        mock_pool = MagicMock()
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool.imap_unordered.return_value = iter([fake_result])
        mock_pool_cls.return_value = mock_pool

        input_path = tmp_path / "input.txt"
        loci_info.to_csv(str(input_path), sep="\t", index=False)

        outdir = str(tmp_path / "output")
        meta_loci(str(input_path), outdir, threads=1)

        assert os.path.exists(f"{outdir}/loci_info.txt")

    @patch("credtools.meta.Pool")
    @patch("credtools.meta.check_loci_info")
    def test_creates_global_heterogeneity(self, mock_check, mock_pool_cls, tmp_path):
        """Should create heterogeneity.txt.gz when het_summary is present."""
        loci_info = pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "locus_id": ["chr1_1000_3000"],
            }
        )
        mock_check.return_value = loci_info

        het_summary = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "ld_4th_moment_mean": [0.5],
                "locus_id": ["chr1_1000_3000"],
            }
        )
        fake_result = (
            [[1, 1000, 3000, "EUR", 10000, "UKB", "/fake/prefix", "chr1_1000_3000"]],
            het_summary,
        )
        mock_pool = MagicMock()
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool.imap_unordered.return_value = iter([fake_result])
        mock_pool_cls.return_value = mock_pool

        input_path = tmp_path / "input.txt"
        loci_info.to_csv(str(input_path), sep="\t", index=False)

        outdir = str(tmp_path / "output")
        meta_loci(str(input_path), outdir, threads=1)

        assert os.path.exists(f"{outdir}/heterogeneity.txt.gz")
        loaded = pd.read_csv(
            f"{outdir}/heterogeneity.txt.gz", sep="\t", compression="gzip"
        )
        assert "locus_id" in loaded.columns

    @patch("credtools.meta.Pool")
    @patch("credtools.meta.check_loci_info")
    def test_no_het_summary_no_global_file(self, mock_check, mock_pool_cls, tmp_path):
        """Should not create heterogeneity.txt.gz when het_summary is empty."""
        loci_info = pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "locus_id": ["chr1_1000_3000"],
            }
        )
        mock_check.return_value = loci_info

        # Empty het_summary
        fake_result = (
            [[1, 1000, 3000, "EUR", 10000, "UKB", "/fake/prefix", "chr1_1000_3000"]],
            pd.DataFrame(),  # empty
        )
        mock_pool = MagicMock()
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool.imap_unordered.return_value = iter([fake_result])
        mock_pool_cls.return_value = mock_pool

        input_path = tmp_path / "input.txt"
        loci_info.to_csv(str(input_path), sep="\t", index=False)

        outdir = str(tmp_path / "output")
        meta_loci(str(input_path), outdir, threads=1)

        assert not os.path.exists(f"{outdir}/heterogeneity.txt.gz")


# ===========================================================================
# TestHeterogeneitySummaryEdgeCases
# ===========================================================================


class TestHeterogeneitySummaryEdgeCases:
    """Edge case tests for heterogeneity_summary (lines 416/425-427)."""

    def test_ld_decay_no_matching_cohort(self):
        """When ld_decay has no matching cohort data, decay_rate should be NaN."""
        locus = _make_locus("EUR", "UKB", seed=42)
        locus_set = LocusSet([locus])

        # ld_decay with a different cohort label
        ld_decay_df = pd.DataFrame(
            {
                "cohort": ["UNKNOWN_COHORT"],
                "distance_kb": [1.0],
                "r2_avg": [0.5],
                "decay_rate": [0.1],
            }
        )
        het_metrics = {
            "ld_4th_moment": pd.DataFrame({"EUR_UKB": [0.5, 0.6]}),
            "ld_decay": ld_decay_df,
        }
        summary = heterogeneity_summary(het_metrics, locus_set)
        assert pd.isna(summary.iloc[0]["ld_decay_rate"])

    def test_ld_decay_none(self):
        """When ld_decay is None, decay_rate should be NaN."""
        locus = _make_locus("EUR", "UKB", seed=42)
        locus_set = LocusSet([locus])

        het_metrics = {
            "ld_4th_moment": pd.DataFrame({"EUR_UKB": [0.5, 0.6]}),
            "ld_decay": None,
        }
        summary = heterogeneity_summary(het_metrics, locus_set)
        assert pd.isna(summary.iloc[0]["ld_decay_rate"])

    def test_ld_4th_moment_no_matching_column(self):
        """When ld_4th_moment has no matching column, value should be NaN."""
        locus = _make_locus("EUR", "UKB", seed=42)
        locus_set = LocusSet([locus])

        het_metrics = {
            "ld_4th_moment": pd.DataFrame({"OTHER_COL": [0.5]}),
            "ld_decay": pd.DataFrame(
                {
                    "cohort": ["EUR_UKB"],
                    "distance_kb": [1.0],
                    "r2_avg": [0.5],
                    "decay_rate": [0.1],
                }
            ),
        }
        summary = heterogeneity_summary(het_metrics, locus_set)
        assert pd.isna(summary.iloc[0]["ld_4th_moment_mean"])


# ===========================================================================
# TestRecoverCompletedLocus
# ===========================================================================


class TestRecoverCompletedLocus:
    """Tests for recover_completed_locus function."""

    def _make_prev_loci_info(self):
        """Create a prev_loci_info DataFrame with one locus having two prefixes."""
        return pd.DataFrame(
            {
                "chr": [1, 1],
                "start": [1000, 1000],
                "end": [3000, 3000],
                "popu": ["EUR", "AFR"],
                "sample_size": [10000, 8000],
                "cohort": ["UKB", "MVP"],
                "prefix": [
                    "/out/chr1_1000_3000/EUR_UKB",
                    "/out/chr1_1000_3000/AFR_MVP",
                ],
                "locus_id": ["chr1_1000_3000", "chr1_1000_3000"],
            }
        )

    def test_returns_none_when_no_prev_info(self):
        """prev_loci_info=None should return None."""
        result = recover_completed_locus("chr1_1000_3000", "/fake", None)
        assert result is None

    def test_returns_none_when_dir_missing(self, tmp_path):
        """Locus directory does not exist should return None."""
        prev = self._make_prev_loci_info()
        result = recover_completed_locus(
            "chr1_1000_3000", str(tmp_path / "nonexistent"), prev
        )
        assert result is None

    def test_returns_none_when_files_incomplete(self, tmp_path):
        """Missing one of the 3 required files should return None."""
        outdir = str(tmp_path)
        locus_dir = tmp_path / "chr1_1000_3000"
        locus_dir.mkdir()

        # Update prefixes to point to tmp_path
        prev = pd.DataFrame(
            {
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "popu": ["EUR"],
                "sample_size": [10000],
                "cohort": ["UKB"],
                "prefix": [str(locus_dir / "EUR_UKB")],
                "locus_id": ["chr1_1000_3000"],
            }
        )

        # Create only sumstats.gz, missing ld.npz and ldmap.gz
        (locus_dir / "EUR_UKB.sumstats.gz").write_bytes(b"fake")

        result = recover_completed_locus("chr1_1000_3000", outdir, prev)
        assert result is None

    def test_returns_none_when_locus_id_not_in_prev(self, tmp_path):
        """Locus id not found in prev should return None."""
        prev = self._make_prev_loci_info()
        locus_dir = tmp_path / "chr2_5000_6000"
        locus_dir.mkdir()

        result = recover_completed_locus("chr2_5000_6000", str(tmp_path), prev)
        assert result is None

    def test_returns_results_when_complete(self, tmp_path):
        """All 3 files present should return (results, het_summary)."""
        outdir = str(tmp_path)
        locus_dir = tmp_path / "chr1_1000_3000"
        locus_dir.mkdir()

        prefix = str(locus_dir / "EUR_UKB")
        prev = pd.DataFrame(
            {
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "popu": ["EUR"],
                "sample_size": [10000],
                "cohort": ["UKB"],
                "prefix": [prefix],
                "locus_id": ["chr1_1000_3000"],
            }
        )

        # Create all 3 required files
        for ext in [".sumstats.gz", ".ld.npz", ".ldmap.gz"]:
            (locus_dir / f"EUR_UKB{ext}").write_bytes(b"fake")

        # Create heterogeneity file
        het_df = pd.DataFrame(
            {"popu": ["EUR"], "cohort": ["UKB"], "locus_id": ["chr1_1000_3000"]}
        )
        het_df.to_csv(
            str(locus_dir / "heterogeneity.txt.gz"),
            sep="\t",
            index=False,
            compression="gzip",
        )

        result = recover_completed_locus("chr1_1000_3000", outdir, prev)
        assert result is not None
        results, het_summary = result
        assert len(results) == 1
        assert results[0][0] == 1  # chr
        assert results[0][3] == "EUR"  # popu
        assert results[0][4] == 10000  # sample_size
        assert results[0][7] == "chr1_1000_3000"  # locus_id
        assert not het_summary.empty

    def test_returns_empty_het_when_no_het_file(self, tmp_path):
        """No heterogeneity.txt.gz should return empty DataFrame."""
        outdir = str(tmp_path)
        locus_dir = tmp_path / "chr1_1000_3000"
        locus_dir.mkdir()

        prefix = str(locus_dir / "EUR_UKB")
        prev = pd.DataFrame(
            {
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "popu": ["EUR"],
                "sample_size": [10000],
                "cohort": ["UKB"],
                "prefix": [prefix],
                "locus_id": ["chr1_1000_3000"],
            }
        )

        for ext in [".sumstats.gz", ".ld.npz", ".ldmap.gz"]:
            (locus_dir / f"EUR_UKB{ext}").write_bytes(b"fake")

        result = recover_completed_locus("chr1_1000_3000", outdir, prev)
        assert result is not None
        _, het_summary = result
        assert het_summary.empty


# ===========================================================================
# TestMetaLoci – skip tests
# ===========================================================================


class TestMetaLociSkip:
    """Tests for meta_loci skip parameter."""

    def _setup_fake_pool(self, mock_pool_cls, fake_results):
        """Configure mock Pool to return fake_results."""
        mock_pool = MagicMock()
        mock_pool.__enter__ = MagicMock(return_value=mock_pool)
        mock_pool.__exit__ = MagicMock(return_value=False)
        mock_pool.imap_unordered.return_value = iter(fake_results)
        mock_pool_cls.return_value = mock_pool
        return mock_pool

    @patch("credtools.meta.meta_locus")
    @patch("credtools.meta.recover_completed_locus")
    @patch("credtools.meta.Pool")
    @patch("credtools.meta.check_loci_info")
    def test_skip_completed_locus(
        self, mock_check, mock_pool_cls, mock_recover, mock_meta_locus, tmp_path
    ):
        """skip=True: completed locus should not go through Pool."""
        loci_info = pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB", "/fake/AFR_MVP"],
                "popu": ["EUR", "AFR"],
                "cohort": ["UKB", "MVP"],
                "sample_size": [10000, 8000],
                "chr": [1, 2],
                "start": [1000, 5000],
                "end": [3000, 7000],
                "locus_id": ["chr1_1000_3000", "chr2_5000_7000"],
            }
        )
        mock_check.return_value = loci_info

        # First locus is recovered, second is not
        recovered_results = [
            [1, 1000, 3000, "EUR", 10000, "UKB", "/fake/prefix1", "chr1_1000_3000"]
        ]
        recovered_het = pd.DataFrame(
            {"popu": ["EUR"], "cohort": ["UKB"], "locus_id": ["chr1_1000_3000"]}
        )
        mock_recover.side_effect = [
            (recovered_results, recovered_het),  # chr1 recovered
            None,  # chr2 not recovered
        ]

        # Pool processes the second locus
        pool_result = (
            [[2, 5000, 7000, "AFR", 8000, "MVP", "/fake/prefix2", "chr2_5000_7000"]],
            pd.DataFrame(
                {"popu": ["AFR"], "cohort": ["MVP"], "locus_id": ["chr2_5000_7000"]}
            ),
        )
        mock_pool = self._setup_fake_pool(mock_pool_cls, [pool_result])

        # Write prev loci_info.txt so skip can read it
        outdir = str(tmp_path / "output")
        os.makedirs(outdir, exist_ok=True)
        loci_info.to_csv(f"{outdir}/loci_info.txt", sep="\t", index=False)

        input_path = tmp_path / "input.txt"
        loci_info.to_csv(str(input_path), sep="\t", index=False)

        meta_loci(str(input_path), outdir, threads=1, skip=True)

        # Verify loci_info.txt has both loci
        result_df = pd.read_csv(f"{outdir}/loci_info.txt", sep="\t")
        assert len(result_df) == 2
        assert set(result_df["locus_id"]) == {"chr1_1000_3000", "chr2_5000_7000"}

    @patch("credtools.meta.recover_completed_locus")
    @patch("credtools.meta.Pool")
    @patch("credtools.meta.check_loci_info")
    def test_skip_false_processes_all(
        self, mock_check, mock_pool_cls, mock_recover, tmp_path
    ):
        """skip=False (default): recover should not be called."""
        loci_info = pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "locus_id": ["chr1_1000_3000"],
            }
        )
        mock_check.return_value = loci_info

        fake_result = (
            [[1, 1000, 3000, "EUR", 10000, "UKB", "/fake/prefix", "chr1_1000_3000"]],
            pd.DataFrame(
                {"popu": ["EUR"], "cohort": ["UKB"], "locus_id": ["chr1_1000_3000"]}
            ),
        )
        self._setup_fake_pool(mock_pool_cls, [fake_result])

        input_path = tmp_path / "input.txt"
        loci_info.to_csv(str(input_path), sep="\t", index=False)

        outdir = str(tmp_path / "output")
        meta_loci(str(input_path), outdir, threads=1, skip=False)

        mock_recover.assert_not_called()

    @patch("credtools.meta.recover_completed_locus")
    @patch("credtools.meta.Pool")
    @patch("credtools.meta.check_loci_info")
    def test_skip_no_prev_loci_info(
        self, mock_check, mock_pool_cls, mock_recover, tmp_path
    ):
        """skip=True but no prev loci_info.txt: all loci processed via Pool."""
        loci_info = pd.DataFrame(
            {
                "prefix": ["/fake/EUR_UKB"],
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "sample_size": [10000],
                "chr": [1],
                "start": [1000],
                "end": [3000],
                "locus_id": ["chr1_1000_3000"],
            }
        )
        mock_check.return_value = loci_info

        fake_result = (
            [[1, 1000, 3000, "EUR", 10000, "UKB", "/fake/prefix", "chr1_1000_3000"]],
            pd.DataFrame(
                {"popu": ["EUR"], "cohort": ["UKB"], "locus_id": ["chr1_1000_3000"]}
            ),
        )
        self._setup_fake_pool(mock_pool_cls, [fake_result])

        input_path = tmp_path / "input.txt"
        loci_info.to_csv(str(input_path), sep="\t", index=False)

        outdir = str(tmp_path / "output_fresh")
        # No prev loci_info.txt exists
        meta_loci(str(input_path), outdir, threads=1, skip=True)

        # recover should not be called because prev_loci_info is None
        mock_recover.assert_not_called()
        assert os.path.exists(f"{outdir}/loci_info.txt")
