"""Tests for credtools.locus module."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet, check_loci_info, intersect_sumstat_ld


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_simple_locus(
    popu="EUR",
    cohort="test",
    n_snps=10,
    chrom=1,
    start=1000,
    end=2000,
    sample_size=10000,
):
    """Create a simple Locus for testing."""
    bps = [start + i * 100 for i in range(n_snps)]
    snpids = [f"{chrom}-{bp}-A-G" for bp in bps]

    sumstats = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: [chrom] * n_snps,
            ColName.BP: bps,
            ColName.RSID: snpids,
            ColName.EA: ["A"] * n_snps,
            ColName.NEA: ["G"] * n_snps,
            ColName.EAF: [0.3] * n_snps,
            ColName.BETA: [0.1] * n_snps,
            ColName.SE: [0.01] * n_snps,
            ColName.P: [1e-8] * n_snps,
        }
    )

    r = np.eye(n_snps)
    ld_map = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: [chrom] * n_snps,
            ColName.BP: bps,
            ColName.A1: ["A"] * n_snps,
            ColName.A2: ["G"] * n_snps,
        }
    )
    ld = LDMatrix(ld_map, r)
    return Locus(popu, cohort, sample_size, sumstats, start, end, ld=ld)


# ---------------------------------------------------------------------------
# TestLocusPrefix
# ---------------------------------------------------------------------------
class TestLocusPrefix:
    """Tests for Locus.prefix property."""

    def test_normal_prefix(self):
        locus = _make_simple_locus(popu="EUR", cohort="UKB")
        assert locus.prefix == "EUR_UKB"

    def test_meta_prefix_with_plus(self):
        locus = _make_simple_locus(cohort="UKB+BBJ+APCDR")
        prefix = locus.prefix
        assert prefix.startswith("EUR_meta3cohorts_")
        assert len(prefix.split("_")) == 3

    def test_meta_prefix_hash_deterministic(self):
        locus1 = _make_simple_locus(cohort="UKB+BBJ")
        locus2 = _make_simple_locus(cohort="UKB+BBJ")
        assert locus1.prefix == locus2.prefix


# ---------------------------------------------------------------------------
# TestLocusSetProperties
# ---------------------------------------------------------------------------
class TestLocusSetProperties:
    """Tests for LocusSet property validation."""

    def test_different_chrom_raises(self):
        l1 = _make_simple_locus(chrom=1)
        l2 = _make_simple_locus(chrom=2)
        ls = LocusSet([l1, l2])
        with pytest.raises(ValueError, match="chromosomes"):
            ls.chrom

    def test_different_start_raises(self):
        l1 = _make_simple_locus(start=1000, end=2000)
        l2 = _make_simple_locus(start=1500, end=2000)
        ls = LocusSet([l1, l2])
        with pytest.raises(ValueError, match="start position"):
            ls.start

    def test_different_end_raises(self):
        l1 = _make_simple_locus(start=1000, end=2000)
        l2 = _make_simple_locus(start=1000, end=2500)
        ls = LocusSet([l1, l2])
        with pytest.raises(ValueError, match="end position"):
            ls.end

    def test_same_properties(self):
        l1 = _make_simple_locus(popu="EUR", cohort="C1")
        l2 = _make_simple_locus(popu="EAS", cohort="C2")
        ls = LocusSet([l1, l2])
        assert ls.chrom == 1
        assert ls.start == 1000
        assert ls.end == 2000
        assert ls.n_loci == 2

    def test_locus_id(self):
        l1 = _make_simple_locus()
        ls = LocusSet([l1])
        assert ls.locus_id == "1:1000-2000"


# ---------------------------------------------------------------------------
# TestCheckLociInfo
# ---------------------------------------------------------------------------
class TestCheckLociInfo:
    """Tests for check_loci_info function."""

    def _make_valid_df(self, n_rows=1):
        return pd.DataFrame(
            {
                "prefix": [f"EUR_test{i}" for i in range(n_rows)],
                "popu": ["EUR"] * n_rows,
                "cohort": [f"test{i}" for i in range(n_rows)],
                "sample_size": [10000] * n_rows,
                "chr": [1] * n_rows,
                "start": [1000] * n_rows,
                "end": [2000] * n_rows,
                "locus_id": ["locus_1"] * n_rows,
            }
        )

    def test_valid_df_passes(self):
        df = self._make_valid_df()
        result = check_loci_info(df)
        assert len(result) == 1

    def test_missing_column_raises(self):
        df = self._make_valid_df()
        df = df.drop(columns=["chr"])
        with pytest.raises(ValueError, match="Missing required columns"):
            check_loci_info(df)

    def test_negative_sample_size_raises(self):
        df = self._make_valid_df()
        df["sample_size"] = -100
        with pytest.raises(ValueError, match="Sample size must be positive"):
            check_loci_info(df)

    def test_chr_out_of_range_raises(self):
        df = self._make_valid_df()
        df["chr"] = 30
        with pytest.raises(ValueError, match="Chromosome must be between"):
            check_loci_info(df)

    def test_end_before_start_raises(self):
        df = self._make_valid_df()
        df["end"] = 500
        with pytest.raises(ValueError, match="End position must be greater"):
            check_loci_info(df)

    def test_duplicate_raises(self):
        df = self._make_valid_df(n_rows=2)
        df["cohort"] = "same"
        with pytest.raises(ValueError, match="unique"):
            check_loci_info(df)

    def test_start_zero_raises(self):
        df = self._make_valid_df()
        df["start"] = 0
        with pytest.raises(ValueError, match="Start position must be positive"):
            check_loci_info(df)


# ---------------------------------------------------------------------------
# TestLoadLocusErrors
# ---------------------------------------------------------------------------
class TestLoadLocusErrors:
    """Tests for load_locus error paths."""

    def test_sumstats_not_found(self, tmp_path):
        from credtools.locus import load_locus

        with pytest.raises(ValueError, match="Sumstats file not found"):
            load_locus(
                str(tmp_path / "nonexistent"),
                "EUR",
                "test",
                10000,
                1000,
                2000,
            )

    def test_ld_not_found(self, tmp_path):
        from credtools.locus import load_locus

        # Create sumstats file but no LD
        prefix = str(tmp_path / "test")
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: ["1-1000-A-G"],
                ColName.CHR: [1],
                ColName.BP: [1000],
                ColName.EA: ["A"],
                ColName.NEA: ["G"],
                ColName.EAF: [0.3],
                ColName.BETA: [0.1],
                ColName.SE: [0.01],
                ColName.P: [1e-8],
            }
        )
        sumstats.to_csv(f"{prefix}.sumstat", sep="\t", index=False)
        with pytest.raises(ValueError, match="LD matrix file not found"):
            load_locus(prefix, "EUR", "test", 10000, 1000, 2000)

    def test_ldmap_not_found(self, tmp_path):
        from credtools.locus import load_locus

        prefix = str(tmp_path / "test")
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: ["1-1000-A-G"],
                ColName.CHR: [1],
                ColName.BP: [1000],
                ColName.EA: ["A"],
                ColName.NEA: ["G"],
                ColName.EAF: [0.3],
                ColName.BETA: [0.1],
                ColName.SE: [0.01],
                ColName.P: [1e-8],
            }
        )
        sumstats.to_csv(f"{prefix}.sumstat", sep="\t", index=False)
        np.savez_compressed(f"{prefix}.ld.npz", ld=np.eye(1))
        with pytest.raises(ValueError, match="LD map file not found"):
            load_locus(prefix, "EUR", "test", 10000, 1000, 2000)


# ---------------------------------------------------------------------------
# TestIntersectSumstatLd
# ---------------------------------------------------------------------------
class TestIntersectSumstatLd:
    """Tests for intersect_sumstat_ld function."""

    def test_ld_none_raises(self):
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: ["s1"],
                ColName.CHR: [1],
                ColName.BP: [100],
                ColName.EA: ["A"],
                ColName.NEA: ["G"],
                ColName.EAF: [0.3],
                ColName.BETA: [0.1],
                ColName.SE: [0.01],
                ColName.P: [1e-8],
            }
        )
        locus = Locus("EUR", "test", 10000, sumstats, 100, 200)
        # Force ld to None to trigger the ValueError path
        locus.ld = None
        with pytest.raises(ValueError, match="LD matrix not found"):
            intersect_sumstat_ld(locus)

    def test_no_common_variants_raises(self):
        n = 3
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: [f"sum_{i}" for i in range(n)],
                ColName.CHR: [1] * n,
                ColName.BP: [100, 200, 300],
                ColName.EA: ["A"] * n,
                ColName.NEA: ["G"] * n,
                ColName.EAF: [0.3] * n,
                ColName.BETA: [0.1] * n,
                ColName.SE: [0.01] * n,
                ColName.P: [1e-8] * n,
            }
        )
        ld_map = pd.DataFrame(
            {
                ColName.SNPID: [f"ld_{i}" for i in range(n)],
                ColName.CHR: [1] * n,
                ColName.BP: [100, 200, 300],
                ColName.A1: ["A"] * n,
                ColName.A2: ["G"] * n,
            }
        )
        ld = LDMatrix(ld_map, np.eye(n))
        locus = Locus("EUR", "test", 10000, sumstats, 100, 400, ld=ld)
        with pytest.raises(ValueError, match="No common Variant IDs"):
            intersect_sumstat_ld(locus)

    def test_matched_returns_same(self):
        locus = _make_simple_locus()
        result = intersect_sumstat_ld(locus)
        assert result.sumstats[ColName.SNPID].equals(locus.sumstats[ColName.SNPID])

    def test_few_common_variants_warns(self):
        """When only <=10 common variants, a warning should be logged."""
        n = 5
        snpids = [f"1-{1000 + i * 100}-A-G" for i in range(n)]
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n,
                ColName.BP: [1000 + i * 100 for i in range(n)],
                ColName.EA: ["A"] * n,
                ColName.NEA: ["G"] * n,
                ColName.EAF: [0.3] * n,
                ColName.BETA: [0.1] * n,
                ColName.SE: [0.01] * n,
                ColName.P: [1e-8] * n,
            }
        )
        # LD map only has first 3 matching + 2 extra
        ld_snpids = snpids[:3] + ["extra_1", "extra_2"]
        ld_map = pd.DataFrame(
            {
                ColName.SNPID: ld_snpids,
                ColName.CHR: [1] * 5,
                ColName.BP: [1000, 1100, 1200, 1300, 1400],
                ColName.A1: ["A"] * 5,
                ColName.A2: ["G"] * 5,
            }
        )
        ld = LDMatrix(ld_map, np.eye(5))
        locus = Locus("EUR", "test", 10000, sumstats, 1000, 2000, ld=ld)
        result = intersect_sumstat_ld(locus)
        assert result.n_snps == 3


# ---------------------------------------------------------------------------
# TestLocusCopy
# ---------------------------------------------------------------------------
class TestLocusCopy:
    """Tests for Locus.copy method."""

    def test_copy_preserves_properties(self):
        locus = _make_simple_locus()
        copied = locus.copy()
        assert copied.popu == locus.popu
        assert copied.cohort == locus.cohort
        assert copied.sample_size == locus.sample_size
        assert copied.n_snps == locus.n_snps

    def test_copy_is_independent(self):
        locus = _make_simple_locus()
        copied = locus.copy()
        copied.sumstats.iloc[0, 0] = "modified"
        assert locus.sumstats.iloc[0, 0] != "modified"

    def test_copy_preserves_lambda_s(self):
        locus = _make_simple_locus()
        locus.lambda_s = 0.95
        copied = locus.copy()
        assert copied.lambda_s == 0.95


# ---------------------------------------------------------------------------
# TestLocusProperties
# ---------------------------------------------------------------------------
class TestLocusProperties:
    """Tests for basic Locus properties."""

    def test_locus_id(self):
        locus = _make_simple_locus(popu="EUR", cohort="UKB")
        assert locus.locus_id == "EUR_UKB_chr1:1000-2000"

    def test_repr(self):
        locus = _make_simple_locus()
        r = repr(locus)
        assert "Locus(" in r
        assert "EUR" in r

    def test_is_matched_true(self):
        locus = _make_simple_locus()
        assert locus.is_matched

    def test_no_ld_creates_empty(self):
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: ["s1"],
                ColName.CHR: [1],
                ColName.BP: [100],
                ColName.EA: ["A"],
                ColName.NEA: ["G"],
                ColName.EAF: [0.3],
                ColName.BETA: [0.1],
                ColName.SE: [0.01],
                ColName.P: [1e-8],
            }
        )
        locus = Locus("EUR", "test", 10000, sumstats, 100, 200)
        assert locus.ld.r.shape == (0,)
