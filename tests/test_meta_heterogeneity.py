"""Tests for heterogeneity computation in meta module."""

import os

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet


def _make_locus(popu: str, cohort: str, seed: int = 42) -> Locus:
    """Create a test locus with realistic LD structure."""
    rng = np.random.default_rng(seed)
    n_snps = 20
    bps = np.arange(1000, 1000 + n_snps * 100, 100)
    snpids = [f"1-{bp}-A-G" for bp in bps]

    sumstats = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: [1] * n_snps,
            ColName.BP: bps,
            ColName.EA: ["A"] * n_snps,
            ColName.NEA: ["G"] * n_snps,
            ColName.EAF: rng.uniform(0.1, 0.5, n_snps),
            ColName.MAF: rng.uniform(0.1, 0.5, n_snps),
            ColName.A1: ["A"] * n_snps,
            ColName.A2: ["G"] * n_snps,
            ColName.BETA: rng.normal(0, 0.1, n_snps),
            ColName.SE: rng.uniform(0.01, 0.05, n_snps),
            ColName.P: rng.uniform(1e-10, 0.05, n_snps),
        }
    )

    # Create a positive-definite LD matrix via random correlation
    A = rng.normal(size=(n_snps, n_snps))
    r = A @ A.T
    d = np.sqrt(np.diag(r))
    r = r / np.outer(d, d)
    r = r.astype(np.float32)

    ld = LDMatrix(sumstats.copy(), r)

    return Locus(
        popu=popu,
        cohort=cohort,
        sample_size=10000,
        sumstats=sumstats,
        locus_start=1000,
        locus_end=1000 + n_snps * 100,
        ld=ld,
    )


@pytest.fixture
def multi_cohort_locus_set():
    """Create a LocusSet with multiple cohorts."""
    locus1 = _make_locus("EUR", "UKB", seed=42)
    locus2 = _make_locus("AFR", "MVP", seed=123)
    locus3 = _make_locus("EAS", "BBJ", seed=456)
    return LocusSet([locus1, locus2, locus3])


@pytest.fixture
def single_locus_set():
    """Create a LocusSet with a single locus."""
    locus1 = _make_locus("EUR", "UKB", seed=42)
    return LocusSet([locus1])


class TestComputeHeterogeneity:
    """Tests for compute_heterogeneity function."""

    def test_multi_cohort_returns_all_metrics(self, multi_cohort_locus_set):
        """Multi-cohort should return all 4 heterogeneity metrics."""
        from credtools.meta import compute_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)

        assert "ld_4th_moment" in het
        assert "ld_decay" in het
        assert "cochran_q" in het
        assert "snp_missingness" in het
        assert len(het) == 4

    def test_single_locus_only_ld_metrics(self, single_locus_set):
        """Single locus should only return ld_decay and ld_4th_moment."""
        from credtools.meta import compute_heterogeneity

        het = compute_heterogeneity(single_locus_set)

        assert "ld_4th_moment" in het
        assert "ld_decay" in het
        assert "cochran_q" not in het
        assert "snp_missingness" not in het
        assert len(het) == 2

    def test_cochran_q_has_expected_columns(self, multi_cohort_locus_set):
        """Cochran Q result should contain Q, Q_pvalue, I_squared columns."""
        from credtools.meta import compute_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)
        cq = het["cochran_q"]

        assert isinstance(cq, pd.DataFrame)
        assert "Q" in cq.columns
        assert "Q_pvalue" in cq.columns
        assert "I_squared" in cq.columns

    def test_ld_decay_has_expected_columns(self, multi_cohort_locus_set):
        """LD decay result should contain distance_kb and r2_avg columns."""
        from credtools.meta import compute_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)
        ld = het["ld_decay"]

        assert isinstance(ld, pd.DataFrame)
        assert "distance_kb" in ld.columns
        assert "r2_avg" in ld.columns
        assert "cohort" in ld.columns

    def test_all_results_are_dataframes(self, multi_cohort_locus_set):
        """All heterogeneity results should be DataFrames."""
        from credtools.meta import compute_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)

        for name, data in het.items():
            assert isinstance(data, pd.DataFrame), f"{name} is not a DataFrame"


class TestSaveHeterogeneity:
    """Tests for save_heterogeneity function."""

    def test_creates_files(self, tmp_path, multi_cohort_locus_set):
        """Save should create .txt.gz files for each metric."""
        from credtools.meta import compute_heterogeneity, save_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)
        out_dir = str(tmp_path / "het_output")
        save_heterogeneity(het, out_dir)

        assert os.path.isdir(out_dir)
        for name in het:
            filepath = os.path.join(out_dir, f"{name}.txt.gz")
            assert os.path.exists(filepath), f"{filepath} not created"

    def test_files_are_readable(self, tmp_path, multi_cohort_locus_set):
        """Saved files should be readable as gzipped TSVs."""
        from credtools.meta import compute_heterogeneity, save_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)
        out_dir = str(tmp_path / "het_output")
        save_heterogeneity(het, out_dir)

        for name in het:
            filepath = os.path.join(out_dir, f"{name}.txt.gz")
            loaded = pd.read_csv(filepath, sep="\t", compression="gzip")
            assert len(loaded) > 0
            assert list(loaded.columns) == list(het[name].columns)

    def test_creates_output_dir(self, tmp_path, single_locus_set):
        """Save should create output directory if it doesn't exist."""
        from credtools.meta import compute_heterogeneity, save_heterogeneity

        het = compute_heterogeneity(single_locus_set)
        out_dir = str(tmp_path / "nested" / "dir")
        save_heterogeneity(het, out_dir)

        assert os.path.isdir(out_dir)


class TestHeterogeneitySummary:
    """Tests for heterogeneity_summary function."""

    def test_multi_cohort_returns_all_columns(self, multi_cohort_locus_set):
        """Multi-cohort summary should contain all expected columns."""
        from credtools.meta import compute_heterogeneity, heterogeneity_summary

        het = compute_heterogeneity(multi_cohort_locus_set)
        summary = heterogeneity_summary(het, multi_cohort_locus_set)

        assert isinstance(summary, pd.DataFrame)
        expected_cols = [
            "popu",
            "cohort",
            "ld_4th_moment_mean",
            "ld_decay_rate",
            "missing_rate",
            "cochran_q_median",
            "i_squared_median",
            "n_het_snps",
        ]
        for col in expected_cols:
            assert col in summary.columns, f"Missing column: {col}"

    def test_multi_cohort_row_count_equals_cohort_count(self, multi_cohort_locus_set):
        """Number of rows should equal number of cohorts."""
        from credtools.meta import compute_heterogeneity, heterogeneity_summary

        het = compute_heterogeneity(multi_cohort_locus_set)
        summary = heterogeneity_summary(het, multi_cohort_locus_set)

        assert len(summary) == len(multi_cohort_locus_set.loci)

    def test_single_locus_has_nan_multi_cohort_metrics(self, single_locus_set):
        """Single locus should have NaN for cochran_q/i_squared/n_het_snps/missing_rate."""
        from credtools.meta import compute_heterogeneity, heterogeneity_summary

        het = compute_heterogeneity(single_locus_set)
        summary = heterogeneity_summary(het, single_locus_set)

        assert len(summary) == 1
        row = summary.iloc[0]
        assert pd.isna(row["cochran_q_median"])
        assert pd.isna(row["i_squared_median"])
        assert pd.isna(row["n_het_snps"])
        assert pd.isna(row["missing_rate"])
        # LD metrics should be present
        assert not pd.isna(row["ld_4th_moment_mean"])
        assert not pd.isna(row["ld_decay_rate"])

    def test_correct_cohorts(self, multi_cohort_locus_set):
        """Popu/cohort columns should match input loci."""
        from credtools.meta import compute_heterogeneity, heterogeneity_summary

        het = compute_heterogeneity(multi_cohort_locus_set)
        summary = heterogeneity_summary(het, multi_cohort_locus_set)

        expected_popus = [loc.popu for loc in multi_cohort_locus_set.loci]
        expected_cohorts = [loc.cohort for loc in multi_cohort_locus_set.loci]
        assert list(summary["popu"]) == expected_popus
        assert list(summary["cohort"]) == expected_cohorts

    def test_cochran_q_values_same_across_cohorts(self, multi_cohort_locus_set):
        """Cochran Q summary values should be identical for all cohorts in same locus."""
        from credtools.meta import compute_heterogeneity, heterogeneity_summary

        het = compute_heterogeneity(multi_cohort_locus_set)
        summary = heterogeneity_summary(het, multi_cohort_locus_set)

        # All rows should have the same cochran_q_median value
        assert summary["cochran_q_median"].nunique() == 1
        assert summary["i_squared_median"].nunique() == 1
        assert summary["n_het_snps"].nunique() == 1


class TestSaveHeterogeneitySummary:
    """Tests for save_heterogeneity with summary parameter."""

    def test_saves_summary_file(self, tmp_path, multi_cohort_locus_set):
        """Save with summary should create heterogeneity.txt.gz."""
        from credtools.meta import (
            compute_heterogeneity,
            heterogeneity_summary,
            save_heterogeneity,
        )

        het = compute_heterogeneity(multi_cohort_locus_set)
        summary = heterogeneity_summary(het, multi_cohort_locus_set)
        out_dir = str(tmp_path / "het_output")
        save_heterogeneity(het, out_dir, summary=summary)

        filepath = os.path.join(out_dir, "heterogeneity.txt.gz")
        assert os.path.exists(filepath)
        loaded = pd.read_csv(filepath, sep="\t", compression="gzip")
        assert len(loaded) == len(multi_cohort_locus_set.loci)

    def test_no_summary_no_file(self, tmp_path, multi_cohort_locus_set):
        """Save without summary should not create heterogeneity.txt.gz."""
        from credtools.meta import compute_heterogeneity, save_heterogeneity

        het = compute_heterogeneity(multi_cohort_locus_set)
        out_dir = str(tmp_path / "het_output")
        save_heterogeneity(het, out_dir)

        filepath = os.path.join(out_dir, "heterogeneity.txt.gz")
        assert not os.path.exists(filepath)


@pytest.fixture
def multi_cohort_same_pop_locus_set():
    """Create a LocusSet with EUR having 2 cohorts and AFR having 1 cohort."""
    locus_eur1 = _make_locus("EUR", "UKB", seed=42)
    locus_eur2 = _make_locus("EUR", "GWAS2", seed=99)
    locus_afr = _make_locus("AFR", "MVP", seed=123)
    return LocusSet([locus_eur1, locus_eur2, locus_afr])


@pytest.fixture
def two_pop_multi_cohort_locus_set():
    """Create a LocusSet where both EUR and AFR have 2 cohorts each."""
    locus_eur1 = _make_locus("EUR", "UKB", seed=42)
    locus_eur2 = _make_locus("EUR", "GWAS2", seed=99)
    locus_afr1 = _make_locus("AFR", "MVP", seed=123)
    locus_afr2 = _make_locus("AFR", "GWAS3", seed=456)
    return LocusSet([locus_eur1, locus_eur2, locus_afr1, locus_afr2])


class TestComputeHeterogeneityByPopulation:
    """Tests for compute_heterogeneity_by_population function."""

    def test_ld_metrics_computed_globally(self, multi_cohort_same_pop_locus_set):
        """LD metrics should be computed across all 3 cohorts."""
        from credtools.meta import compute_heterogeneity_by_population

        het = compute_heterogeneity_by_population(multi_cohort_same_pop_locus_set)

        assert "ld_4th_moment" in het
        assert "ld_decay" in het
        # All 3 cohorts should appear in ld_decay
        cohorts = het["ld_decay"]["cohort"].unique()
        assert len(cohorts) == 3

    def test_cochran_q_only_for_multi_cohort_populations(
        self, multi_cohort_same_pop_locus_set
    ):
        """Cochran Q should only be computed for populations with >1 cohort (EUR)."""
        from credtools.meta import compute_heterogeneity_by_population

        het = compute_heterogeneity_by_population(multi_cohort_same_pop_locus_set)

        assert "cochran_q" in het
        cq = het["cochran_q"]
        assert "population" in cq.columns
        # Only EUR should appear (AFR has only 1 cohort)
        assert set(cq["population"].unique()) == {"EUR"}

    def test_snp_missingness_only_same_population(
        self, multi_cohort_same_pop_locus_set
    ):
        """SNP missingness should only compare cohorts within the same population."""
        from credtools.meta import compute_heterogeneity_by_population

        het = compute_heterogeneity_by_population(multi_cohort_same_pop_locus_set)

        assert "snp_missingness" in het
        miss = het["snp_missingness"]
        assert "population" in miss.columns
        # Only EUR should appear
        assert set(miss["population"].unique()) == {"EUR"}

    def test_single_cohort_population_no_cochran_q(self, single_locus_set):
        """When all populations have only 1 cohort, cochran_q should not be produced."""
        from credtools.meta import compute_heterogeneity_by_population

        het = compute_heterogeneity_by_population(single_locus_set)

        assert "cochran_q" not in het
        assert "snp_missingness" not in het

    def test_summary_eur_has_q_afr_nan(self, multi_cohort_same_pop_locus_set):
        """In summary, EUR should have Q values, AFR should have NaN."""
        from credtools.meta import (
            compute_heterogeneity_by_population,
            heterogeneity_summary,
        )

        het = compute_heterogeneity_by_population(multi_cohort_same_pop_locus_set)
        summary = heterogeneity_summary(het, multi_cohort_same_pop_locus_set)

        eur_rows = summary[summary["popu"] == "EUR"]
        afr_rows = summary[summary["popu"] == "AFR"]

        # EUR should have non-NaN cochran_q values
        assert not eur_rows["cochran_q_median"].isna().any()
        # AFR should have NaN cochran_q values
        assert afr_rows["cochran_q_median"].isna().all()

    def test_two_populations_produce_different_q_values(
        self, two_pop_multi_cohort_locus_set
    ):
        """When both populations have multiple cohorts, each should have its own Q values."""
        from credtools.meta import (
            compute_heterogeneity_by_population,
            heterogeneity_summary,
        )

        het = compute_heterogeneity_by_population(two_pop_multi_cohort_locus_set)
        summary = heterogeneity_summary(het, two_pop_multi_cohort_locus_set)

        eur_q = summary[summary["popu"] == "EUR"]["cochran_q_median"].iloc[0]
        afr_q = summary[summary["popu"] == "AFR"]["cochran_q_median"].iloc[0]

        # Both should be non-NaN
        assert not pd.isna(eur_q)
        assert not pd.isna(afr_q)
        # They should differ (different seeds produce different data)
        assert eur_q != afr_q

    def test_save_snp_missingness_no_population_column(
        self, tmp_path, multi_cohort_same_pop_locus_set
    ):
        """Saved snp_missingness file should not contain population column."""
        from credtools.meta import (
            compute_heterogeneity_by_population,
            save_heterogeneity,
        )

        het = compute_heterogeneity_by_population(multi_cohort_same_pop_locus_set)
        out_dir = str(tmp_path / "het_output")
        save_heterogeneity(het, out_dir)

        filepath = os.path.join(out_dir, "snp_missingness.txt.gz")
        loaded = pd.read_csv(filepath, sep="\t", compression="gzip")
        assert "population" not in loaded.columns

    def test_save_cochran_q_has_population_column(
        self, tmp_path, multi_cohort_same_pop_locus_set
    ):
        """Saved cochran_q file should retain population column."""
        from credtools.meta import (
            compute_heterogeneity_by_population,
            save_heterogeneity,
        )

        het = compute_heterogeneity_by_population(multi_cohort_same_pop_locus_set)
        out_dir = str(tmp_path / "het_output")
        save_heterogeneity(het, out_dir)

        filepath = os.path.join(out_dir, "cochran_q.txt.gz")
        loaded = pd.read_csv(filepath, sep="\t", compression="gzip")
        assert "population" in loaded.columns


class TestMetaLocusIntegration:
    """Integration test for heterogeneity in meta_locus flow."""

    def test_meta_locus_outputs_heterogeneity(self, tmp_path, multi_cohort_locus_set):
        """Compute and save heterogeneity simulating the meta_locus flow."""
        from credtools.meta import (
            compute_heterogeneity,
            heterogeneity_summary,
            meta,
            save_heterogeneity,
        )

        # Compute heterogeneity BEFORE meta
        het = compute_heterogeneity(multi_cohort_locus_set)
        het_summary = heterogeneity_summary(het, multi_cohort_locus_set)
        het_summary["locus_id"] = "chr1_1000_3000"

        # Perform meta-analysis
        meta_result = meta(multi_cohort_locus_set, "meta_all")

        # Save heterogeneity with summary
        out_dir = str(tmp_path / "locus_output")
        save_heterogeneity(het, out_dir, summary=het_summary)

        # Verify detail files exist
        assert os.path.exists(os.path.join(out_dir, "ld_4th_moment.txt.gz"))
        assert os.path.exists(os.path.join(out_dir, "ld_decay.txt.gz"))
        assert os.path.exists(os.path.join(out_dir, "cochran_q.txt.gz"))
        assert os.path.exists(os.path.join(out_dir, "snp_missingness.txt.gz"))

        # Verify summary file exists and has locus_id
        summary_path = os.path.join(out_dir, "heterogeneity.txt.gz")
        assert os.path.exists(summary_path)
        loaded = pd.read_csv(summary_path, sep="\t", compression="gzip")
        assert "locus_id" in loaded.columns
        assert (loaded["locus_id"] == "chr1_1000_3000").all()

        # Verify meta produced a valid result
        assert meta_result.n_loci == 1
