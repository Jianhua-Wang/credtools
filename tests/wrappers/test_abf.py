"""Tests for the ABF fine-mapping wrapper."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.abf import run_abf

from .conftest import _make_locus


class TestRunAbfBasic:
    """Basic ABF functionality tests."""

    def test_basic_call(self, locus_significant):
        """ABF returns a CredibleSet for significant locus."""
        result = run_abf(locus_significant)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.ABF

    def test_output_structure(self, locus_significant):
        """ABF result has correct structure."""
        result = run_abf(locus_significant)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1
        assert len(result.cs_sizes) == 1
        assert result.coverage == 0.95

    def test_pip_sum_close_to_one(self, locus_significant):
        """PIP values should sum to approximately 1."""
        result = run_abf(locus_significant)
        assert abs(result.pips.sum() - 1.0) < 1e-6

    def test_cs_coverage_met(self, locus_significant):
        """Credible set should meet the requested coverage."""
        coverage = 0.95
        result = run_abf(locus_significant, coverage=coverage)
        cs_snps = result.snps[0]
        cs_pip_sum = result.pips[cs_snps].sum()
        assert cs_pip_sum >= coverage

    def test_parameters_stored(self, locus_significant):
        """Parameters should be stored in the result."""
        result = run_abf(locus_significant, coverage=0.90, var_prior=0.15)
        assert result.parameters["coverage"] == 0.90
        assert result.parameters["var_prior"] == 0.15
        assert result.parameters["max_causal"] == 1


class TestRunAbfEdgeCases:
    """Edge case tests for ABF."""

    def test_no_significant_snps(self, locus_no_significant):
        """Should return empty result when no SNPs pass threshold."""
        result = run_abf(locus_no_significant)
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []
        assert result.cs_sizes == []
        assert (result.pips == 0).all()

    def test_max_causal_greater_than_one_warns(self, locus_significant):
        """Should warn and reset max_causal to 1 when max_causal > 1."""
        result = run_abf(locus_significant, max_causal=5)
        assert result.parameters["max_causal"] == 1
        assert result.n_cs == 1

    def test_custom_var_prior(self, locus_significant):
        """Different var_prior should produce different results."""
        r1 = run_abf(locus_significant, var_prior=0.15)
        r2 = run_abf(locus_significant, var_prior=0.5)
        # PIPs should differ with different priors (use strict rtol)
        assert not np.allclose(r1.pips.values, r2.pips.values, atol=0, rtol=1e-2)

    def test_custom_coverage(self, locus_significant):
        """Higher coverage should produce larger credible sets."""
        r_low = run_abf(locus_significant, coverage=0.50)
        r_high = run_abf(locus_significant, coverage=0.99)
        assert r_low.cs_sizes[0] <= r_high.cs_sizes[0]

    def test_lead_snp_has_smallest_p_in_cs(self, locus_significant):
        """Lead SNP should be the one with smallest p-value in the CS."""
        result = run_abf(locus_significant)
        cs_snps = result.snps[0]
        # lead_snps is a list of lists for ABF; extract the inner list
        lead_snp = result.lead_snps[0]
        if isinstance(lead_snp, list):
            lead_snp = lead_snp[0]
        sumstats = locus_significant.original_sumstats
        cs_mask = sumstats[ColName.SNPID].isin(cs_snps)
        expected_lead = sumstats.loc[sumstats[cs_mask][ColName.P].idxmin(), ColName.SNPID]
        assert lead_snp == expected_lead

    def test_purity_calculated_with_ld(self, locus_significant):
        """Purity should be calculated when LD matrix is available."""
        result = run_abf(locus_significant)
        assert result.purity is not None
        assert len(result.purity) == 1
        assert 0.0 <= result.purity[0] <= 1.0

    def test_single_snp_locus(self):
        """ABF should work with a single SNP locus."""
        locus = _make_locus(n_snps=1, p_range=(1e-12, 1e-11))
        result = run_abf(locus)
        assert result.n_cs == 1
        assert result.pips.iloc[0] == pytest.approx(1.0, abs=1e-6)

    def test_numerical_stability_small_se(self):
        """ABF should handle very small SE values without overflow."""
        locus = _make_locus(n_snps=10, p_range=(1e-12, 1e-9), seed=99)
        # Set very small SE values
        locus.sumstats[ColName.SE] = 1e-15
        locus._original_sumstats[ColName.SE] = 1e-15
        result = run_abf(locus)
        assert not np.any(np.isinf(result.pips.values))
        assert not np.any(np.isnan(result.pips.values))

    def test_pip_index_matches_snpids(self, locus_significant):
        """PIP index should match locus SNPID list."""
        result = run_abf(locus_significant)
        expected_snpids = locus_significant.original_sumstats[ColName.SNPID].tolist()
        assert list(result.pips.index) == expected_snpids

    def test_custom_significance_threshold(self):
        """Custom significance threshold should be respected."""
        # Create locus with p-values around 1e-5
        locus = _make_locus(p_range=(1e-6, 1e-5))
        # Default threshold (5e-8) should find nothing
        result_strict = run_abf(locus)
        assert result_strict.n_cs == 0
        # Relaxed threshold should find something
        result_relaxed = run_abf(locus, significant_threshold=1e-4)
        assert result_relaxed.n_cs == 1
