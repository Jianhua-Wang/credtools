"""Tests for the SuSiE fine-mapping wrapper."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.susie import run_susie

from .conftest import _make_locus


def _mock_susie_rss_with_cs(**kwargs):
    """Mock susie_rss that returns results with credible sets."""
    n = len(kwargs["bhat"])
    pip = np.zeros(n)
    pip[0] = 0.8
    pip[1] = 0.15
    pip[2] = 0.05
    return {
        "pip": pip,
        "sets": {
            "cs": {"L1": [0, 1, 2]},
            "purity": {"min_abs_corr": [0.85]},
        },
        "converged": True,
        "niter": 12,
    }


def _mock_susie_rss_not_converged(**kwargs):
    """Mock susie_rss that did not converge but still returned a CS."""
    n = len(kwargs["bhat"])
    pip = np.zeros(n)
    pip[0] = 0.6
    pip[1] = 0.3
    return {
        "pip": pip,
        "sets": {
            "cs": {"L1": [0, 1]},
            "purity": {"min_abs_corr": [0.5]},
        },
        "converged": False,
        "niter": 100,
    }


def _mock_susie_rss_no_cs(**kwargs):
    """Mock susie_rss that returns results without credible sets."""
    n = len(kwargs["bhat"])
    pip = np.full(n, 1.0 / n)
    return {
        "pip": pip,
        "sets": {
            "cs": None,
            "purity": None,
        },
    }


def _mock_susie_rss_multi_cs(**kwargs):
    """Mock susie_rss that returns results with multiple credible sets."""
    n = len(kwargs["bhat"])
    pip = np.zeros(n)
    pip[0] = 0.7
    pip[1] = 0.2
    pip[5] = 0.6
    pip[6] = 0.3
    return {
        "pip": pip,
        "sets": {
            "cs": {"L1": [0, 1], "L2": [5, 6]},
            "purity": {"min_abs_corr": [0.9, 0.85]},
        },
    }


def _mock_susie_rss_no_purity(**kwargs):
    """Mock susie_rss that returns CS without purity info."""
    n = len(kwargs["bhat"])
    pip = np.zeros(n)
    pip[0] = 0.9
    return {
        "pip": pip,
        "sets": {
            "cs": {"L1": [0]},
        },
    }


class TestRunSusieBasic:
    """Basic SuSiE functionality tests."""

    def test_basic_call(self, locus_significant, monkeypatch):
        """Verify SuSiE returns a CredibleSet for significant locus."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_significant)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.SUSIE

    def test_output_structure(self, locus_significant, monkeypatch):
        """Verify SuSiE result has correct structure."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_significant)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1
        assert result.cs_sizes == [3]

    def test_parameters_stored(self, locus_significant, monkeypatch):
        """Parameters should be stored in the result."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(
            locus_significant,
            max_causal=3,
            coverage=0.99,
            max_iter=200,
        )
        assert result.parameters["max_causal"] == 3
        assert result.parameters["coverage"] == 0.99
        assert result.parameters["max_iter"] == 200


class TestRunSusieNoSignificant:
    """Tests for when no SNPs pass significance threshold."""

    def test_no_significant_snps(self, locus_no_significant):
        """Should return empty result (no mock needed, early return)."""
        result = run_susie(locus_no_significant)
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []
        assert (result.pips == 0).all()


class TestRunSusieUnmatched:
    """Tests for unmatched sumstats/LD."""

    def test_unmatched_ld_gets_intersected(self, locus_unmatched, monkeypatch):
        """Unmatched locus should be intersected before running."""
        assert not locus_unmatched.is_matched
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_unmatched)
        assert isinstance(result, CredibleSet)


class TestRunSusieCS:
    """Tests for credible set handling."""

    def test_no_cs_found(self, locus_significant, monkeypatch):
        """Should return n_cs=0 when susie_rss finds no credible sets."""
        monkeypatch.setattr("credtools.wrappers.susie.susie_rss", _mock_susie_rss_no_cs)
        result = run_susie(locus_significant)
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []

    def test_multiple_cs(self, locus_significant, monkeypatch):
        """Should correctly handle multiple credible sets."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_multi_cs
        )
        result = run_susie(locus_significant)
        assert result.n_cs == 2
        assert len(result.snps) == 2
        assert len(result.lead_snps) == 2

    def test_purity_extraction(self, locus_significant, monkeypatch):
        """Purity values should be extracted from susie_rss results."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_significant)
        assert result.purity is not None
        assert result.purity == [0.85]

    def test_purity_none_when_missing(self, locus_significant, monkeypatch):
        """Purity should be None when susie_rss doesn't provide it."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_no_purity
        )
        result = run_susie(locus_significant)
        assert result.purity is None

    def test_lead_snp_is_max_pip_in_cs(self, locus_significant, monkeypatch):
        """Lead SNP should be the one with highest PIP in CS."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_significant)
        cs_snps = result.snps[0]
        lead_snp = result.lead_snps[0]
        assert lead_snp == result.pips[cs_snps].idxmax()


class TestRunSusieConvergence:
    """Tests for convergence reporting and empty-on-nonconvergence."""

    def test_converged_propagated_when_true(self, locus_significant, monkeypatch):
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_significant)
        assert result.converged is True
        assert result.n_iter == 12

    def test_converged_propagated_when_false(self, locus_significant, monkeypatch):
        """Without empty_on_nonconvergence the CS is preserved but converged=False."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_not_converged
        )
        result = run_susie(locus_significant)
        assert result.converged is False
        assert result.n_iter == 100
        assert result.n_cs == 1  # CS still returned by default

    def test_empty_on_nonconvergence_zeros_cs(self, locus_significant, monkeypatch):
        """With empty_on_nonconvergence=True, non-converged run returns n_cs=0."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_not_converged
        )
        result = run_susie(locus_significant, empty_on_nonconvergence=True)
        assert result.converged is False
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []
        assert (result.pips == 0).all()

    def test_empty_on_nonconvergence_keeps_converged_cs(
        self, locus_significant, monkeypatch
    ):
        """When converged=True the flag has no effect — CS preserved."""
        monkeypatch.setattr(
            "credtools.wrappers.susie.susie_rss", _mock_susie_rss_with_cs
        )
        result = run_susie(locus_significant, empty_on_nonconvergence=True)
        assert result.converged is True
        assert result.n_cs == 1

    def test_no_significant_keeps_converged_none(self, locus_no_significant):
        """Early-return path (no significant SNPs) leaves converged=None."""
        result = run_susie(locus_no_significant)
        assert result.converged is None
        assert result.n_iter is None
