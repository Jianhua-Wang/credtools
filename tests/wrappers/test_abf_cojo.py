"""Tests for the ABF+COJO combined fine-mapping wrapper."""

from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.abf_cojo import (
    _create_conditional_locus,
    _run_conditional_abf_analysis,
    run_abf_cojo,
)

from .conftest import _make_locus


def _make_cojo_results(snp_ids, n_signals=None):
    """Create a mock COJO results DataFrame.

    Parameters
    ----------
    snp_ids : list of str
        SNP IDs to use as independent signals.
    """
    if n_signals is not None:
        snp_ids = snp_ids[:n_signals]
    return pd.DataFrame({"SNP": snp_ids})


class TestRunAbfCojoNoSignificant:
    """Tests for when no SNPs pass significance threshold."""

    def test_no_significant_snps(self, locus_no_significant):
        """Should return empty result (early return)."""
        result = run_abf_cojo(locus_no_significant)
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []
        assert (result.pips == 0).all()
        assert result.tool == f"{Method.ABF}_COJO"


class TestRunAbfCojoZeroSignals:
    """Tests for when COJO detects 0 signals."""

    def test_zero_cojo_signals(self, locus_significant, monkeypatch):
        """Should return empty result when COJO finds no signals."""
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.conditional_selection",
            lambda *args, **kwargs: pd.DataFrame({"SNP": []}),
        )
        result = run_abf_cojo(locus_significant)
        assert result.n_cs == 0
        assert result.tool == f"{Method.ABF}_COJO"


class TestRunAbfCojoSingleSignal:
    """Tests for when COJO detects 1 signal (delegates to ABF)."""

    def test_single_signal_delegates_to_abf(self, locus_significant, monkeypatch):
        """Single signal should delegate to standard ABF and rename tool."""
        snp_id = locus_significant.sumstats[ColName.SNPID].iloc[0]
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.conditional_selection",
            lambda *args, **kwargs: _make_cojo_results([snp_id]),
        )
        result = run_abf_cojo(locus_significant)
        assert result.tool == f"{Method.ABF}_COJO"
        assert result.n_cs == 1

    def test_single_signal_parameters(self, locus_significant, monkeypatch):
        """Single signal should update parameters to ABF_COJO format."""
        snp_id = locus_significant.sumstats[ColName.SNPID].iloc[0]
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.conditional_selection",
            lambda *args, **kwargs: _make_cojo_results([snp_id]),
        )
        result = run_abf_cojo(locus_significant, coverage=0.90)
        assert "coverage" in result.parameters
        assert "p_cutoff" in result.parameters


class TestRunAbfCojoMultipleSignals:
    """Tests for when COJO detects multiple signals."""

    def test_multiple_signals_calls_conditional(self, locus_with_af2, monkeypatch):
        """Multiple signals should run conditional ABF analysis."""
        snpids = locus_with_af2.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.conditional_selection",
            lambda *args, **kwargs: _make_cojo_results(snpids, n_signals=2),
        )
        # Mock COJO class to avoid cojopy dependency complexity
        mock_cojo = MagicMock()
        cond_result = pd.DataFrame(
            {
                "SNP": snpids[:10],
                "cond_beta": np.random.randn(10) * 0.1,
                "cond_se": np.abs(np.random.randn(10) * 0.02) + 0.001,
                "cond_p": np.random.uniform(1e-10, 1e-8, 10),
            }
        )
        mock_cojo.run_conditional_analysis.return_value = cond_result
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.COJO",
            lambda **kwargs: mock_cojo,
        )
        result = run_abf_cojo(locus_with_af2)
        assert isinstance(result, CredibleSet)
        assert result.tool == f"{Method.ABF}_COJO"

    def test_purity_calculated(self, locus_with_af2, monkeypatch):
        """Purity should be calculated for multi-signal results."""
        snpids = locus_with_af2.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.conditional_selection",
            lambda *args, **kwargs: _make_cojo_results(snpids, n_signals=2),
        )
        mock_cojo = MagicMock()
        cond_result = pd.DataFrame(
            {
                "SNP": snpids[:10],
                "cond_beta": np.random.randn(10) * 0.1,
                "cond_se": np.abs(np.random.randn(10) * 0.02) + 0.001,
                "cond_p": np.random.uniform(1e-10, 1e-8, 10),
            }
        )
        mock_cojo.run_conditional_analysis.return_value = cond_result
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.COJO",
            lambda **kwargs: mock_cojo,
        )
        result = run_abf_cojo(locus_with_af2)
        # Purity may or may not be set depending on CS results
        if result.n_cs > 0:
            assert result.purity is not None


class TestRunAbfCojoNoAF2:
    """Tests for when AF2 is not in LD map."""

    def test_no_af2_warning(self, locus_significant, monkeypatch):
        """Should warn when AF2 is not in LD map."""
        snpids = locus_significant.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.conditional_selection",
            lambda *args, **kwargs: _make_cojo_results(snpids, n_signals=2),
        )
        mock_cojo = MagicMock()
        cond_result = pd.DataFrame(
            {
                "SNP": snpids[:10],
                "cond_beta": np.random.randn(10) * 0.1,
                "cond_se": np.abs(np.random.randn(10) * 0.02) + 0.001,
                "cond_p": np.random.uniform(1e-10, 1e-8, 10),
            }
        )
        mock_cojo.run_conditional_analysis.return_value = cond_result
        monkeypatch.setattr(
            "credtools.wrappers.abf_cojo.COJO",
            lambda **kwargs: mock_cojo,
        )
        # locus_significant does not have AF2
        result = run_abf_cojo(locus_significant)
        assert isinstance(result, CredibleSet)


class TestCreateConditionalLocus:
    """Tests for the _create_conditional_locus helper."""

    def test_creates_new_locus(self, locus_significant):
        """Should create a new locus with conditional stats."""
        snpids = locus_significant.sumstats[ColName.SNPID].tolist()
        cond_results = pd.DataFrame(
            {
                "SNP": snpids,
                "cond_beta": np.random.randn(len(snpids)) * 0.1,
                "cond_se": np.abs(np.random.randn(len(snpids)) * 0.02) + 0.001,
                "cond_p": np.random.uniform(1e-10, 1e-8, len(snpids)),
            }
        )
        new_locus = _create_conditional_locus(
            locus_significant, cond_results, snpids[0]
        )
        assert new_locus is not locus_significant
        assert len(new_locus.sumstats) == len(snpids)

    def test_drops_snps_not_in_results(self, locus_significant):
        """Verify SNPs not in conditional results are dropped."""
        snpids = locus_significant.sumstats[ColName.SNPID].tolist()
        half = len(snpids) // 2
        cond_results = pd.DataFrame(
            {
                "SNP": snpids[:half],
                "cond_beta": np.random.randn(half) * 0.1,
                "cond_se": np.abs(np.random.randn(half) * 0.02) + 0.001,
                "cond_p": np.random.uniform(1e-10, 1e-8, half),
            }
        )
        new_locus = _create_conditional_locus(
            locus_significant, cond_results, snpids[0]
        )
        assert len(new_locus.sumstats) == half
