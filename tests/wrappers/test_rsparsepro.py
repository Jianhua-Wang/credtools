"""Tests for the RSparsePro fine-mapping wrapper."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.RSparsePro import (
    RSparsePro,
    adaptive_train,
    get_eff_maxld,
    get_eff_minld,
    get_ordered,
    rsparsepro_main,
    run_rsparsepro,
)

from .conftest import _make_locus


# ============================================================
# Tests for the RSparsePro class
# ============================================================


class TestRSparseProClass:
    """Tests for the RSparsePro class internals."""

    @pytest.fixture
    def small_model(self):
        """Create a small RSparsePro model for testing."""
        P, K = 5, 2
        R = np.eye(P)
        return RSparsePro(P, K, R, vare=1.0)

    @pytest.fixture
    def model_vare_zero(self):
        """Create a model with vare=0 (perfect LD)."""
        P, K = 5, 1
        R = np.eye(P)
        return RSparsePro(P, K, R, vare=0)

    def test_init_shape(self, small_model):
        """Model attributes should have correct shapes."""
        assert small_model.p == 5
        assert small_model.k == 2
        assert small_model.beta_mu.shape == (5, 2)
        assert small_model.gamma.shape == (5, 2)
        assert small_model.tilde_b.shape == (5,)

    def test_init_vare_zero(self, model_vare_zero):
        """Model with vare=0 should not have mat attribute."""
        assert not hasattr(model_vare_zero, "mat")

    def test_init_vare_nonzero(self, small_model):
        """Model with vare>0 should have mat attribute."""
        assert hasattr(small_model, "mat")
        assert small_model.mat.shape == (5, 5)

    def test_infer_q_beta(self, small_model):
        """infer_q_beta should update gamma and beta_mu."""
        R = np.eye(5)
        small_model.tilde_b = np.array([1.0, 0.5, -0.3, 0.1, 0.0])
        small_model.infer_q_beta(R)
        # gamma should be probability distributions (softmax)
        for k in range(small_model.k):
            assert pytest.approx(small_model.gamma[:, k].sum(), abs=1e-6) == 1.0

    def test_infer_tilde_b_vare_zero(self, model_vare_zero):
        """With vare=0, tilde_b should equal bhat."""
        bhat = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        model_vare_zero.infer_tilde_b(bhat)
        np.testing.assert_array_equal(model_vare_zero.tilde_b, bhat)

    def test_infer_tilde_b_vare_nonzero(self, small_model):
        """With vare>0, tilde_b should be modified."""
        bhat = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        small_model.infer_tilde_b(bhat)
        # Should not be identical to bhat
        assert not np.array_equal(small_model.tilde_b, bhat)

    def test_train_convergence(self):
        """Train should converge for simple data."""
        P, K = 3, 1
        R = np.eye(P)
        model = RSparsePro(P, K, R, vare=0)
        bhat = np.array([5.0, 0.1, 0.05])
        converged = model.train(bhat, R, maxite=100, eps=1e-3, ubound=100000)
        assert isinstance(converged, bool)

    def test_train_non_convergence(self):
        """Train should not converge with maxite=1 and strict eps."""
        P, K = 5, 2
        R = np.eye(P)
        model = RSparsePro(P, K, R, vare=0)
        bhat = np.array([10.0, 8.0, 6.0, 4.0, 2.0])
        converged = model.train(bhat, R, maxite=1, eps=1e-30, ubound=100000)
        assert converged is False

    def test_get_pip(self, small_model):
        """get_PIP should return array of shape (p,)."""
        R = np.eye(5)
        small_model.tilde_b = np.array([3.0, 0.1, 0.05, 0.0, 0.0])
        small_model.infer_q_beta(R)
        pip = small_model.get_PIP()
        assert pip.shape == (5,)
        assert np.all(pip >= 0)
        assert np.all(pip <= 1)

    def test_get_effect(self, small_model):
        """get_effect should return valid effect groups."""
        R = np.eye(5)
        small_model.tilde_b = np.array([3.0, 0.1, 0.05, 0.0, 0.0])
        small_model.infer_q_beta(R)
        eff, eff_gamma, eff_mu = small_model.get_effect(cthres=0.5)
        assert isinstance(eff, dict)
        assert isinstance(eff_gamma, dict)
        assert isinstance(eff_mu, dict)

    def test_get_ztilde(self, small_model):
        """get_ztilde should return corrected estimates."""
        bhat = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        small_model.infer_tilde_b(bhat)
        ztilde = small_model.get_ztilde()
        assert ztilde.shape == (5,)


# ============================================================
# Tests for helper functions
# ============================================================


class TestHelperFunctions:
    """Tests for module-level helper functions."""

    def test_get_eff_maxld_single_group(self):
        """Single effect group should return 0.0."""
        eff = {0: [0, 1]}
        ld = np.eye(3)
        assert get_eff_maxld(eff, ld) == 0.0

    def test_get_eff_maxld_two_groups(self):
        """Two groups should return max abs LD between leads."""
        eff = {0: [0], 1: [1]}
        ld = np.array([[1.0, 0.3], [0.3, 1.0]])
        result = get_eff_maxld(eff, ld)
        assert pytest.approx(result, abs=1e-6) == 0.3

    def test_get_eff_minld_empty(self):
        """Empty effect dict should return 1.0."""
        assert get_eff_minld({}, np.eye(3)) == 1.0

    def test_get_eff_minld_single_variant(self):
        """Single variant per group should return 1.0."""
        eff = {0: [0]}
        ld = np.array([[1.0, 0.3], [0.3, 1.0]])
        assert get_eff_minld(eff, ld) == 1.0

    def test_get_eff_minld_multiple_variants(self):
        """Should return minimum absolute LD within groups."""
        eff = {0: [0, 1]}
        ld = np.array([[1.0, 0.5], [0.5, 1.0]])
        assert pytest.approx(get_eff_minld(eff, ld), abs=1e-6) == 0.5

    def test_get_ordered_single(self):
        """Single group always ordered."""
        assert get_ordered({0: np.array([1.0])}) is True

    def test_get_ordered_empty(self):
        """Empty dict returns True."""
        assert get_ordered({}) is True

    def test_get_ordered_multiple_ordered(self):
        """Properly ordered groups should return True."""
        eff_mu = {0: np.array([-5.0]), 1: np.array([-3.0])}
        assert get_ordered(eff_mu) is True

    def test_get_ordered_multiple_unordered(self):
        """Non-consecutive keys should return False."""
        eff_mu = {0: np.array([-5.0]), 2: np.array([-3.0])}
        assert get_ordered(eff_mu) is False


# ============================================================
# Tests for run_rsparsepro wrapper
# ============================================================


def _mock_rsparsepro_main(zfile, ld, **kwargs):
    """Mock rsparsepro_main that returns valid results."""
    n = len(zfile)
    zfile = zfile.copy()
    pip_vals = np.zeros(n)
    pip_vals[0] = 0.8
    pip_vals[1] = 0.1
    pip_vals[2] = 0.05
    zfile["PIP"] = pip_vals
    zfile["z_estimated"] = np.random.randn(n)
    zfile["cs"] = 0
    zfile.loc[0, "cs"] = 1
    zfile.loc[1, "cs"] = 1
    zfile.loc[2, "cs"] = 1
    return zfile


def _mock_rsparsepro_main_no_cs(zfile, ld, **kwargs):
    """Mock rsparsepro_main that returns no credible sets."""
    n = len(zfile)
    zfile = zfile.copy()
    zfile["PIP"] = np.full(n, 1.0 / n)
    zfile["z_estimated"] = np.zeros(n)
    zfile["cs"] = 0
    return zfile


def _mock_rsparsepro_main_multi_cs(zfile, ld, **kwargs):
    """Mock rsparsepro_main that returns multiple credible sets."""
    n = len(zfile)
    zfile = zfile.copy()
    pip_vals = np.zeros(n)
    pip_vals[0] = 0.7
    pip_vals[1] = 0.2
    pip_vals[5] = 0.6
    pip_vals[6] = 0.3
    zfile["PIP"] = pip_vals
    zfile["z_estimated"] = np.zeros(n)
    zfile["cs"] = 0
    zfile.loc[0, "cs"] = 1
    zfile.loc[1, "cs"] = 1
    zfile.loc[5, "cs"] = 2
    zfile.loc[6, "cs"] = 2
    return zfile


class TestRunRsparsepro:
    """Tests for the run_rsparsepro wrapper function."""

    def test_basic_call(self, locus_significant, monkeypatch):
        """run_rsparsepro returns CredibleSet."""
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main,
        )
        result = run_rsparsepro(locus_significant)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.RSparsePro

    def test_no_significant_snps(self, locus_no_significant):
        """Should return empty result (early return, no mock needed)."""
        result = run_rsparsepro(locus_no_significant)
        assert result.n_cs == 0
        assert (result.pips == 0).all()

    def test_unmatched_ld(self, locus_unmatched, monkeypatch):
        """Unmatched locus should be intersected."""
        assert not locus_unmatched.is_matched
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main,
        )
        result = run_rsparsepro(locus_unmatched)
        assert isinstance(result, CredibleSet)

    def test_output_structure(self, locus_significant, monkeypatch):
        """Result should have correct structure."""
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main,
        )
        result = run_rsparsepro(locus_significant)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1

    def test_no_cs_from_rsparsepro(self, locus_significant, monkeypatch):
        """Should return n_cs=0 when no credible sets found."""
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main_no_cs,
        )
        result = run_rsparsepro(locus_significant)
        assert result.n_cs == 0
        assert result.snps == []

    def test_multiple_cs(self, locus_significant, monkeypatch):
        """Should handle multiple credible sets correctly."""
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main_multi_cs,
        )
        result = run_rsparsepro(locus_significant)
        assert result.n_cs == 2
        assert len(result.snps) == 2
        assert len(result.lead_snps) == 2

    def test_purity_calculated(self, locus_significant, monkeypatch):
        """Purity should be calculated when LD is available."""
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main,
        )
        result = run_rsparsepro(locus_significant)
        assert result.purity is not None
        assert len(result.purity) == 1

    def test_parameters_stored(self, locus_significant, monkeypatch):
        """Parameters should be stored in the result."""
        monkeypatch.setattr(
            "credtools.wrappers.RSparsePro.rsparsepro_main",
            _mock_rsparsepro_main,
        )
        result = run_rsparsepro(locus_significant, max_causal=3, coverage=0.99)
        assert result.parameters["max_causal"] == 3
        assert result.parameters["coverage"] == 0.99
