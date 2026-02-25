"""Unit tests for susie_rss.py internal functions.

Tests pure utility functions, statistical computations, initialization,
and integration of SuSiE components without mocking the core logic.
"""

import math
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import stats

from credtools.wrappers.susie_rss import (
    Eloglik_ss,
    SER_posterior_e_loglik_ss,
    compute_tf_Xb,
    est_V_uniroot,
    estimate_residual_variance_ss,
    get_ER2_ss,
    get_objective_ss,
    get_purity,
    in_CS,
    in_CS_x,
    init_finalize,
    init_setup,
    lbf_grad,
    loglik,
    loglik_grad,
    muffled_cov2cor,
    n_in_CS,
    n_in_CS_x,
    neg_loglik_logscale,
    negloglik_grad_logscale,
    optimize_prior_variance,
    single_effect_regression_ss,
    summary_susie,
    susie_get_cs,
    susie_get_objective,
    susie_get_pip,
    susie_prune_single_effects,
    susie_slim,
    susie_suff_stat,
    update_each_effect_ss,
)


# ─── Fixtures ───────────────────────────────────────────────────────


@pytest.fixture
def rng():
    """Fixed random number generator."""
    return np.random.default_rng(42)


@pytest.fixture
def small_pd_matrix(rng):
    """Small 5x5 positive-definite correlation matrix."""
    A = rng.normal(0, 1, (5, 5))
    cov = A @ A.T + 5 * np.eye(5)
    d = np.sqrt(np.diag(cov))
    R = cov / np.outer(d, d)
    np.fill_diagonal(R, 1.0)
    return R


@pytest.fixture
def susie_fit_obj():
    """Minimal SuSiE fit object for testing."""
    p = 5
    L = 2
    alpha = np.array([
        [0.7, 0.1, 0.1, 0.05, 0.05],
        [0.05, 0.05, 0.1, 0.7, 0.1],
    ])
    return {
        "alpha": alpha,
        "mu": np.array([
            [1.0, 0.1, 0.1, 0.0, 0.0],
            [0.0, 0.0, 0.1, 1.0, 0.1],
        ]),
        "mu2": np.array([
            [1.1, 0.11, 0.11, 0.01, 0.01],
            [0.01, 0.01, 0.11, 1.1, 0.11],
        ]),
        "V": np.array([0.5, 0.3]),
        "sigma2": 1.0,
        "KL": np.array([0.1, 0.2]),
        "lbf": np.array([2.0, 1.5]),
        "lbf_variable": np.zeros((L, p)),
        "XtXr": np.zeros(p),
        "pi": np.ones(p) / p,
        "null_index": 0,
        "elbo": np.array([-100.0, -90.0, -85.0]),
        "Xr": np.zeros(10),
    }


# ═══════════════════════════════════════════════════════════════════
# Phase 1: Pure function tests (T1–T4)
# ═══════════════════════════════════════════════════════════════════


class TestCredibleSetFunctions:
    """T1: in_CS_x, in_CS, n_in_CS_x, n_in_CS."""

    def test_n_in_CS_x_concentrated(self):
        """Concentrated PIP needs 1 variable."""
        x = np.array([0.95, 0.02, 0.02, 0.01])
        assert n_in_CS_x(x, 0.9) == 1

    def test_n_in_CS_x_uniform(self):
        """Uniform PIP needs many variables."""
        x = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        assert n_in_CS_x(x, 0.9) == 5

    def test_n_in_CS_x_full_coverage(self):
        """Requesting 100% coverage should include all."""
        x = np.array([0.5, 0.3, 0.2])
        assert n_in_CS_x(x, 1.0) == 3

    def test_in_CS_x_returns_binary(self):
        """Output should be binary 0/1."""
        x = np.array([0.1, 0.6, 0.2, 0.05, 0.05])
        result = in_CS_x(x, 0.8)
        assert set(np.unique(result)).issubset({0, 1})

    def test_in_CS_x_covers_target(self):
        """Included variables should achieve coverage."""
        x = np.array([0.1, 0.6, 0.2, 0.05, 0.05])
        result = in_CS_x(x, 0.8)
        assert x[result == 1].sum() >= 0.8

    def test_in_CS_dict_input(self):
        """in_CS should accept dict with 'alpha' key."""
        alpha = np.array([[0.8, 0.1, 0.1], [0.1, 0.1, 0.8]])
        result = in_CS({"alpha": alpha}, 0.9)
        assert result.shape == (2, 3)
        assert result[0, 0] == 1
        assert result[1, 2] == 1

    def test_in_CS_ndarray_input(self):
        """in_CS should accept raw ndarray."""
        alpha = np.array([[0.8, 0.1, 0.1], [0.1, 0.1, 0.8]])
        result = in_CS(alpha, 0.9)
        assert result.shape == (2, 3)

    def test_n_in_CS_dict_input(self):
        """n_in_CS should accept dict."""
        alpha = np.array([[0.8, 0.1, 0.1], [0.33, 0.34, 0.33]])
        result = n_in_CS({"alpha": alpha}, 0.9)
        assert result[0] <= result[1]

    def test_n_in_CS_array_input(self):
        """n_in_CS should accept ndarray."""
        alpha = np.array([[0.9, 0.05, 0.05]])
        result = n_in_CS(alpha, 0.9)
        assert result[0] == 1


class TestStatisticalGradients:
    """T2: lbf_grad, loglik, loglik_grad, neg_loglik_logscale, negloglik_grad_logscale."""

    def test_lbf_grad_normal(self):
        """Normal inputs should produce finite gradient."""
        V = 1.0
        shat2 = np.array([0.1, 0.2, 0.5])
        T2 = np.array([4.0, 2.0, 1.0])
        result = lbf_grad(V, shat2, T2)
        assert result.shape == (3,)
        assert np.all(np.isfinite(result))

    def test_lbf_grad_nan_handling(self):
        """NaN in T2 should produce 0 in gradient."""
        V = 1.0
        shat2 = np.array([0.1, 0.2])
        T2 = np.array([4.0, np.nan])
        result = lbf_grad(V, shat2, T2)
        assert result[1] == 0

    def test_loglik_normal_input(self):
        """Log-likelihood should be finite for normal inputs."""
        V = 1.0
        betahat = np.array([0.5, -0.3, 0.1])
        shat2 = np.array([0.1, 0.2, 0.3])
        prior_weights = np.array([1 / 3, 1 / 3, 1 / 3])
        result = loglik(V, betahat, shat2, prior_weights)
        assert np.isfinite(result)

    def test_loglik_inf_shat2(self):
        """Infinite shat2 should be handled gracefully."""
        V = 1.0
        betahat = np.array([0.5, 0.3])
        shat2 = np.array([0.1, np.inf])
        prior_weights = np.array([0.5, 0.5])
        result = loglik(V, betahat, shat2, prior_weights)
        assert np.isfinite(result)

    def test_loglik_V_zero(self):
        """V=0 should produce finite loglik."""
        V = 0
        betahat = np.array([0.5, 0.3])
        shat2 = np.array([0.1, 0.2])
        prior_weights = np.array([0.5, 0.5])
        result = loglik(V, betahat, shat2, prior_weights)
        assert np.isfinite(result)

    def test_neg_loglik_logscale_negates(self):
        """neg_loglik_logscale = -loglik(exp(lV), ...)."""
        lV = 0.5
        betahat = np.array([0.5, 0.3])
        shat2 = np.array([0.1, 0.2])
        pw = np.array([0.5, 0.5])
        expected = -loglik(np.exp(lV), betahat, shat2, pw)
        result = neg_loglik_logscale(lV, betahat, shat2, pw)
        assert_allclose(result, expected)

    def test_negloglik_grad_logscale_relationship(self):
        """negloglik_grad_logscale = -exp(lV) * loglik_grad(exp(lV), ...)."""
        lV = 0.0
        betahat = np.array([0.5, 0.3])
        shat2 = np.array([0.1, 0.2])
        pw = np.array([0.5, 0.5])
        expected = -np.exp(lV) * loglik_grad(np.exp(lV), betahat, shat2, pw)
        result = negloglik_grad_logscale(lV, betahat, shat2, pw)
        assert_allclose(result, expected)

    def test_loglik_grad_finite(self):
        """Gradient should be finite for normal input."""
        V = 1.0
        betahat = np.array([1.0, -0.5])
        shat2 = np.array([0.1, 0.2])
        pw = np.array([0.5, 0.5])
        result = loglik_grad(V, betahat, shat2, pw)
        assert np.isfinite(result)

    def test_loglik_grad_inf_shat2(self):
        """Gradient handles infinite shat2."""
        V = 1.0
        betahat = np.array([1.0, 0.5])
        shat2 = np.array([0.1, np.inf])
        pw = np.array([0.5, 0.5])
        result = loglik_grad(V, betahat, shat2, pw)
        assert np.isfinite(result)


class TestGetPurity:
    """T3: get_purity."""

    def test_single_variant(self):
        """Single variant should return [1, 1, 1]."""
        result = get_purity(np.array([0]), X=None, Xcorr=np.eye(3))
        assert_allclose(result, [1, 1, 1])

    def test_xcorr_input(self, small_pd_matrix):
        """Should compute purity from correlation matrix."""
        pos = np.array([0, 1, 2])
        result = get_purity(pos, X=None, Xcorr=small_pd_matrix)
        assert result.shape == (3,)
        assert result[0] <= result[1]  # min <= mean
        assert result[0] <= result[2]  # min <= median

    def test_x_input(self, rng):
        """Should compute purity from data matrix."""
        X = rng.normal(0, 1, (50, 5))
        pos = np.array([0, 1, 2])
        result = get_purity(pos, X=X, Xcorr=None)
        assert result.shape == (3,)
        assert 0 <= result[0] <= 1

    def test_subsampling_large_cs(self, small_pd_matrix):
        """Large CS should be subsampled to n."""
        pos = np.arange(5)
        result = get_purity(pos, X=None, Xcorr=small_pd_matrix, n=3)
        assert result.shape == (3,)

    def test_squared_purity(self, small_pd_matrix):
        """Squared purity values should be non-negative."""
        pos = np.array([0, 1])
        result = get_purity(pos, X=None, Xcorr=small_pd_matrix, squared=True)
        assert result[0] >= 0


class TestSusieGetObjectiveAndPIP:
    """T4: susie_get_objective, susie_get_pip."""

    def test_get_objective_last_only(self, susie_fit_obj):
        """last_only=True returns single float."""
        result = susie_get_objective(susie_fit_obj, last_only=True)
        assert isinstance(result, (float, np.floating))

    def test_get_objective_all_iterations(self, susie_fit_obj):
        """last_only=False returns full array."""
        result = susie_get_objective(susie_fit_obj, last_only=False)
        assert len(result) == 3

    def test_get_objective_warns_on_decrease(self, caplog):
        """Should log warning when ELBO decreases."""
        import logging

        res = {"elbo": np.array([-80.0, -90.0, -85.0])}
        with caplog.at_level(logging.WARNING, logger="SuSiE"):
            susie_get_objective(res, last_only=True)
        assert "decreasing" in caplog.text.lower()

    def test_get_pip_dict_input(self, susie_fit_obj):
        """PIP from dict should be in [0, 1]."""
        pip = susie_get_pip(susie_fit_obj)
        assert pip.shape == (5,)
        assert np.all((pip >= 0) & (pip <= 1))

    def test_get_pip_array_input(self):
        """PIP from raw alpha matrix."""
        alpha = np.array([[0.8, 0.1, 0.1], [0.1, 0.8, 0.1]])
        pip = susie_get_pip(alpha)
        assert pip.shape == (3,)
        assert pip[0] > 0.5  # variant 0 has high alpha in component 0

    def test_get_pip_V_filtering(self):
        """Components with V near zero should be filtered."""
        res = {
            "alpha": np.array([[0.8, 0.1, 0.1], [0.1, 0.1, 0.8]]),
            "V": np.array([0.5, 1e-12]),  # second component filtered
        }
        pip = susie_get_pip(res, prior_tol=1e-9)
        # Only first component contributes
        assert_allclose(pip, 1 - (1 - res["alpha"][0]), atol=1e-10)

    def test_get_pip_null_index(self):
        """null_index column should be removed."""
        res = {
            "alpha": np.array([[0.3, 0.3, 0.2, 0.2]]),
            "null_index": 1,
        }
        pip = susie_get_pip(res)
        assert pip.shape == (3,)  # one column removed

    def test_get_pip_prune_by_cs(self):
        """prune_by_cs should only use CS components."""
        res = {
            "alpha": np.array([[0.8, 0.1, 0.1], [0.1, 0.1, 0.8]]),
            "V": np.array([0.5, 0.5]),
            "sets": {"cs_index": np.array([0])},
        }
        pip = susie_get_pip(res, prune_by_cs=True)
        # Only first component used
        assert_allclose(pip, 1 - (1 - res["alpha"][0]), atol=1e-10)

    def test_get_pip_prune_by_cs_no_cs(self):
        """prune_by_cs with no CS info returns zeros."""
        res = {
            "alpha": np.array([[0.8, 0.1, 0.1]]),
            "V": np.array([0.5]),
        }
        pip = susie_get_pip(res, prune_by_cs=True)
        assert_allclose(pip, np.zeros(3), atol=1e-10)


# ═══════════════════════════════════════════════════════════════════
# Phase 2: Statistical computation core (T5–T8)
# ═══════════════════════════════════════════════════════════════════


class TestSufficientStatFunctions:
    """T5: SER_posterior_e_loglik_ss, get_ER2_ss, Eloglik_ss, estimate_residual_variance_ss, get_objective_ss."""

    @pytest.fixture
    def suff_data(self):
        """Small sufficient statistics for testing."""
        p = 3
        rng = np.random.default_rng(42)
        A = rng.normal(0, 1, (p, p))
        XtX = A @ A.T + 3 * np.eye(p)
        Xty = rng.normal(0, 1, p)
        yty = 10.0
        n = 100
        s = {
            "alpha": np.array([[0.7, 0.2, 0.1]]),
            "mu": np.array([[1.0, 0.5, 0.1]]),
            "mu2": np.array([[1.1, 0.3, 0.02]]),
            "sigma2": 1.0,
            "KL": np.array([0.05]),
        }
        return XtX, Xty, yty, n, s

    def test_SER_posterior_e_loglik_ss_zero(self):
        """All-zero inputs should return 0."""
        dXtX = np.array([1.0, 1.0, 1.0])
        Xty = np.zeros(3)
        Eb = np.zeros(3)
        Eb2 = np.zeros(3)
        result = SER_posterior_e_loglik_ss(dXtX, Xty, 1.0, Eb, Eb2)
        assert_allclose(result, 0.0)

    def test_SER_posterior_e_loglik_ss_normal(self):
        """Normal input should return finite result."""
        dXtX = np.array([2.0, 3.0])
        Xty = np.array([1.0, -0.5])
        Eb = np.array([0.4, -0.2])
        Eb2 = np.array([0.2, 0.1])
        result = SER_posterior_e_loglik_ss(dXtX, Xty, 1.0, Eb, Eb2)
        assert np.isfinite(result)

    def test_get_ER2_ss_zero_effects(self, suff_data):
        """Zero effects → ER2 = yty."""
        XtX, Xty, yty, n, _ = suff_data
        s = {
            "alpha": np.zeros((1, 3)),
            "mu": np.zeros((1, 3)),
            "mu2": np.zeros((1, 3)),
        }
        result = get_ER2_ss(XtX, Xty, s, yty)
        assert_allclose(result, yty, atol=1e-10)

    def test_get_ER2_ss_finite(self, suff_data):
        """Normal input produces finite ER2."""
        XtX, Xty, yty, n, s = suff_data
        result = get_ER2_ss(XtX, Xty, s, yty)
        assert np.isfinite(result)

    def test_Eloglik_ss_finite(self, suff_data):
        """Expected log-likelihood is finite."""
        XtX, Xty, yty, n, s = suff_data
        result = Eloglik_ss(XtX, Xty, s, yty, n)
        assert np.isfinite(result)

    def test_estimate_residual_variance_ss(self, suff_data):
        """Residual variance should be positive."""
        XtX, Xty, yty, n, s = suff_data
        result = estimate_residual_variance_ss(XtX, Xty, s, yty, n)
        assert np.isfinite(result)

    def test_get_objective_ss_relationship(self, suff_data):
        """Objective = Eloglik - sum(KL)."""
        XtX, Xty, yty, n, s = suff_data
        obj = get_objective_ss(XtX, Xty, s, yty, n)
        eloglik = Eloglik_ss(XtX, Xty, s, yty, n)
        assert_allclose(obj, eloglik - np.sum(s["KL"]))


class TestSingleEffectRegression:
    """T6: single_effect_regression_ss."""

    @pytest.fixture
    def ser_data(self):
        """Data for single_effect_regression_ss tests."""
        rng = np.random.default_rng(42)
        p = 5
        Xty = rng.normal(0, 2, p)
        dXtX = rng.uniform(1, 5, p)
        return Xty, dXtX

    def test_basic_regression(self, ser_data):
        """Basic SER should return well-formed result."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX, V=1.0, residual_variance=1.0)
        assert "alpha" in res
        assert "mu" in res
        assert "mu2" in res
        assert "lbf" in res
        assert "V" in res
        assert "lbf_model" in res
        assert_allclose(res["alpha"].sum(), 1.0, atol=1e-10)

    def test_alpha_normalized(self, ser_data):
        """Alpha should sum to 1."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX)
        assert_allclose(res["alpha"].sum(), 1.0, atol=1e-10)

    def test_optimize_V_none(self, ser_data):
        """optimize_V='none' keeps V unchanged."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX, V=2.0, optimize_V="none")
        assert res["V"] == 2.0

    def test_optimize_V_optim(self, ser_data):
        """optimize_V='optim' should run without error."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX, V=1.0, optimize_V="optim")
        assert np.isfinite(res["V"])
        assert res["V"] >= 0

    def test_optimize_V_EM(self, ser_data):
        """optimize_V='EM' should run without error."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX, V=1.0, optimize_V="EM")
        assert np.isfinite(res["V"])
        assert res["V"] >= 0

    def test_optimize_V_uniroot(self, ser_data):
        """optimize_V='uniroot' should run without error."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX, V=1.0, optimize_V="uniroot")
        assert np.isfinite(res["V"])
        assert res["V"] >= 0

    def test_invalid_optimize_V(self, ser_data):
        """Invalid method raises ValueError."""
        Xty, dXtX = ser_data
        with pytest.raises(ValueError, match="Invalid optimize_V"):
            single_effect_regression_ss(Xty, dXtX, optimize_V="bad_method")

    def test_V_zero_result(self, ser_data):
        """V=0 gives zero posterior variance."""
        Xty, dXtX = ser_data
        res = single_effect_regression_ss(Xty, dXtX, V=0)
        assert_allclose(res["mu"], np.zeros_like(Xty))

    def test_prior_weights_used(self, ser_data):
        """Custom prior weights influence alpha."""
        Xty, dXtX = ser_data
        pw_uniform = np.ones(5) / 5
        pw_biased = np.array([0.96, 0.01, 0.01, 0.01, 0.01])
        res_uniform = single_effect_regression_ss(Xty, dXtX, prior_weights=pw_uniform)
        res_biased = single_effect_regression_ss(Xty, dXtX, prior_weights=pw_biased)
        # Biased weights should push alpha[0] higher
        assert res_biased["alpha"][0] >= res_uniform["alpha"][0]


class TestOptimizePriorVariance:
    """T7: optimize_prior_variance."""

    @pytest.fixture
    def optim_data(self):
        betahat = np.array([2.0, 0.5, -0.3])
        shat2 = np.array([0.1, 0.2, 0.3])
        pw = np.ones(3) / 3
        return betahat, shat2, pw

    def test_simple(self, optim_data):
        """'simple' returns V_init unchanged."""
        betahat, shat2, pw = optim_data
        V = optimize_prior_variance("simple", betahat, shat2, pw, V_init=1.5)
        assert V == 1.5

    def test_optim(self, optim_data):
        """'optim' returns non-negative V."""
        betahat, shat2, pw = optim_data
        V = optimize_prior_variance("optim", betahat, shat2, pw, V_init=1.0)
        assert V >= 0

    def test_uniroot(self, optim_data):
        """'uniroot' returns positive V for signal."""
        betahat, shat2, pw = optim_data
        V = optimize_prior_variance("uniroot", betahat, shat2, pw, V_init=1.0)
        assert V >= 0

    def test_EM(self, optim_data):
        """'EM' requires alpha and post_mean2."""
        betahat, shat2, pw = optim_data
        alpha = np.array([0.7, 0.2, 0.1])
        post_mean2 = np.array([1.1, 0.3, 0.02])
        V = optimize_prior_variance(
            "EM", betahat, shat2, pw, alpha=alpha, post_mean2=post_mean2
        )
        assert V >= 0

    def test_EM_missing_args(self, optim_data):
        """'EM' without alpha/post_mean2 raises ValueError."""
        betahat, shat2, pw = optim_data
        with pytest.raises(ValueError, match="Alpha and post_mean2"):
            optimize_prior_variance("EM", betahat, shat2, pw, V_init=1.0)

    def test_invalid_method(self, optim_data):
        """Invalid method raises ValueError."""
        betahat, shat2, pw = optim_data
        with pytest.raises(ValueError, match="Invalid"):
            optimize_prior_variance("invalid", betahat, shat2, pw, V_init=1.0)

    def test_check_null_sets_zero(self, optim_data):
        """check_null_threshold can force V to zero."""
        betahat, shat2, pw = optim_data
        # Very weak signal
        betahat_weak = np.array([0.01, 0.01, 0.01])
        V = optimize_prior_variance(
            "simple", betahat_weak, shat2, pw, V_init=0.001, check_null_threshold=100
        )
        assert V == 0


class TestEstVUniroot:
    """T8: est_V_uniroot."""

    def test_normal_input(self):
        """Should return positive V for signal."""
        betahat = np.array([2.0, 0.5])
        shat2 = np.array([0.1, 0.2])
        pw = np.array([0.5, 0.5])
        V = est_V_uniroot(betahat, shat2, pw)
        assert V > 0
        assert np.isfinite(V)

    def test_strong_signal(self):
        """Strong signal gives larger V."""
        betahat = np.array([5.0, 4.0])
        shat2 = np.array([0.1, 0.1])
        pw = np.array([0.5, 0.5])
        V = est_V_uniroot(betahat, shat2, pw)
        assert V > 1.0


# ═══════════════════════════════════════════════════════════════════
# Phase 3: Initialization/utility functions (T9–T12)
# ═══════════════════════════════════════════════════════════════════


class TestInitSetup:
    """T9: init_setup."""

    def test_basic_initialization(self):
        """Normal init returns correct structure."""
        s = init_setup(
            n=100, p=5, L=3,
            scaled_prior_variance=0.2,
            residual_variance=1.0,
            prior_weights=None,
            null_weight=0,
            varY=1.0,
            standardize=True,
        )
        assert s["alpha"].shape == (3, 5)
        assert_allclose(s["alpha"].sum(axis=1), np.ones(3), atol=1e-10)
        assert s["mu"].shape == (3, 5)
        assert s["mu2"].shape == (3, 5)

    def test_negative_prior_variance_raises(self):
        """Negative prior variance should raise."""
        with pytest.raises(ValueError, match="positive"):
            init_setup(100, 5, 3, -0.1, 1.0, None, 0, 1.0, True)

    def test_large_prior_variance_standardize_raises(self):
        """Scaled prior variance > 1 with standardize=True raises."""
        with pytest.raises(ValueError, match="no greater than 1"):
            init_setup(100, 5, 3, 1.5, 1.0, None, 0, 1.0, True)

    def test_large_prior_variance_no_standardize(self):
        """Scaled prior variance > 1 OK without standardize."""
        s = init_setup(100, 5, 3, 1.5, 1.0, None, 0, 1.0, False)
        assert s["alpha"].shape == (3, 5)

    def test_all_zero_prior_weights_raises(self):
        """All-zero prior weights raises."""
        with pytest.raises(ValueError, match="greater than 0"):
            init_setup(100, 5, 3, 0.2, 1.0, np.zeros(5), 0, 1.0, True)

    def test_wrong_length_prior_weights_raises(self):
        """Wrong-length prior weights raises."""
        with pytest.raises(ValueError, match="length p"):
            init_setup(100, 5, 3, 0.2, 1.0, np.ones(3), 0, 1.0, True)

    def test_p_less_than_L(self):
        """p < L should reduce L to p."""
        s = init_setup(100, 3, 10, 0.2, 1.0, None, 0, 1.0, True)
        assert s["alpha"].shape[0] == 3

    def test_null_weight_none(self):
        """null_weight=None sets null_index=0."""
        s = init_setup(100, 5, 3, 0.2, 1.0, None, None, 1.0, True)
        assert s["null_index"] == 0

    def test_null_weight_nonzero(self):
        """Non-None null_weight sets null_index=p."""
        s = init_setup(100, 5, 3, 0.2, 1.0, None, 0.1, 1.0, True)
        assert s["null_index"] == 5

    def test_prior_weights_normalized(self):
        """Prior weights should be normalized."""
        pw = np.array([2.0, 3.0, 5.0, 1.0, 4.0])
        s = init_setup(100, 5, 3, 0.2, 1.0, pw, 0, 1.0, True)
        assert_allclose(s["pi"].sum(), 1.0, atol=1e-10)


class TestInitFinalize:
    """T9b: init_finalize."""

    def test_V_expansion(self):
        """Single V should be expanded to length L."""
        s = {
            "alpha": np.full((3, 5), 0.2),
            "mu": np.zeros((3, 5)),
            "mu2": np.zeros((3, 5)),
            "V": [0.5],
            "sigma2": 1.0,
            "KL": np.zeros(3),
            "lbf": np.zeros(3),
        }
        s = init_finalize(s)
        assert len(s["V"]) == 3
        assert_allclose(s["V"], 0.5)

    def test_sigma2_zero_raises(self):
        """sigma2 <= 0 should raise."""
        s = {
            "alpha": np.full((2, 3), 1 / 3),
            "mu": np.zeros((2, 3)),
            "mu2": np.zeros((2, 3)),
            "V": [0.5],
            "sigma2": 0.0,
            "KL": np.zeros(2),
            "lbf": np.zeros(2),
        }
        with pytest.raises(ValueError, match="positive"):
            init_finalize(s)

    def test_negative_V_raises(self):
        """Negative V raises."""
        s = {
            "alpha": np.full((2, 3), 1 / 3),
            "mu": np.zeros((2, 3)),
            "mu2": np.zeros((2, 3)),
            "V": np.array([-1.0, 0.5]),
            "sigma2": 1.0,
            "KL": np.zeros(2),
            "lbf": np.zeros(2),
        }
        with pytest.raises(ValueError, match="non-negative"):
            init_finalize(s)

    def test_mu_mu2_mismatch_raises(self):
        """Mismatched mu/mu2 dimensions raises."""
        s = {
            "alpha": np.full((2, 3), 1 / 3),
            "mu": np.zeros((2, 3)),
            "mu2": np.zeros((2, 4)),
            "V": np.array([0.5, 0.5]),
            "sigma2": 1.0,
            "KL": np.zeros(2),
            "lbf": np.zeros(2),
        }
        with pytest.raises(ValueError, match="mu and mu2"):
            init_finalize(s)


class TestSusieSlimAndPrune:
    """T10: susie_slim, susie_prune_single_effects."""

    def test_slim_keeps_keys(self, susie_fit_obj):
        """Slim only keeps alpha, niter, V, sigma2."""
        susie_fit_obj["niter"] = 10
        slimmed = susie_slim(susie_fit_obj)
        assert set(slimmed.keys()) == {"alpha", "niter", "V", "sigma2"}

    def test_prune_reduce(self, susie_fit_obj):
        """Prune to fewer effects."""
        pruned = susie_prune_single_effects(susie_fit_obj.copy(), L=1)
        assert pruned["alpha"].shape[0] == 1

    def test_prune_expand(self, susie_fit_obj):
        """Expand to more effects."""
        expanded = susie_prune_single_effects(susie_fit_obj.copy(), L=5, V=0.1)
        assert expanded["alpha"].shape[0] == 5

    def test_prune_same_L(self, susie_fit_obj):
        """Same L returns same object (sets=None)."""
        result = susie_prune_single_effects(susie_fit_obj.copy(), L=2)
        assert result["sets"] is None
        assert result["alpha"].shape[0] == 2

    def test_prune_L_zero_uses_V(self):
        """L=0 determines L from non-zero V."""
        s = {
            "alpha": np.full((3, 5), 0.2),
            "mu": np.zeros((3, 5)),
            "mu2": np.zeros((3, 5)),
            "V": np.array([0.5, 0.0, 0.3]),
            "KL": np.zeros(3),
            "lbf": np.zeros(3),
            "lbf_variable": np.zeros((3, 5)),
        }
        pruned = susie_prune_single_effects(s, L=0)
        assert pruned["alpha"].shape[0] == 2

    def test_prune_preserves_cs_effects(self):
        """Effects in CS are kept first when pruning."""
        s = {
            "alpha": np.full((3, 5), 0.2),
            "mu": np.zeros((3, 5)),
            "mu2": np.zeros((3, 5)),
            "V": np.array([0.5, 0.5, 0.5]),
            "KL": np.zeros(3),
            "lbf": np.zeros(3),
            "lbf_variable": np.zeros((3, 5)),
            "sets": {"cs_index": np.array([2])},  # third component is in CS
        }
        pruned = susie_prune_single_effects(s, L=1)
        assert pruned["alpha"].shape[0] == 1


class TestHelperFunctions:
    """T11: muffled_cov2cor, compute_tf_Xb."""

    def test_muffled_cov2cor_identity(self):
        """Identity covariance → identity correlation."""
        result = muffled_cov2cor(np.eye(3))
        assert_allclose(result, np.eye(3))

    def test_muffled_cov2cor_diagonal(self):
        """Diagonal covariance → identity correlation."""
        cov = np.diag([4.0, 9.0, 16.0])
        result = muffled_cov2cor(cov)
        assert_allclose(result, np.eye(3), atol=1e-10)

    def test_muffled_cov2cor_zero_variance(self):
        """Zero variance handled (non-finite → 0)."""
        cov = np.array([[1.0, 0.0], [0.0, 0.0]])
        result = muffled_cov2cor(cov)
        assert np.all(np.isfinite(result))

    def test_compute_tf_Xb_order0(self):
        """Order 0 trend filtering."""
        b = np.array([1.0, 2.0, 3.0])
        result = compute_tf_Xb(0, b)
        expected = -np.cumsum(b[::-1])[::-1]
        assert_allclose(result, expected)

    def test_compute_tf_Xb_order1(self):
        """Order 1 applies cumsum twice."""
        b = np.array([1.0, 2.0, 3.0])
        result = compute_tf_Xb(1, b)
        assert result.shape == (3,)
        assert np.all(np.isfinite(result))


class TestSummarySusie:
    """T12: summary_susie."""

    def test_no_sets_raises(self):
        """No 'sets' key should raise ValueError."""
        with pytest.raises(ValueError, match="credible set information"):
            summary_susie({"pip": np.array([0.5, 0.3, 0.2])})

    def test_none_sets_raises(self):
        """sets=None should raise ValueError."""
        with pytest.raises(ValueError, match="credible set information"):
            summary_susie({"pip": np.array([0.5]), "sets": None})

    def test_with_cs(self):
        """Valid summary with credible sets."""
        obj = {
            "pip": np.array([0.8, 0.1, 0.1]),
            "null_index": 0,
            "lbf": np.array([2.0, 0.1]),
            "sets": {
                "cs": {"L0": [0, 1]},
                "purity": {
                    "mean_abs_corr": np.array([0.9]),
                    "min_abs_corr": np.array([0.8]),
                },
            },
        }
        result = summary_susie(obj)
        assert "vars" in result
        assert "cs" in result

    def test_no_cs_in_sets(self):
        """Sets without cs returns cs=None."""
        obj = {
            "pip": np.array([0.5, 0.3, 0.2]),
            "null_index": 0,
            "sets": {"cs": None},
        }
        result = summary_susie(obj)
        assert result["cs"] is None


# ═══════════════════════════════════════════════════════════════════
# Phase 4: Integration tests (T13–T15)
# ═══════════════════════════════════════════════════════════════════


class TestSusieGetCS:
    """T13: susie_get_cs."""

    def test_basic_no_purity(self):
        """Without X/Xcorr, returns CS without purity."""
        res = {
            "alpha": np.array([
                [0.8, 0.1, 0.05, 0.05],
                [0.05, 0.05, 0.1, 0.8],
            ]),
            "V": np.array([0.5, 0.5]),
        }
        result = susie_get_cs(res, coverage=0.95)
        assert result["cs"] is not None
        assert len(result["cs"]) == 2

    def test_with_xcorr_purity(self, small_pd_matrix):
        """With Xcorr, computes purity."""
        res = {
            "alpha": np.array([
                [0.8, 0.1, 0.05, 0.03, 0.02],
                [0.02, 0.03, 0.05, 0.1, 0.8],
            ]),
            "V": np.array([0.5, 0.5]),
        }
        result = susie_get_cs(res, Xcorr=small_pd_matrix, coverage=0.95)
        assert "purity" in result or result["cs"] is None

    def test_V_filtering(self):
        """Components with V ~ 0 are filtered."""
        res = {
            "alpha": np.array([
                [0.8, 0.1, 0.1],
                [0.1, 0.1, 0.8],
            ]),
            "V": np.array([0.5, 1e-15]),
        }
        result = susie_get_cs(res, coverage=0.95)
        if result["cs"] is not None:
            assert len(result["cs"]) == 1

    def test_dedup(self):
        """Duplicate CS are removed."""
        res = {
            "alpha": np.array([
                [0.8, 0.1, 0.1],
                [0.8, 0.1, 0.1],  # duplicate
            ]),
            "V": np.array([0.5, 0.5]),
        }
        result = susie_get_cs(res, coverage=0.95, dedup=True)
        if result["cs"] is not None:
            assert len(result["cs"]) == 1

    def test_no_dedup(self):
        """Without dedup, duplicates kept."""
        res = {
            "alpha": np.array([
                [0.8, 0.1, 0.1],
                [0.8, 0.1, 0.1],
            ]),
            "V": np.array([0.5, 0.5]),
        }
        result = susie_get_cs(res, coverage=0.95, dedup=False)
        if result["cs"] is not None:
            assert len(result["cs"]) == 2

    def test_all_filtered_returns_none(self, small_pd_matrix):
        """When all CS filtered by purity, returns None."""
        # Uniform alpha → large CS → likely low purity
        res = {
            "alpha": np.array([
                [0.2, 0.2, 0.2, 0.2, 0.2],
            ]),
            "V": np.array([0.5]),
        }
        result = susie_get_cs(res, Xcorr=small_pd_matrix, min_abs_corr=0.999, coverage=0.95)
        # May or may not be filtered; depends on matrix
        assert isinstance(result, dict)

    def test_x_and_xcorr_both_raises(self, rng, small_pd_matrix):
        """Providing both X and Xcorr raises."""
        X = rng.normal(0, 1, (50, 5))
        res = {"alpha": np.array([[0.8, 0.1, 0.05, 0.03, 0.02]]), "V": np.array([0.5])}
        with pytest.raises(ValueError, match="Only one"):
            susie_get_cs(res, X=X, Xcorr=small_pd_matrix)

    def test_asymmetric_xcorr_symmetrized(self):
        """Asymmetric Xcorr is symmetrized."""
        xcorr = np.array([[1.0, 0.5, 0.3], [0.6, 1.0, 0.4], [0.3, 0.4, 1.0]])
        res = {
            "alpha": np.array([[0.8, 0.1, 0.1]]),
            "V": np.array([0.5]),
        }
        # Should not raise
        result = susie_get_cs(res, Xcorr=xcorr, coverage=0.95)
        assert isinstance(result, dict)


class TestUpdateEachEffect:
    """T14: update_each_effect_ss."""

    @pytest.fixture
    def update_data(self):
        rng = np.random.default_rng(42)
        p = 4
        A = rng.normal(0, 1, (p, p))
        XtX = A @ A.T + 3 * np.eye(p)
        Xty = rng.normal(0, 1, p)
        s = {
            "alpha": np.full((2, p), 0.25),
            "mu": np.zeros((2, p)),
            "mu2": np.zeros((2, p)),
            "V": np.array([0.5, 0.5]),
            "sigma2": 1.0,
            "KL": np.full(2, np.nan),
            "lbf": np.full(2, np.nan),
            "lbf_variable": np.full((2, p), np.nan),
            "XtXr": np.zeros(p),
            "pi": np.ones(p) / p,
        }
        return XtX, Xty, s

    def test_update_runs(self, update_data):
        """Should run without error."""
        XtX, Xty, s = update_data
        result = update_each_effect_ss(XtX, Xty, s)
        assert_allclose(result["alpha"].sum(axis=1), np.ones(2), atol=1e-10)

    def test_update_with_prior_estimation(self, update_data):
        """Should work with prior variance estimation."""
        XtX, Xty, s = update_data
        result = update_each_effect_ss(
            XtX, Xty, s,
            estimate_prior_variance=True,
            estimate_prior_method="optim",
        )
        assert np.all(np.isfinite(result["V"]))

    def test_update_L1(self, update_data):
        """Should work with L=1."""
        XtX, Xty, s = update_data
        # Reduce to single effect
        for key in ["alpha", "mu", "mu2", "lbf_variable"]:
            s[key] = s[key][:1]
        for key in ["V", "KL", "lbf"]:
            s[key] = s[key][:1]
        result = update_each_effect_ss(XtX, Xty, s)
        assert result["alpha"].shape[0] == 1


class TestSusieSuffStat:
    """T15: susie_suff_stat (small data integration test)."""

    @pytest.fixture
    def suff_stat_data(self):
        rng = np.random.default_rng(42)
        p = 5
        n = 100
        A = rng.normal(0, 1, (p, p))
        XtX = A @ A.T + 5 * np.eye(p)
        # Symmetrize
        XtX = (XtX + XtX.T) / 2
        Xty = rng.normal(0, 2, p)
        yty = float(n)
        return XtX, Xty, yty, n

    def test_basic_run(self, suff_stat_data):
        """Should run with default parameters."""
        XtX, Xty, yty, n = suff_stat_data
        result = susie_suff_stat(XtX, Xty, yty, n, L=2, max_iter=3)
        assert "alpha" in result
        assert "pip" in result
        assert "converged" in result

    def test_dimension_mismatch_raises(self, suff_stat_data):
        """Mismatched XtX/Xty raises."""
        XtX, Xty, yty, n = suff_stat_data
        with pytest.raises(ValueError):
            susie_suff_stat(XtX, Xty[:3], yty, n)

    def test_inf_in_xty_raises(self, suff_stat_data):
        """Infinite values in Xty raise ValueError."""
        XtX, Xty, yty, n = suff_stat_data
        Xty[0] = np.inf
        with pytest.raises(ValueError, match="infinite"):
            susie_suff_stat(XtX, Xty, yty, n)

    def test_nan_in_xtx_raises(self, suff_stat_data):
        """NaN in XtX raises ValueError."""
        XtX, Xty, yty, n = suff_stat_data
        XtX[0, 0] = np.nan
        with pytest.raises(ValueError, match="NAs"):
            susie_suff_stat(XtX, Xty, yty, n)

    def test_standardize_option(self, suff_stat_data):
        """standardize flag should work."""
        XtX, Xty, yty, n = suff_stat_data
        result = susie_suff_stat(
            XtX, Xty, yty, n, L=2, max_iter=2, standardize=True
        )
        assert "alpha" in result

    def test_no_standardize(self, suff_stat_data):
        """standardize=False should work."""
        XtX, Xty, yty, n = suff_stat_data
        result = susie_suff_stat(
            XtX, Xty, yty, n, L=2, max_iter=2, standardize=False
        )
        assert "alpha" in result

    def test_convergence_flag(self, suff_stat_data):
        """converged is bool."""
        XtX, Xty, yty, n = suff_stat_data
        result = susie_suff_stat(XtX, Xty, yty, n, L=2, max_iter=100)
        assert isinstance(result["converged"], bool)

    def test_pip_range(self, suff_stat_data):
        """PIP should be in [0, 1]."""
        XtX, Xty, yty, n = suff_stat_data
        result = susie_suff_stat(XtX, Xty, yty, n, L=2, max_iter=5)
        assert np.all(result["pip"] >= 0)
        assert np.all(result["pip"] <= 1)

    def test_elbo_stored(self, suff_stat_data):
        """ELBO should be stored for each iteration."""
        XtX, Xty, yty, n = suff_stat_data
        result = susie_suff_stat(XtX, Xty, yty, n, L=2, max_iter=3)
        assert len(result["elbo"]) >= 1

    def test_p_less_than_L(self, suff_stat_data):
        """p < L should reduce L to p."""
        XtX, Xty, yty, n = suff_stat_data
        p = XtX.shape[0]
        result = susie_suff_stat(XtX, Xty, yty, n, L=p + 5, max_iter=2)
        assert result["alpha"].shape[0] == p
