"""Unit tests for multisusie_rss.py internal functions.

Tests classes, pure utility functions, statistical computations,
credible set logic, and integration for MultiSuSiE.
"""

import functools

import numpy as np
import pytest
from numpy.testing import assert_allclose

from credtools.wrappers.multisusie_rss import (
    S,
    SER_RESULTS,
    Eloglik,
    SER_posterior_e_loglik,
    compute_lbf,
    compute_lbf_and_moments,
    compute_lbf_and_moments_safe,
    compute_lbf_no_moments,
    estimate_residual_variance_func,
    get_ER2,
    get_objective,
    get_purity_x,
    in_CS,
    in_CS_x,
    multisusie_rss,
    n_in_CS,
    n_in_CS_x,
    recover_R_from_XTX,
    recover_XTX_and_XTY,
    recover_XTX_and_XTY_from_Z,
    susie_get_cs,
    susie_get_pip,
    susie_multi_ss,
    update_ER2,
)


# ─── Fixtures ───────────────────────────────────────────────────────


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def small_K2_data(rng):
    """Small test data: K=2 populations, p=4 variants, L=2 effects."""
    K, p, L = 2, 4, 2
    pop_sizes = [1000, 800]

    XTX_list = []
    XTY_list = []
    R_list = []
    YTY_list = []
    for k in range(K):
        A = rng.normal(0, 1, (p, p))
        XtX = A @ A.T + 5 * np.eye(p)
        XtX = (XtX + XtX.T) / 2
        # Convert to correlation-like
        d = np.sqrt(np.diag(XtX))
        R = XtX / np.outer(d, d)
        np.fill_diagonal(R, 1.0)
        # Scale like real XTX
        XtX = R * (pop_sizes[k] - 1)
        XTX_list.append(XtX.astype(np.float32))
        XTY_list.append(rng.normal(0, 2, p).astype(np.float32))
        R_list.append(R.astype(np.float32))
        YTY_list.append(np.float32(pop_sizes[k] - 1))

    rho = np.array([[1.0, 0.8], [0.8, 1.0]], dtype=np.float32)

    return {
        "K": K, "p": p, "L": L,
        "pop_sizes": pop_sizes,
        "XTX_list": XTX_list,
        "XTY_list": XTY_list,
        "R_list": R_list,
        "YTY_list": YTY_list,
        "rho": rho,
    }


@pytest.fixture
def small_S_obj(small_K2_data):
    """Create a small S object for testing."""
    d = small_K2_data
    s = S(
        pop_sizes=d["pop_sizes"],
        L=d["L"],
        XTX_list=d["XTX_list"],
        scaled_prior_variance=0.2,
        residual_variance=np.ones(d["K"], dtype=np.float32),
        varY=np.ones(d["K"], dtype=np.float32),
        prior_weights=np.ones(d["p"], dtype=np.float32) / d["p"],
    )
    return s


@pytest.fixture
def X_l2_arr(small_K2_data):
    """Diagonal of XTX for each population, shape (K, p)."""
    d = small_K2_data
    return np.array(
        [np.diag(XTX) for XTX in d["XTX_list"]], dtype=np.float32
    )


# ═══════════════════════════════════════════════════════════════════
# Phase 1: Pure function tests (M1–M4)
# ═══════════════════════════════════════════════════════════════════


class TestSClass:
    """M1: S class initialization and shape validation."""

    def test_alpha_shape(self, small_S_obj, small_K2_data):
        """alpha should be L x p."""
        d = small_K2_data
        assert small_S_obj.alpha.shape == (d["L"], d["p"])

    def test_mu_shape(self, small_S_obj, small_K2_data):
        """mu should be K x L x p."""
        d = small_K2_data
        assert small_S_obj.mu.shape == (d["K"], d["L"], d["p"])

    def test_mu2_shape(self, small_S_obj, small_K2_data):
        """mu2 should be K x K x L x p."""
        d = small_K2_data
        assert small_S_obj.mu2.shape == (d["K"], d["K"], d["L"], d["p"])

    def test_alpha_uniform(self, small_S_obj, small_K2_data):
        """alpha initialized to 1/p."""
        d = small_K2_data
        assert_allclose(small_S_obj.alpha, 1.0 / d["p"], atol=1e-6)

    def test_converged_initial(self, small_S_obj):
        """converged starts as False."""
        assert small_S_obj.converged is False

    def test_V_shape(self, small_S_obj, small_K2_data):
        """V should be K x L."""
        d = small_K2_data
        assert small_S_obj.V.shape == (d["K"], d["L"])

    def test_float_type(self, small_K2_data):
        """Custom float type should be respected."""
        d = small_K2_data
        s = S(
            pop_sizes=d["pop_sizes"], L=d["L"], XTX_list=d["XTX_list"],
            scaled_prior_variance=0.2,
            residual_variance=np.ones(d["K"]),
            varY=np.ones(d["K"]),
            prior_weights=np.ones(d["p"]) / d["p"],
            float_type=np.float64,
        )
        assert s.alpha.dtype == np.float64
        assert s.mu.dtype == np.float64

    def test_ER2_initialized_nan(self, small_S_obj):
        """ER2 initialized to NaN."""
        assert np.all(np.isnan(small_S_obj.ER2))

    def test_KL_initialized_nan(self, small_S_obj):
        """KL initialized to NaN."""
        assert np.all(np.isnan(small_S_obj.KL))


class TestSERResultsClass:
    """M1b: SER_RESULTS class."""

    def test_attributes_stored(self):
        """All attributes should be stored."""
        alpha = np.array([0.5, 0.3, 0.2])
        mu = np.array([[1.0, 0.5, 0.1], [0.5, 1.0, 0.1]])
        mu2 = np.array([[[1.1, 0.3, 0.02], [0.5, 1.0, 0.1]],
                         [[0.5, 1.0, 0.1], [1.1, 0.3, 0.02]]])
        lbf = np.array([0.5, 0.3, 0.1])
        lbf_model = 1.0
        V = np.array([0.5, 0.3])
        ser = SER_RESULTS(alpha, mu, mu2, lbf, lbf_model, V)
        assert_allclose(ser.alpha, alpha)
        assert_allclose(ser.mu, mu)
        assert_allclose(ser.mu2, mu2)
        assert_allclose(ser.lbf, lbf)
        assert ser.lbf_model == lbf_model
        assert_allclose(ser.V, V)


class TestMultisusieCSFunctions:
    """M2: in_CS_x, in_CS, n_in_CS_x (multisusie versions)."""

    def test_n_in_CS_x_concentrated(self):
        """Concentrated PIP needs 1 variable."""
        alpha = np.array([0.95, 0.02, 0.02, 0.01])
        assert n_in_CS_x(alpha, 0.9) == 1

    def test_n_in_CS_x_uniform(self):
        """Uniform PIP needs many variables."""
        alpha = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        assert n_in_CS_x(alpha, 0.9) == 5

    def test_in_CS_x_binary(self):
        """Output should be binary 0/1."""
        alpha = np.array([0.6, 0.2, 0.1, 0.1])
        result = in_CS_x(alpha, 0.8)
        assert set(np.unique(result)).issubset({0, 1})

    def test_in_CS_x_coverage(self):
        """Included variables achieve coverage."""
        alpha = np.array([0.6, 0.2, 0.1, 0.1])
        result = in_CS_x(alpha, 0.8)
        assert alpha[result == 1].sum() >= 0.8

    def test_in_CS_matrix(self):
        """in_CS works on LxP matrix."""
        alpha = np.array([
            [0.8, 0.1, 0.1],
            [0.1, 0.1, 0.8],
        ])
        result = in_CS(alpha, 0.9)
        assert result.shape == (2, 3)
        assert result[0, 0] == 1
        assert result[1, 2] == 1

    def test_n_in_CS_array(self):
        """n_in_CS returns per-component counts."""
        alpha = np.array([
            [0.9, 0.05, 0.05],
            [0.33, 0.34, 0.33],
        ])
        result = n_in_CS(alpha, 0.9)
        assert result[0] == 1
        assert result[1] == 3


class TestRecoverFunctions:
    """M3: recover_XTX_and_XTY, recover_XTX_and_XTY_from_Z."""

    def test_recover_from_bhat_shat(self, rng):
        """recover_XTX_and_XTY returns correct shapes."""
        p = 5
        b = rng.normal(0, 0.1, p)
        s = np.abs(rng.normal(0.02, 0.005, p)) + 1e-6
        R = np.eye(p, dtype=np.float32)
        YTY = 1.0
        n = 1000
        R_copy = R.copy()
        XTX, XTY = recover_XTX_and_XTY(b, s, R_copy, YTY, n)
        assert XTX.shape == (p, p)
        assert XTY.shape == (p,)
        # R is mutated in-place
        assert not np.allclose(R_copy, R)

    def test_recover_from_z(self, rng):
        """recover_XTX_and_XTY_from_Z returns correct shapes."""
        p = 5
        z = rng.normal(0, 2, p)
        R = np.eye(p, dtype=np.float32) * 999  # will be scaled
        R_copy = R.copy()
        n = 1000
        XTX, XTY = recover_XTX_and_XTY_from_Z(z, R_copy, n)
        assert XTX.shape == (p, p)
        assert XTY.shape == (p,)

    def test_recover_nan_handling(self):
        """NaN in input handled gracefully."""
        p = 3
        b = np.array([0.1, np.nan, 0.3])
        s = np.array([0.02, 0.02, 0.02])
        R = np.eye(p, dtype=np.float32)
        XTX, XTY = recover_XTX_and_XTY(b, s, R.copy(), 1.0, 1000)
        assert np.all(np.isfinite(XTX))
        assert np.all(np.isfinite(XTY))

    def test_recover_z_nan_handling(self):
        """NaN in z scores handled."""
        z = np.array([1.0, np.nan, 2.0])
        R = np.eye(3, dtype=np.float32)
        XTX, XTY = recover_XTX_and_XTY_from_Z(z, R.copy(), 1000)
        assert np.all(np.isfinite(XTY))

    def test_recover_z_preserves_original(self, rng):
        """Original R should not be used after mutation."""
        p = 3
        z = rng.normal(0, 2, p)
        R_orig = np.eye(p, dtype=np.float32)
        R_test = R_orig.copy()
        recover_XTX_and_XTY_from_Z(z, R_test, 1000)
        # R_test is mutated
        assert not np.allclose(R_test, R_orig)


class TestRecoverRFromXTX:
    """M4: recover_R_from_XTX."""

    def test_recover_identity(self):
        """XTX proportional to identity → R = identity."""
        p = 3
        XTX = np.eye(p, dtype=np.float64) * 5.0
        X_l2 = np.diag(XTX).copy()
        recover_R_from_XTX(XTX, X_l2)
        assert_allclose(XTX, np.eye(p), atol=1e-10)

    def test_recover_with_zeros(self):
        """Handles X_l2 containing zeros."""
        XTX = np.array([
            [5.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 3.0],
        ])
        X_l2 = np.diag(XTX).copy()
        recover_R_from_XTX(XTX, X_l2)
        assert_allclose(XTX[0, 0], 1.0, atol=1e-10)
        assert_allclose(XTX[2, 2], 1.0, atol=1e-10)


# ═══════════════════════════════════════════════════════════════════
# Phase 2: Statistical computation core (M5–M9)
# ═══════════════════════════════════════════════════════════════════


class TestGetER2Multi:
    """M5: get_ER2, SER_posterior_e_loglik, update_ER2, estimate_residual_variance_func."""

    def test_get_ER2_zero_effects(self, small_K2_data):
        """Zero effects → ER2 = YTY."""
        d = small_K2_data
        p = d["p"]
        alpha = np.zeros((2, p), dtype=np.float32)
        mu = np.zeros((2, p), dtype=np.float32)
        mu2 = np.zeros((2, p), dtype=np.float32)
        X_l2 = np.diag(d["XTX_list"][0])
        result = get_ER2(d["XTX_list"][0], d["XTY_list"][0], d["YTY_list"][0],
                         alpha, mu, mu2, X_l2)
        assert_allclose(result, d["YTY_list"][0], atol=1e-3)

    def test_get_ER2_finite(self, small_K2_data, small_S_obj, X_l2_arr):
        """Normal input produces finite ER2."""
        d = small_K2_data
        result = get_ER2(
            d["XTX_list"][0], d["XTY_list"][0], d["YTY_list"][0],
            small_S_obj.alpha, small_S_obj.mu[0], small_S_obj.mu2[0, 0],
            X_l2_arr[0],
        )
        assert np.isfinite(result)

    def test_update_ER2_shape(self, small_K2_data, small_S_obj, X_l2_arr):
        """update_ER2 returns K-length array."""
        d = small_K2_data
        result = update_ER2(
            d["XTX_list"], d["XTY_list"], d["YTY_list"], small_S_obj, X_l2_arr
        )
        assert result.shape == (d["K"],)
        assert np.all(np.isfinite(result))

    def test_SER_posterior_e_loglik_zero(self, small_K2_data, X_l2_arr):
        """Zero effects give zero loglik contribution."""
        d = small_K2_data
        K, p = d["K"], d["p"]
        Eb = [np.zeros(p) for _ in range(K)]
        Eb2 = [np.zeros(p) for _ in range(K)]
        s2 = np.ones(K)
        result = SER_posterior_e_loglik(X_l2_arr, d["XTY_list"], s2, Eb, Eb2)
        assert_allclose(result, 0.0, atol=1e-6)

    def test_SER_posterior_e_loglik_finite(self, small_K2_data, X_l2_arr, rng):
        """Normal input produces finite value."""
        d = small_K2_data
        K, p = d["K"], d["p"]
        Eb = [rng.normal(0, 0.1, p).astype(np.float32) for _ in range(K)]
        Eb2 = [np.abs(rng.normal(0, 0.1, p)).astype(np.float32) for _ in range(K)]
        s2 = np.ones(K, dtype=np.float32)
        result = SER_posterior_e_loglik(X_l2_arr, d["XTY_list"], s2, Eb, Eb2)
        assert np.isfinite(result)

    def test_estimate_residual_variance(self, small_K2_data, small_S_obj, X_l2_arr):
        """Residual variance should be positive."""
        d = small_K2_data
        # Set ER2 to something valid
        small_S_obj.ER2 = update_ER2(
            d["XTX_list"], d["XTY_list"], d["YTY_list"], small_S_obj, X_l2_arr
        )
        result = estimate_residual_variance_func(
            d["XTX_list"], d["XTY_list"], d["YTY_list"],
            small_S_obj, X_l2_arr, d["pop_sizes"],
        )
        assert result.shape == (d["K"],)
        assert np.all(np.isfinite(result))


class TestEloglikAndObjective:
    """M6: Eloglik, get_objective."""

    def test_Eloglik_finite(self, small_K2_data, small_S_obj, X_l2_arr):
        """Eloglik should be finite."""
        d = small_K2_data
        small_S_obj.ER2 = update_ER2(
            d["XTX_list"], d["XTY_list"], d["YTY_list"], small_S_obj, X_l2_arr
        )
        result = Eloglik(
            d["XTX_list"], d["XTY_list"], small_S_obj, d["YTY_list"], X_l2_arr
        )
        assert np.isfinite(result)

    def test_get_objective_relationship(self, small_K2_data, small_S_obj, X_l2_arr):
        """Objective = Eloglik - sum(KL)."""
        d = small_K2_data
        small_S_obj.ER2 = update_ER2(
            d["XTX_list"], d["XTY_list"], d["YTY_list"], small_S_obj, X_l2_arr
        )
        small_S_obj.KL = np.array([0.1, 0.2], dtype=np.float32)
        obj = get_objective(
            d["XTX_list"], d["XTY_list"], small_S_obj, d["YTY_list"], X_l2_arr
        )
        eloglik = Eloglik(
            d["XTX_list"], d["XTY_list"], small_S_obj, d["YTY_list"], X_l2_arr
        )
        assert_allclose(obj, eloglik - np.sum(small_S_obj.KL), atol=1e-5)


class TestSusieGetPIPMulti:
    """M7: susie_get_pip (multisusie version)."""

    def test_pip_range(self, small_S_obj):
        """PIP in [0, 1]."""
        pip = susie_get_pip(small_S_obj)
        assert np.all(pip >= 0)
        assert np.all(pip <= 1)

    def test_pip_V_filtering(self, small_K2_data):
        """Components with V near zero filtered."""
        d = small_K2_data
        s = S(
            pop_sizes=d["pop_sizes"], L=d["L"], XTX_list=d["XTX_list"],
            scaled_prior_variance=0.2,
            residual_variance=np.ones(d["K"], dtype=np.float32),
            varY=np.ones(d["K"], dtype=np.float32),
            prior_weights=np.ones(d["p"], dtype=np.float32) / d["p"],
        )
        # Zero out V for second effect
        s.V[:, 1] = 0
        pip = susie_get_pip(s)
        assert pip.shape == (d["p"],)

    def test_pip_all_zero_V(self, small_K2_data):
        """All zero V returns zero PIP."""
        d = small_K2_data
        s = S(
            pop_sizes=d["pop_sizes"], L=d["L"], XTX_list=d["XTX_list"],
            scaled_prior_variance=0.0,
            residual_variance=np.ones(d["K"], dtype=np.float32),
            varY=np.ones(d["K"], dtype=np.float32),
            prior_weights=np.ones(d["p"], dtype=np.float32) / d["p"],
        )
        pip = susie_get_pip(s)
        assert_allclose(pip, np.zeros(d["p"]), atol=1e-6)


class TestSusieGetCSMulti:
    """M8: susie_get_cs (multisusie version)."""

    def test_basic_cs(self, small_S_obj, small_K2_data):
        """Basic CS extraction."""
        d = small_K2_data
        cs, purity, coverage, mask = susie_get_cs(
            small_S_obj, d["R_list"], coverage=0.95, purity=0.0,
        )
        assert isinstance(cs, list)
        assert isinstance(mask, np.ndarray)

    def test_no_purity_calc(self, small_S_obj, small_K2_data):
        """calculate_purity=False skips purity."""
        d = small_K2_data
        cs, purity, coverage, mask = susie_get_cs(
            small_S_obj, d["R_list"], coverage=0.95,
            calculate_purity=False,
        )
        assert np.all(np.isnan(purity[mask]))

    def test_dedup(self, small_K2_data):
        """Duplicate CS are removed."""
        d = small_K2_data
        s = S(
            pop_sizes=d["pop_sizes"], L=2, XTX_list=d["XTX_list"],
            scaled_prior_variance=0.2,
            residual_variance=np.ones(d["K"], dtype=np.float32),
            varY=np.ones(d["K"], dtype=np.float32),
            prior_weights=np.ones(d["p"], dtype=np.float32) / d["p"],
        )
        # Make both effects identical → same CS
        s.alpha[1] = s.alpha[0]
        cs, purity, coverage, mask = susie_get_cs(
            s, d["R_list"], dedup=True, calculate_purity=False,
        )
        # At most 1 CS should be included due to dedup
        assert np.sum(mask) <= 1

    def test_zero_V_excluded(self, small_K2_data):
        """Effects with V near zero excluded."""
        d = small_K2_data
        s = S(
            pop_sizes=d["pop_sizes"], L=2, XTX_list=d["XTX_list"],
            scaled_prior_variance=0.2,
            residual_variance=np.ones(d["K"], dtype=np.float32),
            varY=np.ones(d["K"], dtype=np.float32),
            prior_weights=np.ones(d["p"], dtype=np.float32) / d["p"],
        )
        s.V[:, 1] = 0
        cs, purity, coverage, mask = susie_get_cs(
            s, d["R_list"], calculate_purity=False,
        )
        assert mask[1] is False or mask[1] == False  # noqa: E712


class TestGetPurityX:
    """M9: get_purity_x."""

    def test_single_variant(self, small_K2_data):
        """Single variant → purity 1 (minimum of identity matrix)."""
        d = small_K2_data
        result = get_purity_x(np.array([0]), d["R_list"], 0.5, 100, None)
        assert result == 1.0

    def test_two_variants(self, small_K2_data):
        """Two variants should have finite purity."""
        d = small_K2_data
        result = get_purity_x(np.array([0, 1]), d["R_list"], 0.5, 100, None)
        assert np.isfinite(result)
        assert result >= 0

    def test_subsampling(self, small_K2_data):
        """Large CS uses subsampling."""
        d = small_K2_data
        result = get_purity_x(np.array([0, 1, 2, 3]), d["R_list"], 0.5, 2, None)
        assert np.isfinite(result)

    def test_from_X_list(self, rng):
        """Compute purity from X matrices."""
        X_list = [rng.normal(0, 1, (50, 5)) for _ in range(2)]
        result = get_purity_x(np.array([0, 1, 2]), None, 0.5, 100, X_list)
        assert np.isfinite(result)


# ═══════════════════════════════════════════════════════════════════
# Phase 3: Compute LBF and optimization (M10–M11)
# ═══════════════════════════════════════════════════════════════════


class TestComputeLBF:
    """M10: compute_lbf, compute_lbf_no_moments, compute_lbf_and_moments, compute_lbf_and_moments_safe."""

    @pytest.fixture
    def lbf_data(self, small_K2_data, X_l2_arr):
        """Data for compute_lbf tests."""
        d = small_K2_data
        V = np.array([0.1, 0.1], dtype=np.float64)
        residual_variance = np.ones(d["K"], dtype=np.float64)
        XTY_list = [x.astype(np.float64) for x in d["XTY_list"]]
        XTX_list = [x.astype(np.float64) for x in d["XTX_list"]]
        X_l2 = X_l2_arr.astype(np.float64)
        return V, XTY_list, XTX_list, X_l2, d["rho"].astype(np.float64), residual_variance

    def test_compute_lbf_all_zero_V(self, lbf_data):
        """All-zero V returns zero LBF."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        V_zero = np.zeros(2, dtype=np.float64)
        lbf = compute_lbf(V_zero, XTY_list, XTX_list, X_l2, rho, resvar)
        assert_allclose(lbf, np.zeros(4), atol=1e-6)

    def test_compute_lbf_all_zero_V_with_moments(self, lbf_data):
        """All-zero V with moments returns zeros."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        V_zero = np.zeros(2, dtype=np.float64)
        lbf, pm, pm2 = compute_lbf(
            V_zero, XTY_list, XTX_list, X_l2, rho, resvar, return_moments=True,
        )
        assert_allclose(lbf, np.zeros(4), atol=1e-6)
        assert_allclose(pm, np.zeros((2, 4)), atol=1e-6)

    def test_compute_lbf_partial_zero_V(self, lbf_data):
        """One pop with zero V, recursive call."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        V_partial = np.array([0.1, 0.0], dtype=np.float64)
        lbf = compute_lbf(V_partial, XTY_list, XTX_list, X_l2, rho, resvar)
        assert lbf.shape == (4,)
        assert np.all(np.isfinite(lbf))

    def test_compute_lbf_nonzero_V(self, lbf_data):
        """Non-zero V gives finite LBF."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        lbf = compute_lbf(V, XTY_list, XTX_list, X_l2, rho, resvar)
        assert lbf.shape == (4,)
        assert np.all(np.isfinite(lbf))

    def test_compute_lbf_with_moments(self, lbf_data):
        """compute_lbf with return_moments=True."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        lbf, pm, pm2 = compute_lbf(
            V, XTY_list, XTX_list, X_l2, rho, resvar, return_moments=True,
        )
        assert lbf.shape == (4,)
        assert pm.shape == (2, 4)
        assert pm2.shape == (2, 2, 4)

    def test_no_moments_vs_moments_lbf_agree(self, lbf_data):
        """LBF should agree whether moments are computed or not."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        lbf_no = compute_lbf(V, XTY_list, XTX_list, X_l2, rho, resvar, return_moments=False)
        lbf_yes, _, _ = compute_lbf(V, XTY_list, XTX_list, X_l2, rho, resvar, return_moments=True)
        assert_allclose(lbf_no, lbf_yes, atol=1e-5, rtol=1e-4)

    def test_compute_lbf_and_moments_safe_consistent(self, lbf_data):
        """Safe version should produce similar results."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        XTY = np.ascontiguousarray(np.stack(XTY_list, axis=1))

        lbf_std, pm_std, pm2_std = compute_lbf_and_moments(
            V, XTY, X_l2, rho, resvar, np.float64,
        )
        lbf_safe, pm_safe, pm2_safe = compute_lbf_and_moments_safe(
            V, XTY, X_l2, rho, resvar, np.float64,
        )
        assert_allclose(lbf_std, lbf_safe, atol=1e-5, rtol=1e-4)
        assert_allclose(pm_std, pm_safe, atol=1e-5, rtol=1e-4)

    def test_compute_lbf_no_moments_direct(self, lbf_data):
        """Direct call to compute_lbf_no_moments."""
        V, XTY_list, XTX_list, X_l2, rho, resvar = lbf_data
        XTY = np.ascontiguousarray(np.stack(XTY_list, axis=1))
        lbf = compute_lbf_no_moments(V, XTY, X_l2, rho, resvar, np.float64)
        assert lbf.shape == (4,)
        assert np.all(np.isfinite(lbf))


class TestMultisusieOptimizePriorVariance:
    """M11: optimize_prior_variance (multisusie version), loglik (multisusie version)."""

    @pytest.fixture
    def optim_setup(self, small_K2_data, X_l2_arr):
        """Setup for optimization tests."""
        d = small_K2_data
        K, p = d["K"], d["p"]
        prior_weights = np.ones(p, dtype=np.float32) / p
        V = np.array([0.1, 0.1], dtype=np.float64)
        residual_variance = np.ones(K, dtype=np.float64)
        XTY_list = [x.astype(np.float64) for x in d["XTY_list"]]
        XTX_list = [x.astype(np.float64) for x in d["XTX_list"]]
        X_l2 = X_l2_arr.astype(np.float64)
        rho = d["rho"].astype(np.float64)
        compute_lbf_params = (XTY_list, XTX_list, X_l2, rho, residual_variance)
        # Mock alpha/post_mean2 for EM
        alpha = np.ones((1, p), dtype=np.float32) / p
        post_mean2 = np.zeros((K, K, p), dtype=np.float32)
        post_mean2[0, 0] = 0.01
        post_mean2[1, 1] = 0.01
        w_pop = np.array([0.5, 0.5], dtype=np.float32)
        return {
            "prior_weights": prior_weights, "K": K,
            "compute_lbf_params": compute_lbf_params,
            "alpha": alpha[0], "post_mean2": post_mean2,
            "w_pop": w_pop, "current_V": V,
        }

    def test_EM_method(self, optim_setup):
        """EM method returns positive V."""
        from credtools.wrappers.multisusie_rss import optimize_prior_variance as opt_pv
        s = optim_setup
        V = opt_pv(
            "EM", s["prior_weights"], s["K"],
            compute_lbf_params=s["compute_lbf_params"],
            alpha=s["alpha"], post_mean2=s["post_mean2"],
            w_pop=s["w_pop"], current_V=s["current_V"],
        )
        assert np.all(np.isfinite(np.atleast_1d(V)))

    def test_early_EM_method(self, optim_setup):
        """early_EM method works."""
        from credtools.wrappers.multisusie_rss import optimize_prior_variance as opt_pv
        s = optim_setup
        V = opt_pv(
            "early_EM", s["prior_weights"], s["K"],
            compute_lbf_params=s["compute_lbf_params"],
            alpha=s["alpha"], post_mean2=s["post_mean2"],
            w_pop=s["w_pop"], current_V=s["current_V"],
        )
        assert np.all(np.isfinite(np.atleast_1d(V)))

    def test_optim_method(self, optim_setup):
        """optim method works."""
        from credtools.wrappers.multisusie_rss import optimize_prior_variance as opt_pv
        s = optim_setup
        V = opt_pv(
            "optim", s["prior_weights"], s["K"],
            compute_lbf_params=s["compute_lbf_params"],
            alpha=s["alpha"], post_mean2=s["post_mean2"],
            w_pop=s["w_pop"], current_V=s["current_V"],
        )
        V_arr = np.atleast_1d(V)
        assert np.all(V_arr >= 0)

    def test_grid_method(self, optim_setup):
        """grid method works."""
        from credtools.wrappers.multisusie_rss import optimize_prior_variance as opt_pv
        s = optim_setup
        V = opt_pv(
            "grid", s["prior_weights"], s["K"],
            compute_lbf_params=s["compute_lbf_params"],
            alpha=s["alpha"], post_mean2=s["post_mean2"],
            w_pop=s["w_pop"], current_V=s["current_V"],
        )
        V_arr = np.atleast_1d(V)
        assert np.all(V_arr >= 0)

    def test_invalid_method(self, optim_setup):
        """Invalid method raises ValueError."""
        from credtools.wrappers.multisusie_rss import optimize_prior_variance as opt_pv
        s = optim_setup
        with pytest.raises(ValueError, match="unknown"):
            opt_pv(
                "invalid_method", s["prior_weights"], s["K"],
                compute_lbf_params=s["compute_lbf_params"],
                alpha=s["alpha"], post_mean2=s["post_mean2"],
                w_pop=s["w_pop"], current_V=s["current_V"],
            )


# ═══════════════════════════════════════════════════════════════════
# Phase 4: Integration tests (M12–M13)
# ═══════════════════════════════════════════════════════════════════


class TestMultisusieRSSInputValidation:
    """M12: multisusie_rss input validation."""

    def test_z_and_b_both_raises(self, small_K2_data):
        """Providing both z_list and b_list raises."""
        d = small_K2_data
        z_list = [np.zeros(d["p"]) for _ in range(d["K"])]
        b_list = [np.zeros(d["p"]) for _ in range(d["K"])]
        s_list = [np.ones(d["p"]) * 0.01 for _ in range(d["K"])]
        with pytest.raises(Exception):
            multisusie_rss(
                R_list=d["R_list"], population_sizes=d["pop_sizes"],
                z_list=z_list, b_list=b_list, s_list=s_list,
                rho=d["rho"],
            )

    def test_neither_z_nor_b_raises(self, small_K2_data):
        """Providing neither z_list nor b_list raises."""
        d = small_K2_data
        with pytest.raises(Exception):
            multisusie_rss(
                R_list=d["R_list"], population_sizes=d["pop_sizes"],
                rho=d["rho"],
            )


class TestSusieMultiSS:
    """M13: susie_multi_ss integration test."""

    def test_basic_run(self, small_K2_data):
        """Should run with small data and few iterations."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=3,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert isinstance(s, S)
        assert s.alpha.shape == (2, d["p"])

    def test_standardize_off(self, small_K2_data):
        """standardize=False should work."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=2,
            standardize=False,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert isinstance(s, S)

    def test_p_less_than_L(self, small_K2_data):
        """p < L clips L to p."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=100,
            max_iter=2,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert s.alpha.shape[0] == d["p"]

    def test_prior_weights_normalized(self, small_K2_data):
        """Custom prior weights are normalized."""
        d = small_K2_data
        pw = np.array([2.0, 3.0, 5.0, 1.0], dtype=np.float32)
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=2,
            prior_weights=pw,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert_allclose(s.pi.sum(), 1.0, atol=1e-5)

    def test_with_prior_variance_estimation(self, small_K2_data):
        """estimate_prior_variance=True should work."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=3,
            estimate_prior_variance=True,
            estimate_prior_method="early_EM",
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert isinstance(s, S)

    def test_converged_attribute(self, small_K2_data):
        """converged attribute should be bool."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=2,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert isinstance(s.converged, bool)

    def test_sets_attribute(self, small_K2_data):
        """Result should have sets attribute from susie_get_cs."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=3,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert hasattr(s, "sets")

    def test_pip_attribute(self, small_K2_data):
        """Result should have pip attribute."""
        d = small_K2_data
        s = susie_multi_ss(
            XTX_list=d["XTX_list"],
            XTY_list=d["XTY_list"],
            YTY_list=d["YTY_list"],
            rho=d["rho"],
            population_sizes=d["pop_sizes"],
            L=2,
            max_iter=3,
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            R_list=d["R_list"],
        )
        assert hasattr(s, "pip")
        assert np.all(s.pip >= 0)
        assert np.all(s.pip <= 1)
