"""Test wrapper of the `susie_rss` function."""

from typing import Any, Dict

import numpy as np
import pytest
from numpy.testing import assert_allclose

# Import the function to be tested
from credtools.wrappers.susie_rss import susie_rss


# Mock the susie_suff_stat function
def mock_susie_suff_stat(**kwargs):
    """Mock the susie_suff_stat function."""
    return {
        "alpha": np.array([[0.1, 0.2, 0.3]]),
        "mu": np.array([[1.0, 2.0, 3.0]]),
        "mu2": np.array([[1.1, 4.1, 9.1]]),
        "lbf": np.array([0.5]),
        "lbf_variable": np.array([[0.1, 0.2, 0.3]]),
        "V": 1.0,
        "elbo": np.array([10.0]),
        "sets": {"cs": [0, 1, 2]},
        "pip": np.array([0.1, 0.2, 0.3]),
        "niter": 10,
        "converged": True,
    }


# Patch the susie_suff_stat function
@pytest.fixture(autouse=True)
def patch_susie_suff_stat(monkeypatch):
    """Patch the susie_suff_stat function."""
    monkeypatch.setattr(
        "credtools.wrappers.susie_rss.susie_suff_stat", mock_susie_suff_stat
    )


def test_susie_rss_with_z():
    """Test the susie_rss function with the `z` parameter."""
    z = np.array([1.0, 2.0, 3.0])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100

    result = susie_rss(z=z, R=R, n=n)

    assert isinstance(result, dict)
    assert "alpha" in result
    assert "mu" in result
    assert "mu2" in result
    assert "lbf" in result
    assert "lbf_variable" in result
    assert "V" in result
    assert "elbo" in result
    assert "sets" in result
    assert "pip" in result
    assert "niter" in result
    assert "converged" in result


def test_susie_rss_with_bhat_shat():
    """Test the susie_rss function with beta and se."""
    bhat = np.array([0.1, 0.2, 0.3])
    shat = np.array([0.01, 0.02, 0.03])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100

    result = susie_rss(bhat=bhat, shat=shat, R=R, n=n)

    assert isinstance(result, dict)
    assert "alpha" in result


def test_susie_rss_with_var_y():
    """Test susie_rss with var_y."""
    z = np.array([1.0, 2.0, 3.0])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100
    var_y = 1.5

    result = susie_rss(z=z, R=R, n=n, var_y=var_y)

    assert isinstance(result, dict)
    assert "alpha" in result


def test_susie_rss_estimate_residual_variance():
    """Test susie_rss estimate_residual_variance."""
    z = np.array([1.0, 2.0, 3.0])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100

    result = susie_rss(z=z, R=R, n=n, estimate_residual_variance=True)

    assert isinstance(result, dict)
    assert "alpha" in result


def test_susie_rss_invalid_inputs():
    """Test susie_rss invalid inputes."""
    with pytest.raises(ValueError):
        susie_rss(z=np.array([1.0, 2.0]), R=np.array([[1.0, 0.5], [0.5, 1.0]]), n=1)

    with pytest.raises(ValueError):
        susie_rss(
            z=np.array([1.0, 2.0]),
            R=np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]]),
        )

    with pytest.raises(ValueError):
        susie_rss(
            bhat=np.array([1.0, 2.0]),
            shat=np.array([0.1]),
            R=np.array([[1.0, 0.5], [0.5, 1.0]]),
        )

    with pytest.raises(ValueError):
        susie_rss(
            bhat=np.array([1.0, 2.0]),
            shat=np.array([0.1, -0.1]),
            R=np.array([[1.0, 0.5], [0.5, 1.0]]),
        )


def test_susie_rss_z_ld_weight():
    """Test susie_rss z_ld_weight."""
    z = np.array([1.0, 2.0, 3.0])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100

    result = susie_rss(z=z, R=R, n=n, z_ld_weight=0.1)

    assert isinstance(result, dict)
    assert "alpha" in result


def test_susie_rss_check_prior():
    """Test susie_rss check_prior."""
    z = np.array([1.0, 2.0, 3.0])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100

    result = susie_rss(z=z, R=R, n=n, check_prior=False)

    assert isinstance(result, dict)
    assert "alpha" in result


def test_susie_rss_additional_kwargs():
    """Test susie_rss addition kwargs."""
    z = np.array([1.0, 2.0, 3.0])
    R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])
    n = 100

    result = susie_rss(z=z, R=R, n=n, L=5, tol=1e-5)

    assert isinstance(result, dict)
    assert "alpha" in result
