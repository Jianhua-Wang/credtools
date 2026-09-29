"""Regression tests for finite, nonnegative ELBO convergence in both engines."""

import importlib

import numpy as np
import pytest


@pytest.mark.parametrize("engine", ["susie", "multisusie"])
@pytest.mark.parametrize(
    "trace, converged, niter",
    [
        ([10.0, 0.0, 0.0005], True, 3),
        ([10.0, 9.9999], False, 2),
        ([10.0, 10.0], True, 2),
        ([10.0, 10.01], False, 2),
    ],
)
def test_decrease_is_not_convergence(monkeypatch, engine, trace, converged, niter):
    """A decreasing step must not end IBSS, even when its magnitude is tiny."""
    result = _run(monkeypatch, engine, trace)
    actual = result["converged"] if engine == "susie" else result.converged
    iterations = result["niter"] if engine == "susie" else result.niter
    assert actual is converged
    assert iterations == niter


@pytest.mark.parametrize("engine", ["susie", "multisusie"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_elbo_is_error(monkeypatch, engine, value):
    """Computed infinities/NaNs differ from the initial -Inf sentinel."""
    with pytest.raises(ValueError, match="nonfinite"):
        _run(monkeypatch, engine, [value])


def _run(monkeypatch, engine, trace):
    p = 4
    if engine == "susie":
        mod = importlib.import_module("credtools.wrappers.susie_rss")
        values = iter(np.repeat(trace, 2))  # Debug log and stored objective.
        monkeypatch.setattr(mod, "get_objective_ss", lambda *args: next(values))
        return mod.susie_suff_stat(
            np.eye(p) * 100,
            np.array([5.0, 1.0, 0.0, 0.0]),
            100.0,
            100,
            L=1,
            max_iter=len(trace),
            estimate_prior_variance=False,
            estimate_residual_variance=False,
            tol=0.001,
            coverage=None,  # Isolate stopping from credible-set purity handling.
        )
    mod = importlib.import_module("credtools.wrappers.multisusie_rss")
    values = iter(trace)
    monkeypatch.setattr(mod, "get_objective", lambda *args: next(values))
    return mod.susie_multi_ss(
        XTX_list=[np.eye(p) * 100 for _ in range(2)],
        XTY_list=[np.array([5.0, 1.0, 0.0, 0.0]) for _ in range(2)],
        YTY_list=[100.0, 100.0],
        population_sizes=[100, 100],
        rho=np.array([[1.0, 0.75], [0.75, 1.0]]),
        L=1,
        max_iter=len(trace),
        estimate_prior_variance=False,
        estimate_residual_variance=False,
        iter_before_zeroing_effects=-1,
        tol=0.001,
        float_type=np.float64,
        R_list=[np.eye(p), np.eye(p)],
    )


def test_multisusie_wrapper_uses_float64(locus_set_two_pop, monkeypatch):
    """Guard the actual wrapper boundary, not just LD storage dtype."""
    from .test_multisusie import _make_mock_multisusie_rss
    from credtools.constants import ColName
    from credtools.wrappers.multisusie import run_multisusie

    snps = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
    # Use the existing mock result builder, with an explicit boundary assertion.
    mock = _make_mock_multisusie_rss(len(snps))

    def fit(**kwargs):
        assert kwargs["float_type"] is np.float64
        return mock(**kwargs)

    monkeypatch.setattr("credtools.wrappers.multisusie.multisusie_rss", fit)
    result = run_multisusie(locus_set_two_pop)
    assert result.parameters["inference_dtype"] == "float64"
