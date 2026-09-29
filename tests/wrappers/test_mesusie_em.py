"""Exercise the bundled adapter against an installed MESuSiE, when available."""

import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def rscript_with_mesusie():
    """Skip only when the optional R backend is genuinely unavailable."""
    r = shutil.which("Rscript")
    if not r:
        pytest.skip("Rscript is unavailable")
    probe = subprocess.run(
        [
            r,
            "--vanilla",
            "-e",
            'quit(status=if(all(vapply(c("MESuSiE","Rcpp","RcppArmadillo"),requireNamespace,logical(1),quietly=TRUE)))0 else 1)',
        ],
        capture_output=True,
        text=True,
    )
    if probe.returncode:
        pytest.skip("MESuSiE/Rcpp/RcppArmadillo are unavailable")
    return r


def test_real_mesusie_em_permutations(tmp_path, rscript_with_mesusie):
    r = rscript_with_mesusie
    root = Path(__file__).parents[2]
    result = subprocess.run(
        [
            r,
            "--vanilla",
            str(Path(__file__).with_suffix(".R")),
            str(root / "credtools/wrappers"),
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "MESUSIE_EM_INTEGRATION_PASS" in result.stdout


def test_python_to_r_em_adapter(tmp_path, monkeypatch, rscript_with_mesusie):
    """Exercise actual subprocess parameters, CSs and diagnostics end-to-end."""
    import numpy as np
    from scipy.stats import norm
    from credtools.constants import ColName
    from credtools.wrappers.mesusie import run_mesusie
    from .conftest import _make_locus_set

    monkeypatch.chdir(tmp_path)
    loci = _make_locus_set(n_pop=5, n_snps=12)
    for i, locus in enumerate(loci.loci):
        z = np.array([8.0 + i / 10] + [0.1] * 11)
        locus.sumstats[ColName.BETA] = z / 100
        locus.sumstats[ColName.SE] = 0.01
        locus.sumstats[ColName.P] = 2 * norm.sf(abs(z))
        locus.ld.r[:] = np.eye(12)
    result = run_mesusie(loci, max_causal=2, coverage=0.9, max_iter=50, tol=0.002)
    assert result.converged is True
    assert result.n_iter is not None and 1 < result.n_iter <= 50
    assert result.n_cs == 1
    assert result.parameters["mesusie_optimizer"] == "em"
    status = result.parameters["mesusie_runtime"]
    assert float(status["coverage"]) == 0.9
    assert float(status["outer_tol"]) == 0.002
    assert int(status["inner_calls"]) > 0
    assert np.isfinite(result.pips).all()
