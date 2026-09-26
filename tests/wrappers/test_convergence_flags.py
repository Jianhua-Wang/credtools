"""Convergence metadata must reflect native evidence without changing outputs."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.wrappers.susiex import run_susiex

from .test_susiex import _make_mock_susiex_run_tool


@pytest.mark.parametrize("marker, expected", [("FAIL", False), ("NULL", True)])
def test_susiex_native_markers(locus_set_two_pop, monkeypatch, marker, expected):
    """Both native no-CS states keep empty CS/PIPs but differ in convergence."""
    snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
    original = _make_mock_susiex_run_tool(snpids, no_cs=True)

    def run(*args, **kwargs):
        original(*args, **kwargs)
        Path(args[3][1]).write_text(marker + "\n")

    monkeypatch.setattr("credtools.wrappers.susiex.tool_manager.run_tool", run)
    result = run_susiex(locus_set_two_pop)
    assert result.converged is expected
    assert result.n_cs == 0
    assert result.snps == []
    assert result.pips.index.tolist() == snpids
    assert result.pips.isna().all()  # Preserve the existing no-CS PIP convention.


def test_susiex_normal_cs_preserves_outputs(locus_set_two_pop, monkeypatch):
    """Adding the flag must not change normal native CS/PIP/purity parsing."""
    snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
    original = _make_mock_susiex_run_tool(snpids, n_cs=2)
    native = {}

    def run(*args, **kwargs):
        original(*args, **kwargs)
        frame = pd.read_csv(args[3][0], sep="\t").set_index("SNP")
        native["pip"] = frame.max(axis=1)

    monkeypatch.setattr("credtools.wrappers.susiex.tool_manager.run_tool", run)
    result = run_susiex(locus_set_two_pop)
    assert result.converged is True
    assert result.snps == [snpids[:3], snpids[3:6]]
    np.testing.assert_array_equal(result.pips.values, native["pip"].values)
    np.testing.assert_allclose(result.purity, [0.81, 0.82])


@pytest.mark.parametrize(
    "content",
    [
        None,
        "",
        "UNKNOWN\n",
        "FAIL\nextra\n",
        "NULL\nextra\n",
        "CS_ID\tSNP\n",
        "CS_ID\tBAD\n1\trs1\n",
        "CS_ID\tSNP\n1\n",
        "CS_ID\tSNP\n1\t\n",
        "CS_ID\tSNP\nnot-a-number\trs1\n",
        "CS_ID\tSNP\n0\trs1\n",
        "CS_ID\tSNP\n1.5\trs1\n",
        "CS_ID\tSNP\n1\trs1\textra\n",
    ],
)
def test_susiex_bad_cs_is_error(locus_set_two_pop, monkeypatch, content):
    """Do not infer convergence from a missing/malformed CS file or SNP shape."""
    snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
    original = _make_mock_susiex_run_tool(snpids, no_cs=True)

    def run(*args, **kwargs):
        original(*args, **kwargs)
        path = Path(args[3][1])
        if content is None:
            path.unlink()
        else:
            path.write_text(content)

    monkeypatch.setattr("credtools.wrappers.susiex.tool_manager.run_tool", run)
    with pytest.raises((FileNotFoundError, ValueError)):
        run_susiex(locus_set_two_pop)


def test_mesusie_original_elbo_rule(tmp_path):
    """Exercise real R, including unused padding, exhaustion and negative deltas."""
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is required for the R adapter unit test")
    wrapper = Path(__file__).parents[2] / "credtools/wrappers/mesusie_wrapper.R"
    check = tmp_path / "check.R"
    check.write_text(
        r"""
args <- commandArgs(trailingOnly = TRUE)
env <- new.env(parent = baseenv())
for (expr in parse(args[1])) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("mesusie_converged_from_elbo"))) {
    eval(expr, env)
  }
}
stopifnot(exists("mesusie_converged_from_elbo", envir = env, inherits = FALSE))
cases <- list(
  list(c(-Inf, -10, -9.9995, NA, NA), TRUE),
  list(c(-Inf, 0, 0.001), FALSE),
  list(c(-Inf, 0, 0.002), FALSE),
  list(c(-Inf, -10, -20, NA), TRUE),
  list(c(-Inf, 0, 0), TRUE),
  list(c(-Inf, -10, -9, -8), FALSE),
  list(c(-Inf, -10), FALSE),
  list(c(-Inf, NA, NA), FALSE),
  list(numeric(0), FALSE),
  list(NULL, FALSE),
  list(c(-Inf, -10, -9, Inf), FALSE),
  list(c(-Inf, -10, -9, -Inf), FALSE),
  list(c(-Inf, -10, -9, NaN, NA), FALSE),
  list(c(-Inf, -10, NA, -9.9995), FALSE),
  list(c(-Inf, 1e308, -1e308), FALSE)
)
for (case in cases) {
  actual <- env$mesusie_converged_from_elbo(case[[1]])
  stopifnot(identical(actual, case[[2]]))
}
"""
    )
    result = subprocess.run(
        [rscript, "--vanilla", str(check), str(wrapper)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
