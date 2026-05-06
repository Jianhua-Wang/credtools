"""Tests for the SuSiE-inf fine-mapping wrapper."""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.susie_inf import run_susie_inf

from .conftest import _make_locus


def _make_mock_susie_inf_subprocess(
    snpids,
    n_cs=1,
    no_cs=False,
    converged=True,
    n_iter=42,
):
    """Create a mock subprocess.run that writes SuSiE-inf output files."""

    def mock_subprocess_run(cmd, capture_output=True, text=True, check=False, **kw):
        temp_dir = None
        for i, arg in enumerate(cmd):
            if arg == "--temp_dir":
                temp_dir = cmd[i + 1]
                break
        if temp_dir is None:
            raise ValueError("Could not find --temp_dir in command args")

        os.makedirs(temp_dir, exist_ok=True)

        rng = np.random.default_rng(0)
        if no_cs:
            pip_df = pd.DataFrame(
                {
                    "SNP": snpids,
                    "PIP": rng.uniform(0, 0.1, len(snpids)),
                }
            )
            pip_df.to_csv(f"{temp_dir}/susie_inf_pips.csv", index=False)
            pd.DataFrame({"CS_ID": [], "SNP": []}).to_csv(
                f"{temp_dir}/susie_inf_cs.csv", index=False
            )
        else:
            pip_vals = rng.uniform(0, 0.5, len(snpids))
            pip_vals[0] = 0.91
            pip_df = pd.DataFrame({"SNP": snpids, "PIP": pip_vals})
            pip_df.to_csv(f"{temp_dir}/susie_inf_pips.csv", index=False)

            cs_rows = []
            for cs_i in range(1, n_cs + 1):
                start_idx = (cs_i - 1) * 3
                cs_snp_ids = snpids[start_idx : start_idx + 3]
                for snp in cs_snp_ids:
                    cs_rows.append({"CS_ID": cs_i, "SNP": snp})
            pd.DataFrame(cs_rows).to_csv(f"{temp_dir}/susie_inf_cs.csv", index=False)

        with open(f"{temp_dir}/susie_inf_status.txt", "w") as f:
            f.write(f"converged\t{'TRUE' if converged else 'FALSE'}\n")
            f.write(f"niter\t{n_iter}\n")

        return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="", stderr="")

    return mock_subprocess_run


@pytest.fixture
def locus_significant_for_susie_inf():
    return _make_locus(p_range=(1e-12, 1e-9), n_snps=20)


@pytest.fixture
def locus_no_significant_for_susie_inf():
    return _make_locus(p_range=(0.1, 0.9), n_snps=20)


@pytest.fixture
def locus_unmatched_for_susie_inf():
    return _make_locus(matched=False, p_range=(1e-12, 1e-9), n_snps=20)


class TestRunSusieInfBasic:
    """Basic SuSiE-inf functionality tests."""

    def test_basic_call(self, locus_significant_for_susie_inf, monkeypatch):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_significant_for_susie_inf)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.SUSIE_INF

    def test_output_structure(self, locus_significant_for_susie_inf, monkeypatch):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_significant_for_susie_inf)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1
        assert isinstance(result.pips, pd.Series)
        assert len(result.pips) == len(snpids)

    def test_parameters_stored(self, locus_significant_for_susie_inf, monkeypatch):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(
            locus_significant_for_susie_inf,
            max_causal=4,
            coverage=0.99,
            max_iter=200,
            estimate_residual_variance=True,
        )
        assert result.parameters["max_causal"] == 4
        assert result.parameters["coverage"] == 0.99
        assert result.parameters["max_iter"] == 200
        assert result.parameters["estimate_residual_variance"] is True
        assert result.parameters["unmappable_effects"] == "inf"


class TestRunSusieInfSignificanceFilter:
    """Tests for significance threshold short-circuit."""

    def test_no_significant_returns_empty(
        self, locus_no_significant_for_susie_inf, monkeypatch
    ):
        called = {"count": 0}

        def mock_run(cmd, **kw):
            called["count"] += 1
            return subprocess.CompletedProcess(args=cmd, returncode=0)

        monkeypatch.setattr("credtools.wrappers.susie_inf.subprocess.run", mock_run)
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_no_significant_for_susie_inf)
        assert result.n_cs == 0
        assert result.snps == []
        assert (result.pips == 0).all()
        assert called["count"] == 0


class TestRunSusieInfNoCS:
    """Tests for when SuSiE-inf finds no credible sets."""

    def test_no_cs(self, locus_significant_for_susie_inf, monkeypatch):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, no_cs=True),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_significant_for_susie_inf)
        assert result.n_cs == 0
        assert result.snps == []


class TestRunSusieInfMultiCS:
    """Tests for multiple credible sets."""

    def test_multi_cs(self, locus_significant_for_susie_inf, monkeypatch):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=2),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_significant_for_susie_inf, max_causal=5)
        assert result.n_cs == 2
        assert len(result.snps) == 2
        assert len(result.lead_snps) == 2


class TestRunSusieInfLeadSNP:
    """Tests for lead SNP extraction."""

    def test_lead_snp_is_max_pip(self, locus_significant_for_susie_inf, monkeypatch):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_significant_for_susie_inf)
        cs_snps = result.snps[0]
        assert (
            result.lead_snps[0] == result.pips[result.pips.index.isin(cs_snps)].idxmax()
        )


class TestRunSusieInfConvergence:
    """Tests for convergence metadata propagation."""

    def test_converged_true_recorded(
        self, locus_significant_for_susie_inf, monkeypatch
    ):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1, converged=True, n_iter=12),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_significant_for_susie_inf)
        assert result.converged is True
        assert result.n_iter == 12

    def test_nonconverged_returns_empty_when_requested(
        self, locus_significant_for_susie_inf, monkeypatch
    ):
        snpids = locus_significant_for_susie_inf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1, converged=False, n_iter=99),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(
            locus_significant_for_susie_inf, empty_on_nonconvergence=True
        )
        assert result.n_cs == 0
        assert result.converged is False
        assert (result.pips == 0).all()


class TestRunSusieInfUnmatched:
    """Tests for the unmatched LD/sumstat path."""

    def test_unmatched_runs(self, locus_unmatched_for_susie_inf, monkeypatch):
        from credtools.locus import intersect_sumstat_ld

        matched_locus = intersect_sumstat_ld(locus_unmatched_for_susie_inf)
        snpids = matched_locus.sumstats[ColName.SNPID].tolist()

        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run",
            _make_mock_susie_inf_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )
        result = run_susie_inf(locus_unmatched_for_susie_inf)
        assert isinstance(result, CredibleSet)
        assert result.n_cs == 1


class TestRunSusieInfErrorHandling:
    """Tests for error handling."""

    def test_r_not_installed(self, locus_significant_for_susie_inf, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda x: None)
        with pytest.raises(FileNotFoundError, match="Rscript"):
            run_susie_inf(locus_significant_for_susie_inf)

    def test_susier_not_installed(self, locus_significant_for_susie_inf, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda x: "/usr/bin/Rscript")

        def mock_subprocess_run(cmd, **kw):
            return subprocess.CompletedProcess(
                args=cmd,
                returncode=1,
                stdout="",
                stderr="there is no package called 'susieR'",
            )

        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run", mock_subprocess_run
        )
        with pytest.raises(FileNotFoundError, match="susieR"):
            run_susie_inf(locus_significant_for_susie_inf)

    def test_subprocess_error(self, locus_significant_for_susie_inf, monkeypatch):
        monkeypatch.setattr(
            "credtools.wrappers.susie_inf._check_r_and_susie_inf", lambda: None
        )

        def mock_subprocess_run(cmd, **kw):
            return subprocess.CompletedProcess(
                args=cmd, returncode=1, stdout="", stderr="Error in susie_rss"
            )

        monkeypatch.setattr(
            "credtools.wrappers.susie_inf.subprocess.run", mock_subprocess_run
        )
        with pytest.raises(RuntimeError, match="SuSiE-inf R script failed"):
            run_susie_inf(locus_significant_for_susie_inf)
