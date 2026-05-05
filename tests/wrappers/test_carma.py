"""Tests for the CARMA fine-mapping wrapper."""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.carma import run_carma

from .conftest import _make_locus


def _make_mock_carma_subprocess(
    snpids,
    n_cs=1,
    has_outliers=False,
    no_cs=False,
):
    """Create a mock subprocess.run that writes CARMA output files."""

    def mock_subprocess_run(cmd, capture_output=True, text=True, check=False, **kw):
        temp_dir = None
        for i, arg in enumerate(cmd):
            if arg == "--temp_dir":
                temp_dir = cmd[i + 1]
                break
        if temp_dir is None:
            raise ValueError("Could not find --temp_dir in command args")

        os.makedirs(temp_dir, exist_ok=True)

        if no_cs:
            pip_df = pd.DataFrame(
                {
                    "SNP": snpids,
                    "PIP": np.random.uniform(0, 0.1, len(snpids)),
                }
            )
            pip_df.to_csv(f"{temp_dir}/carma_pips.csv", index=False)
            pd.DataFrame({"CS_ID": [], "SNP": []}).to_csv(
                f"{temp_dir}/carma_cs.csv", index=False
            )
            pd.DataFrame({"SNP": []}).to_csv(
                f"{temp_dir}/carma_outliers.csv", index=False
            )
        else:
            pip_vals = np.random.uniform(0, 0.5, len(snpids))
            pip_vals[0] = 0.92
            pip_df = pd.DataFrame({"SNP": snpids, "PIP": pip_vals})
            pip_df.to_csv(f"{temp_dir}/carma_pips.csv", index=False)

            cs_rows = []
            for cs_i in range(1, n_cs + 1):
                start_idx = (cs_i - 1) * 3
                cs_snp_ids = snpids[start_idx : start_idx + 3]
                for snp in cs_snp_ids:
                    cs_rows.append({"CS_ID": cs_i, "SNP": snp})
            cs_df = pd.DataFrame(cs_rows)
            cs_df.to_csv(f"{temp_dir}/carma_cs.csv", index=False)

            if has_outliers:
                pd.DataFrame({"SNP": [snpids[-1]]}).to_csv(
                    f"{temp_dir}/carma_outliers.csv", index=False
                )
            else:
                pd.DataFrame({"SNP": []}).to_csv(
                    f"{temp_dir}/carma_outliers.csv", index=False
                )

        return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="", stderr="")

    return mock_subprocess_run


@pytest.fixture
def locus_significant_for_carma():
    """Locus with significant SNPs and matched LD/sumstats."""
    return _make_locus(p_range=(1e-12, 1e-9), n_snps=20)


@pytest.fixture
def locus_no_significant_for_carma():
    """Locus where no SNPs reach significance threshold."""
    return _make_locus(p_range=(0.1, 0.9), n_snps=20)


@pytest.fixture
def locus_unmatched_for_carma():
    """Locus where sumstats and LD have different orders."""
    return _make_locus(matched=False, p_range=(1e-12, 1e-9), n_snps=20)


class TestRunCarmaBasic:
    """Basic CARMA functionality tests."""

    def test_basic_call(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.CARMA

    def test_output_structure(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1
        assert isinstance(result.pips, pd.Series)
        assert len(result.pips) == len(snpids)

    def test_parameters_stored(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(
            locus_significant_for_carma,
            max_causal=3,
            coverage=0.99,
            outlier_switch=False,
        )
        assert result.parameters["max_causal"] == 3
        assert result.parameters["coverage"] == 0.99
        assert result.parameters["outlier_switch"] is False


class TestRunCarmaSignificanceFilter:
    """Tests for significance threshold short-circuit."""

    def test_no_significant_returns_empty(
        self, locus_no_significant_for_carma, monkeypatch
    ):
        called = {"count": 0}

        def mock_run(cmd, **kw):
            called["count"] += 1
            return subprocess.CompletedProcess(args=cmd, returncode=0)

        monkeypatch.setattr("credtools.wrappers.carma.subprocess.run", mock_run)
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_no_significant_for_carma)
        assert result.n_cs == 0
        assert result.snps == []
        assert (result.pips == 0).all()
        assert called["count"] == 0


class TestRunCarmaNoCS:
    """Tests for when CARMA finds no credible sets."""

    def test_no_cs(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, no_cs=True),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma)
        assert result.n_cs == 0
        assert result.snps == []


class TestRunCarmaMultiCS:
    """Tests for multiple credible sets."""

    def test_multi_cs(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=2),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma, max_causal=5)
        assert result.n_cs == 2
        assert len(result.snps) == 2
        assert len(result.lead_snps) == 2


class TestRunCarmaLeadSNP:
    """Tests for lead SNP extraction."""

    def test_lead_snp_is_max_pip(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma)
        cs_snps = result.snps[0]
        assert (
            result.lead_snps[0] == result.pips[result.pips.index.isin(cs_snps)].idxmax()
        )


class TestRunCarmaOutliers:
    """Tests for outlier detection output."""

    def test_outliers_stored_in_parameters(
        self, locus_significant_for_carma, monkeypatch
    ):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1, has_outliers=True),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma)
        assert "outliers" in result.parameters
        assert snpids[-1] in result.parameters["outliers"]

    def test_no_outliers_recorded(self, locus_significant_for_carma, monkeypatch):
        snpids = locus_significant_for_carma.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1, has_outliers=False),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_significant_for_carma)
        assert result.parameters.get("outliers", []) == []


class TestRunCarmaUnmatched:
    """Tests for the unmatched LD/sumstat path."""

    def test_unmatched_runs(self, locus_unmatched_for_carma, monkeypatch):
        # After intersect_sumstat_ld, snpids will come from the matched sumstats.
        # Compute them after the alignment to mock the right output.
        from credtools.locus import intersect_sumstat_ld

        matched_locus = intersect_sumstat_ld(locus_unmatched_for_carma)
        snpids = matched_locus.sumstats[ColName.SNPID].tolist()

        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run",
            _make_mock_carma_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)
        result = run_carma(locus_unmatched_for_carma)
        assert isinstance(result, CredibleSet)
        assert result.n_cs == 1


class TestRunCarmaErrorHandling:
    """Tests for error handling."""

    def test_r_not_installed(self, locus_significant_for_carma, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda x: None)
        with pytest.raises(FileNotFoundError, match="Rscript"):
            run_carma(locus_significant_for_carma)

    def test_carma_not_installed(self, locus_significant_for_carma, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda x: "/usr/bin/Rscript")

        def mock_subprocess_run(cmd, **kw):
            return subprocess.CompletedProcess(
                args=cmd, returncode=1, stdout="", stderr="there is no package"
            )

        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run", mock_subprocess_run
        )
        with pytest.raises(FileNotFoundError, match="CARMA"):
            run_carma(locus_significant_for_carma)

    def test_subprocess_error(self, locus_significant_for_carma, monkeypatch):
        monkeypatch.setattr("credtools.wrappers.carma._check_r_and_carma", lambda: None)

        def mock_subprocess_run(cmd, **kw):
            return subprocess.CompletedProcess(
                args=cmd, returncode=1, stdout="", stderr="Error in CARMA"
            )

        monkeypatch.setattr(
            "credtools.wrappers.carma.subprocess.run", mock_subprocess_run
        )
        with pytest.raises(RuntimeError, match="CARMA R script failed"):
            run_carma(locus_significant_for_carma)
