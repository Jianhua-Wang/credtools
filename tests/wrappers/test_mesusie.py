"""Tests for the MESuSiE multi-ancestry fine-mapping wrapper."""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.mesusie import run_mesusie

from .conftest import _make_locus_set


def _make_mock_mesusie_subprocess(
    snpids,
    n_cs=1,
    has_purity=True,
    no_cs=False,
):
    """Create a mock subprocess.run that writes MESuSiE output files.

    Parameters
    ----------
    snpids : list of str
        SNP IDs for output.
    n_cs : int
        Number of credible sets.
    has_purity : bool
        Whether to write purity values.
    no_cs : bool
        If True, writes output with no credible sets.
    """

    def mock_subprocess_run(cmd, capture_output=True, text=True, check=False, **kw):
        # Find temp_dir from command args
        temp_dir = None
        for i, arg in enumerate(cmd):
            if arg == "--temp_dir":
                temp_dir = cmd[i + 1]
                break

        if temp_dir is None:
            raise ValueError("Could not find --temp_dir in command args")

        os.makedirs(temp_dir, exist_ok=True)

        if no_cs:
            # Write PIPs with low values, no credible sets
            pip_df = pd.DataFrame(
                {
                    "SNP": snpids,
                    "PIP": np.random.uniform(0, 0.1, len(snpids)),
                }
            )
            pip_df.to_csv(f"{temp_dir}/mesusie_pips.csv", index=False)
            # Empty CS file
            pd.DataFrame({"CS_ID": [], "SNP": []}).to_csv(
                f"{temp_dir}/mesusie_cs.csv", index=False
            )
            # Empty purity file
            pd.DataFrame({"CS_ID": [], "PURITY": [], "CS_TYPE": []}).to_csv(
                f"{temp_dir}/mesusie_purity.csv", index=False
            )
        else:
            # Write PIPs
            pip_vals = np.random.uniform(0, 0.5, len(snpids))
            pip_vals[0] = 0.85
            pip_df = pd.DataFrame({"SNP": snpids, "PIP": pip_vals})
            pip_df.to_csv(f"{temp_dir}/mesusie_pips.csv", index=False)

            # Write CS file
            cs_rows = []
            for cs_i in range(1, n_cs + 1):
                start_idx = (cs_i - 1) * 3
                cs_snp_ids = snpids[start_idx : start_idx + 3]
                for snp in cs_snp_ids:
                    cs_rows.append({"CS_ID": cs_i, "SNP": snp})
            cs_df = pd.DataFrame(cs_rows)
            cs_df.to_csv(f"{temp_dir}/mesusie_cs.csv", index=False)

            # Write purity file
            if has_purity:
                purity_rows = []
                cs_types = ["shared", "EUR_specific", "AFR_specific"]
                for cs_i in range(1, n_cs + 1):
                    purity_rows.append(
                        {
                            "CS_ID": cs_i,
                            "PURITY": 0.8 + cs_i * 0.01,
                            "CS_TYPE": cs_types[(cs_i - 1) % len(cs_types)],
                        }
                    )
                purity_df = pd.DataFrame(purity_rows)
                purity_df.to_csv(f"{temp_dir}/mesusie_purity.csv", index=False)
            else:
                pd.DataFrame({"CS_ID": [], "PURITY": [], "CS_TYPE": []}).to_csv(
                    f"{temp_dir}/mesusie_purity.csv", index=False
                )

        # Write converged file
        with open(f"{temp_dir}/mesusie_converged.txt", "w") as f:
            f.write("TRUE")

        return subprocess.CompletedProcess(
            args=cmd, returncode=0, stdout="", stderr=""
        )

    return mock_subprocess_run


class TestRunMesusieBasic:
    """Basic MESuSiE functionality tests."""

    def test_basic_call(self, locus_set_two_pop, monkeypatch):
        """Verify MESuSiE returns a CredibleSet."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.MESUSIE

    def test_output_structure(self, locus_set_two_pop, monkeypatch):
        """Verify MESuSiE result has correct structure."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1
        assert isinstance(result.pips, pd.Series)

    def test_parameters_stored(self, locus_set_two_pop, monkeypatch):
        """Parameters should be stored in the result."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop, max_causal=3, coverage=0.99)
        assert result.parameters["max_causal"] == 3
        assert result.parameters["coverage"] == 0.99

    def test_three_pop(self, locus_set_three_pop, monkeypatch):
        """Verify MESuSiE works with 3 populations."""
        snpids = locus_set_three_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_three_pop)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.MESUSIE


class TestRunMesusieNoCS:
    """Tests for when no credible set is found."""

    def test_no_cs(self, locus_set_two_pop, monkeypatch):
        """Should handle no credible set found."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, no_cs=True),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert result.n_cs == 0
        assert result.snps == []


class TestRunMesuSiePurity:
    """Tests for purity extraction."""

    def test_purity_from_output(self, locus_set_two_pop, monkeypatch):
        """Purity values should be extracted from output."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1, has_purity=True),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert result.purity is not None
        assert len(result.purity) == 1

    def test_purity_none_without_values(self, locus_set_two_pop, monkeypatch):
        """Purity should be None when no purity data exists."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1, has_purity=False),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        if result.n_cs > 0:
            for p in result.purity:
                assert p is None


class TestRunMesusieMultiCS:
    """Tests for multiple credible sets."""

    def test_multi_cs(self, locus_set_two_pop, monkeypatch):
        """Should handle multiple credible sets correctly."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=2),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert result.n_cs == 2
        assert len(result.snps) == 2
        assert len(result.lead_snps) == 2

    def test_multi_cs_purity(self, locus_set_two_pop, monkeypatch):
        """Multiple CS should have purity values."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=2, has_purity=True),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert result.purity is not None
        assert len(result.purity) == 2

    def test_cs_types_stored(self, locus_set_two_pop, monkeypatch):
        """CS types (shared/specific) should be stored in parameters."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=2, has_purity=True),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        assert "cs_types" in result.parameters
        assert len(result.parameters["cs_types"]) == 2


class TestRunMesusieLeadSNP:
    """Tests for lead SNP extraction."""

    def test_lead_snp_is_max_pip(self, locus_set_two_pop, monkeypatch):
        """Lead SNP should be the one with highest PIP in CS."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            _make_mock_mesusie_subprocess(snpids, n_cs=1),
        )
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )
        result = run_mesusie(locus_set_two_pop)
        if result.n_cs > 0:
            cs_snps = result.snps[0]
            lead = result.lead_snps[0]
            assert lead == result.pips[result.pips.index.isin(cs_snps)].idxmax()


class TestRunMesusieErrorHandling:
    """Tests for error handling."""

    def test_r_not_installed(self, locus_set_two_pop, monkeypatch):
        """Should raise FileNotFoundError when R is not installed."""
        monkeypatch.setattr("shutil.which", lambda x: None)
        with pytest.raises(FileNotFoundError, match="Rscript"):
            run_mesusie(locus_set_two_pop)

    def test_mesusie_not_installed(self, locus_set_two_pop, monkeypatch):
        """Should raise FileNotFoundError when MESuSiE R package is not installed."""
        monkeypatch.setattr("shutil.which", lambda x: "/usr/bin/Rscript")

        def mock_subprocess_run(cmd, **kw):
            return subprocess.CompletedProcess(
                args=cmd, returncode=1, stdout="", stderr="there is no package"
            )

        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            mock_subprocess_run,
        )
        with pytest.raises(FileNotFoundError, match="MESuSiE"):
            run_mesusie(locus_set_two_pop)

    def test_subprocess_error(self, locus_set_two_pop, monkeypatch):
        """Should raise RuntimeError when R script fails."""
        monkeypatch.setattr(
            "credtools.wrappers.mesusie._check_r_and_mesusie", lambda: None
        )

        def mock_subprocess_run(cmd, **kw):
            return subprocess.CompletedProcess(
                args=cmd, returncode=1, stdout="", stderr="Error in meSuSie_core"
            )

        monkeypatch.setattr(
            "credtools.wrappers.mesusie.subprocess.run",
            mock_subprocess_run,
        )
        with pytest.raises(RuntimeError, match="MESuSiE R script failed"):
            run_mesusie(locus_set_two_pop)
