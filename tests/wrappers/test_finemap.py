"""Tests for the FINEMAP fine-mapping wrapper."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.finemap import run_finemap

from .conftest import _make_locus


def _make_mock_finemap_run_tool(temp_dir_ref, snpids, n_causal=1):
    """Create a mock run_tool that writes FINEMAP output files.

    Parameters
    ----------
    temp_dir_ref : list
        Single-element list to capture temp_dir from the call.
    snpids : list of str
        SNP IDs to include in outputs.
    n_causal : int
        Number of causal SNPs to simulate.
    """

    def mock_run_tool(tool_name, cmd, log_file, required_output_files, timeout=None):
        # Determine temp_dir from the log_file path
        td = os.path.dirname(log_file)
        temp_dir_ref.append(td)

        # Write .snp file
        snp_data = pd.DataFrame(
            {
                "index": range(len(snpids)),
                "rsid": snpids,
                "chromosome": [1] * len(snpids),
                "position": [1000 + i * 100 for i in range(len(snpids))],
                "allele1": ["A"] * len(snpids),
                "allele2": ["G"] * len(snpids),
                "maf": [0.3] * len(snpids),
                "beta": [0.1] * len(snpids),
                "se": [0.02] * len(snpids),
                "prob": np.zeros(len(snpids)),
                "log10bf": np.zeros(len(snpids)),
                "mean": np.zeros(len(snpids)),
                "sd": np.zeros(len(snpids)),
                "mean_incl": np.zeros(len(snpids)),
                "sd_incl": np.zeros(len(snpids)),
            }
        )
        snp_data.to_csv(f"{td}/finemap.snp", sep=" ", index=False)

        # Write .config file
        with open(f"{td}/finemap.config", "w") as f:
            f.write("rank config log10bf prob\n")
            f.write("1 1 0.5 0.8\n")

        # Write .cred files for each causal
        for nc in range(1, n_causal + 1):
            n_for_cred = min(nc, len(snpids))
            cred_data = {"index": range(max(3, n_for_cred))}
            for ci in range(1, nc + 1):
                cred_col = [None] * max(3, n_for_cred)
                prob_col = [None] * max(3, n_for_cred)
                # Pick a SNP for each cred set
                snp_idx = min(ci - 1, len(snpids) - 1)
                cred_col[0] = snpids[snp_idx]
                prob_col[0] = 0.9
                if n_for_cred > 1 and ci == 1:
                    next_idx = min(snp_idx + 1, len(snpids) - 1)
                    cred_col[1] = snpids[next_idx]
                    prob_col[1] = 0.05
                cred_data[f"cred{ci}"] = cred_col
                cred_data[f"prob{ci}"] = prob_col

            cred_df = pd.DataFrame(cred_data)
            with open(f"{td}/finemap.cred{nc}", "w") as f:
                f.write(f"# posterior probability of {nc} causal SNP(s) = 0.80\n")
            cred_df.to_csv(f"{td}/finemap.cred{nc}", sep=" ", index=False, mode="a")

    return mock_run_tool


def _make_mock_finemap_empty_cred(temp_dir_ref, snpids):
    """Mock that produces no .cred files."""

    def mock_run_tool(tool_name, cmd, log_file, required_output_files, timeout=None):
        td = os.path.dirname(log_file)
        temp_dir_ref.append(td)

        snp_data = pd.DataFrame(
            {
                "rsid": snpids,
                "chromosome": [1] * len(snpids),
                "position": [1000 + i * 100 for i in range(len(snpids))],
            }
        )
        snp_data.to_csv(f"{td}/finemap.snp", sep=" ", index=False)

        with open(f"{td}/finemap.config", "w") as f:
            f.write("rank config log10bf prob\n")

    return mock_run_tool


class TestRunFinemapBasic:
    """Basic FINEMAP functionality tests."""

    def test_basic_call(self, locus_with_maf, monkeypatch):
        """FINEMAP returns a CredibleSet."""
        temp_dir_ref = []
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids, n_causal=1),
        )
        result = run_finemap(locus_with_maf)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.FINEMAP

    def test_output_structure(self, locus_with_maf, monkeypatch):
        """FINEMAP result has correct structure."""
        temp_dir_ref = []
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids, n_causal=1),
        )
        result = run_finemap(locus_with_maf)
        assert result.n_cs >= 0
        assert isinstance(result.pips, pd.Series)
        assert result.coverage == 0.95

    def test_parameters_stored(self, locus_with_maf, monkeypatch):
        """Parameters should be stored in the result."""
        temp_dir_ref = []
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids),
        )
        result = run_finemap(locus_with_maf, max_causal=3, n_iter=50000)
        assert result.parameters["max_causal"] == 3
        assert result.parameters["n_iter"] == 50000


class TestRunFinemapEdgeCases:
    """Edge case tests for FINEMAP."""

    def test_no_significant_snps(self, locus_no_significant):
        """Should return empty result (early return before tool_manager)."""
        # Need MAF column for this test
        locus_no_significant.sumstats[ColName.MAF] = 0.3
        result = run_finemap(locus_no_significant)
        assert result.n_cs == 0
        assert (result.pips == 0).all()

    def test_missing_maf_raises(self, locus_without_maf, monkeypatch):
        """Should raise ValueError when MAF column is missing."""
        temp_dir_ref = []
        snpids = locus_without_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids),
        )
        with pytest.raises(ValueError, match="MAF"):
            run_finemap(locus_without_maf)

    def test_timeout_zero_raises(self, locus_with_maf):
        """Should raise ValueError for non-positive timeout."""
        with pytest.raises(ValueError, match="timeout_minutes"):
            run_finemap(locus_with_maf, timeout_minutes=0)

    def test_timeout_negative_raises(self, locus_with_maf):
        """Should raise ValueError for negative timeout."""
        with pytest.raises(ValueError, match="timeout_minutes"):
            run_finemap(locus_with_maf, timeout_minutes=-5)

    def test_unmatched_ld(self, monkeypatch):
        """Unmatched locus should be intersected before running."""
        locus = _make_locus(matched=False, add_maf=True, p_range=(1e-12, 1e-9))
        temp_dir_ref = []
        snpids = locus.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids),
        )
        result = run_finemap(locus)
        assert isinstance(result, CredibleSet)

    def test_empty_cred_output(self, locus_with_maf, monkeypatch):
        """Should handle empty FINEMAP output (no cred files)."""
        temp_dir_ref = []
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_empty_cred(temp_dir_ref, snpids),
        )
        result = run_finemap(locus_with_maf)
        assert result.n_cs == 0

    def test_multi_causal(self, locus_with_maf, monkeypatch):
        """Should handle multiple causal SNPs."""
        temp_dir_ref = []
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids, n_causal=2),
        )
        result = run_finemap(locus_with_maf, max_causal=2)
        assert isinstance(result, CredibleSet)

    def test_purity_calculated(self, locus_with_maf, monkeypatch):
        """Purity should be calculated when CS found and LD available."""
        temp_dir_ref = []
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool",
            _make_mock_finemap_run_tool(temp_dir_ref, snpids, n_causal=1),
        )
        result = run_finemap(locus_with_maf)
        if result.n_cs > 0:
            assert result.purity is not None

    def test_maf_zero_replacement(self, locus_with_maf, monkeypatch):
        """MAF=0 should be replaced with 0.00001."""
        temp_dir_ref = []
        locus_with_maf.sumstats.loc[0, ColName.MAF] = 0.0
        snpids = locus_with_maf.sumstats[ColName.SNPID].tolist()

        def mock_and_check(tool_name, cmd, log_file, required_files, timeout=None):
            td = os.path.dirname(log_file)
            # Read the z file to verify MAF replacement
            z_file = f"{td}/finemap.z"
            if os.path.exists(z_file):
                z_data = pd.read_csv(z_file, sep=" ")
                assert (z_data["maf"] > 0).all()
            _make_mock_finemap_run_tool(temp_dir_ref, snpids)(
                tool_name, cmd, log_file, required_files, timeout
            )

        monkeypatch.setattr(
            "credtools.wrappers.finemap.tool_manager.run_tool", mock_and_check
        )
        result = run_finemap(locus_with_maf)
        assert isinstance(result, CredibleSet)
