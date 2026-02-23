"""Tests for the SuSiEx multi-ancestry fine-mapping wrapper."""

import os

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.susiex import run_susiex

from .conftest import _make_locus_set


def _make_mock_susiex_run_tool(
    snpids,
    n_cs=1,
    has_purity_summary=True,
    no_cs=False,
):
    """Create a mock run_tool that writes SuSiEx output files.

    Parameters
    ----------
    snpids : list of str
        SNP IDs for output.
    n_cs : int
        Number of credible sets.
    has_purity_summary : bool
        Whether to write a summary file with purity.
    no_cs : bool
        If True, writes a .snp file with only 2 columns (no CS found).
    """

    def mock_run_tool(tool_name, cmd, log_file, required_output_files, timeout=None):
        td = os.path.dirname(log_file)

        # Parse output prefix from required_output_files
        snp_file = required_output_files[0]
        cs_file = required_output_files[1]
        prefix = snp_file.replace(".snp", "")

        if no_cs:
            # Only SNP and PIP columns (no credible set found)
            pip_df = pd.DataFrame(
                {
                    "SNP": snpids,
                    "PIP": np.random.uniform(0, 0.1, len(snpids)),
                }
            )
            pip_df.to_csv(snp_file, sep="\t", index=False)
            # Still need .cs file to exist
            pd.DataFrame({"CS_ID": [], "SNP": []}).to_csv(
                cs_file, sep="\t", index=False
            )
        else:
            # Write .snp file with PIP columns
            # SuSiEx outputs SNP + PIP_all + PIP1 [+ PIP2 ...] when CS found
            # Must have > 2 columns to not be treated as "no CS"
            pip_data = {"SNP": snpids}
            pip_vals_all = np.random.uniform(0, 0.5, len(snpids))
            pip_vals_all[0] = 0.85
            pip_data["PIP_all"] = pip_vals_all
            for cs_i in range(1, n_cs + 1):
                pip_data[f"PIP{cs_i}"] = np.random.uniform(0, 0.5, len(snpids))
            pip_data["PIP1"][0] = 0.85
            pip_df = pd.DataFrame(pip_data)
            pip_df.to_csv(snp_file, sep="\t", index=False)

            # Write .cs file
            cs_rows = []
            for cs_i in range(1, n_cs + 1):
                start_idx = (cs_i - 1) * 3
                cs_snp_ids = snpids[start_idx : start_idx + 3]
                for snp in cs_snp_ids:
                    cs_rows.append({"CS_ID": cs_i, "SNP": snp})
            cs_df = pd.DataFrame(cs_rows)
            cs_df.to_csv(cs_file, sep="\t", index=False)

            # Write .summary file for purity
            if has_purity_summary:
                summary_rows = []
                for cs_i in range(1, n_cs + 1):
                    summary_rows.append({"CS_ID": cs_i, "CS_PURITY": 0.8 + cs_i * 0.01})
                summary_df = pd.DataFrame(summary_rows)
                summary_df.to_csv(f"{prefix}.summary", sep="\t", index=False)

    return mock_run_tool


class TestRunSusiexBasic:
    """Basic SuSiEx functionality tests."""

    def test_basic_call(self, locus_set_two_pop, monkeypatch):
        """Verify SuSiEx returns a CredibleSet."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=1),
        )
        result = run_susiex(locus_set_two_pop)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.SUSIEX

    def test_output_structure(self, locus_set_two_pop, monkeypatch):
        """Verify SuSiEx result has correct structure."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=1),
        )
        result = run_susiex(locus_set_two_pop)
        assert result.n_cs == 1
        assert len(result.snps) == 1
        assert len(result.lead_snps) == 1
        assert isinstance(result.pips, pd.Series)

    def test_parameters_stored(self, locus_set_two_pop, monkeypatch):
        """Parameters should be stored in the result."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=1),
        )
        result = run_susiex(locus_set_two_pop, max_causal=3, coverage=0.99)
        assert result.parameters["max_causal"] == 3
        assert result.parameters["coverage"] == 0.99


class TestRunSusiexNoCS:
    """Tests for when no credible set is found."""

    def test_no_cs_two_columns(self, locus_set_two_pop, monkeypatch):
        """Should handle 2-column .snp file (no CS found)."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, no_cs=True),
        )
        result = run_susiex(locus_set_two_pop)
        assert result.n_cs == 0
        assert result.snps == []


class TestRunSusiexPurity:
    """Tests for purity extraction."""

    def test_purity_from_summary(self, locus_set_two_pop, monkeypatch):
        """Purity values should be extracted from .summary file."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=1, has_purity_summary=True),
        )
        result = run_susiex(locus_set_two_pop)
        assert result.purity is not None
        assert len(result.purity) == 1

    def test_purity_none_without_summary(self, locus_set_two_pop, monkeypatch):
        """Purity should be None when no summary file exists."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=1, has_purity_summary=False),
        )
        result = run_susiex(locus_set_two_pop)
        if result.n_cs > 0:
            # Purity values should be None for each CS
            for p in result.purity:
                assert p is None


class TestRunSusiexMultiCS:
    """Tests for multiple credible sets."""

    def test_multi_cs(self, locus_set_two_pop, monkeypatch):
        """Should handle multiple credible sets correctly."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=2),
        )
        result = run_susiex(locus_set_two_pop)
        assert result.n_cs == 2
        assert len(result.snps) == 2
        assert len(result.lead_snps) == 2

    def test_multi_cs_purity(self, locus_set_two_pop, monkeypatch):
        """Multiple CS should have purity values."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=2, has_purity_summary=True),
        )
        result = run_susiex(locus_set_two_pop)
        assert result.purity is not None
        assert len(result.purity) == 2


class TestRunSusiexLeadSNP:
    """Tests for lead SNP extraction."""

    def test_lead_snp_is_max_pip(self, locus_set_two_pop, monkeypatch):
        """Lead SNP should be the one with highest PIP in CS."""
        snpids = locus_set_two_pop.loci[0].sumstats[ColName.SNPID].tolist()
        monkeypatch.setattr(
            "credtools.wrappers.susiex.tool_manager.run_tool",
            _make_mock_susiex_run_tool(snpids, n_cs=1),
        )
        result = run_susiex(locus_set_two_pop)
        if result.n_cs > 0:
            cs_snps = result.snps[0]
            lead = result.lead_snps[0]
            assert lead == result.pips[result.pips.index.isin(cs_snps)].idxmax()
