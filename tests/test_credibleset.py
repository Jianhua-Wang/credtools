"""Tests for credtools.credibleset module."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.credibleset import (
    CredibleSet,
    filter_credset_by_purity,
    generate_cs_summary,
)
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_credset(
    n_cs=2,
    tool="susie",
    purity=None,
    per_locus_results=None,
):
    """Create a CredibleSet for testing."""
    snps = [[f"snp_{i}_{j}" for j in range(3)] for i in range(n_cs)]
    lead_snps = [s[0] for s in snps]
    cs_sizes = [len(s) for s in snps]
    all_snpids = [snp for cs in snps for snp in cs]
    pips = pd.Series(
        np.linspace(0.9, 0.1, len(all_snpids)),
        index=all_snpids,
    )
    return CredibleSet(
        tool=tool,
        parameters={"max_causal": 5, "coverage": 0.95},
        coverage=0.95,
        n_cs=n_cs,
        cs_sizes=cs_sizes,
        lead_snps=lead_snps,
        snps=snps,
        pips=pips,
        purity=purity,
        per_locus_results=per_locus_results,
    )


# ---------------------------------------------------------------------------
# TestCredibleSetCopy
# ---------------------------------------------------------------------------
class TestCredibleSetCopy:
    """Tests for CredibleSet.copy method."""

    def test_basic_copy(self):
        cs = _make_credset()
        copied = cs.copy()
        assert copied.tool == cs.tool
        assert copied.n_cs == cs.n_cs
        assert copied.cs_sizes == cs.cs_sizes
        assert copied.lead_snps == cs.lead_snps
        assert copied.snps == cs.snps
        assert copied.pips.equals(cs.pips)
        # Ensure it's a deep copy
        copied._snps[0].append("extra")
        assert "extra" not in cs.snps[0]

    def test_copy_with_per_locus_results(self):
        inner = _make_credset(n_cs=1, tool="abf")
        cs = _make_credset(per_locus_results={"locus_1": inner})
        copied = cs.copy()
        assert "locus_1" in copied.per_locus_results
        assert copied.per_locus_results["locus_1"].tool == "abf"

    def test_self_referencing_copy(self):
        cs = _make_credset()
        cs.set_per_locus_results({"self": cs})
        copied = cs.copy()
        # The self-referencing copy should point to the new copied object
        assert copied.per_locus_results["self"] is copied

    def test_copy_with_purity(self):
        cs = _make_credset(purity=[0.8, 0.6])
        copied = cs.copy()
        assert copied.purity == [0.8, 0.6]
        copied.purity[0] = 0.99
        assert cs.purity[0] == 0.8  # original unchanged


# ---------------------------------------------------------------------------
# TestCredibleSetToDict
# ---------------------------------------------------------------------------
class TestCredibleSetToDict:
    """Tests for CredibleSet.to_dict method."""

    def test_to_dict_keys(self):
        cs = _make_credset(purity=[0.9, 0.7])
        d = cs.to_dict()
        expected_keys = {"tool", "n_cs", "coverage", "lead_snps", "snps", "cs_sizes", "parameters", "purity"}
        assert set(d.keys()) == expected_keys

    def test_roundtrip(self):
        cs = _make_credset(purity=[0.9, 0.7])
        d = cs.to_dict()
        restored = CredibleSet.from_dict(d, cs.pips)
        assert restored.tool == cs.tool
        assert restored.n_cs == cs.n_cs
        assert restored.coverage == cs.coverage
        assert restored.lead_snps == cs.lead_snps
        assert restored.snps == cs.snps
        assert restored.purity == cs.purity

    def test_to_dict_without_purity(self):
        cs = _make_credset()
        d = cs.to_dict()
        assert d["purity"] is None


# ---------------------------------------------------------------------------
# TestFilterCredsetByPurity
# ---------------------------------------------------------------------------
class TestFilterCredsetByPurity:
    """Tests for filter_credset_by_purity function."""

    def test_purity_none_returns_original(self):
        cs = _make_credset()
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == cs.n_cs  # no filtering applied

    def test_min_purity_zero_returns_original(self):
        cs = _make_credset(purity=[0.3, 0.8])
        result = filter_credset_by_purity(cs, min_purity=0.0)
        assert result.n_cs == cs.n_cs

    def test_filter_low_purity(self):
        cs = _make_credset(purity=[0.3, 0.8])
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == 1
        assert result.purity == [0.8]

    def test_all_filtered_returns_empty(self):
        cs = _make_credset(purity=[0.1, 0.2])
        result = filter_credset_by_purity(cs, min_purity=0.9)
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []

    def test_empty_purity_list_returns_original(self):
        cs = CredibleSet(
            tool="test",
            parameters={},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
            purity=[],
        )
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == 0


# ---------------------------------------------------------------------------
# TestGenerateCsSummary
# ---------------------------------------------------------------------------
class TestGenerateCsSummary:
    """Tests for generate_cs_summary function."""

    def _make_locus_set(self):
        """Create a simple LocusSet for testing."""
        n_snps = 5
        snpids = [f"1-{1000 + i * 100}-A-G" for i in range(n_snps)]
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n_snps,
                ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                ColName.EA: ["A"] * n_snps,
                ColName.NEA: ["G"] * n_snps,
                ColName.EAF: [0.3] * n_snps,
                ColName.BETA: [0.1] * n_snps,
                ColName.SE: [0.01] * n_snps,
                ColName.P: [1e-8] * n_snps,
            }
        )
        r = np.eye(n_snps)
        ld_map = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n_snps,
                ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                ColName.A1: ["A"] * n_snps,
                ColName.A2: ["G"] * n_snps,
            }
        )
        ld = LDMatrix(ld_map, r)
        locus = Locus("EUR", "test", 10000, sumstats, 1000, 2000, ld=ld)
        return LocusSet([locus])

    def test_empty_dataframe_returns_empty(self):
        ls = self._make_locus_set()
        result = generate_cs_summary(pd.DataFrame(), "locus_1", ls)
        assert result == []

    def test_single_cs(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1, 1, 1],
                "PIP": [0.8, 0.15, 0.05],
                "SNPID": ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        assert len(result) == 1
        row = result[0]
        assert row["locus_id"] == "locus_1"
        assert row["cs_id"] == 1
        assert row["cs_size"] == 3
        assert row["lead_snp"] == "1-1000-A-G"

    def test_multiple_cs(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1, 1, 2, 2],
                "PIP": [0.8, 0.15, 0.7, 0.3],
                "SNPID": ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G", "1-1300-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        assert len(result) == 2
        assert result[0]["cs_id"] == 1
        assert result[1]["cs_id"] == 2

    def test_pip_thresholds(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1, 1, 1, 1],
                "PIP": [0.95, 0.55, 0.15, 0.05],
                "SNPID": ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G", "1-1300-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        row = result[0]
        assert row["pip_01"] == 3  # >=0.1: 0.95, 0.55, 0.15
        assert row["pip_05"] == 2  # >=0.5: 0.95, 0.55
        assert row["pip_09"] == 1  # >=0.9: 0.95

    def test_required_fields(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1],
                "PIP": [0.8],
                "SNPID": ["1-1000-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        expected_keys = {"locus_id", "cs_id", "lead_snp", "cs_size", "pip_01", "pip_05", "pip_09", "purity"}
        assert set(result[0].keys()) == expected_keys
