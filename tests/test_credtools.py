#!/usr/bin/env python
"""Tests for `credtools` package."""

import os
import warnings
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from typing import List, Optional

from credtools import __version__
from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.credtools import (
    _adaptive_fine_map,
    _adaptive_fine_map_multi,
    _empty_credible_set,
    _generate_run_summary,
    _is_success,
    fine_map,
)
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet, load_locus


def test_version():
    """Test that version is a string."""
    assert isinstance(__version__, str)


@pytest.fixture
def sample_locus():
    """Create a sample locus for testing."""
    # Create sample data with all mandatory columns
    sumstats = pd.DataFrame(
        {
            ColName.SNPID: ["rs1", "rs2", "rs3"],
            ColName.CHR: [1, 1, 1],
            ColName.BP: [1000, 2000, 3000],
            ColName.EA: ["A", "C", "G"],  # Effect allele
            ColName.NEA: ["T", "G", "T"],  # Non-effect allele
            ColName.EAF: [0.1, 0.2, 0.3],  # Effect allele frequency
            ColName.A1: ["A", "C", "G"],
            ColName.A2: ["T", "G", "T"],
            ColName.BETA: [0.1, 0.2, 0.3],
            ColName.SE: [0.01, 0.02, 0.03],
            ColName.P: [0.001, 0.002, 0.003],
            ColName.MAF: [0.1, 0.2, 0.3],
        }
    )

    # Create sample LD matrix
    ld_matrix = np.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.3], [0.2, 0.3, 1.0]])

    # Create LDMatrix object
    ld = LDMatrix(sumstats, ld_matrix)

    return Locus(
        popu="EUR",
        cohort="test",
        sample_size=1000,
        sumstats=sumstats,
        locus_start=1000,
        locus_end=3000,
        ld=ld,
    )


def test_locus_creation(sample_locus):
    """Test that a Locus object can be created."""
    assert isinstance(sample_locus, Locus)
    assert len(sample_locus.sumstats) == 3
    assert sample_locus.ld.r.shape == (3, 3)


def test_credible_set_creation(sample_locus):
    """Test that a CredibleSet can be created."""
    # Create a sample credible set
    cs = CredibleSet(
        tool=Method.SUSIE,
        parameters={"max_causal": 1, "coverage": 0.95},
        coverage=0.95,
        n_cs=1,
        cs_sizes=[2],
        lead_snps=["rs1"],
        snps=[["rs1", "rs2"]],
        pips=pd.Series({"rs1": 0.8, "rs2": 0.2, "rs3": 0.0}),
    )

    assert isinstance(cs, CredibleSet)
    assert len(cs.snps) == 1
    assert cs.coverage == 0.95


def test_locus_set_creation(sample_locus):
    """Test that a LocusSet can be created."""
    locus_set = LocusSet([sample_locus])
    assert isinstance(locus_set, LocusSet)
    assert locus_set.n_loci == 1


@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory for file operations."""
    return tmp_path


def test_file_io(temp_dir, sample_locus):
    """Test that loci can be saved and loaded."""
    # Save locus
    output_prefix = str(temp_dir / "test_locus")
    sample_locus.sumstats.to_csv(
        f"{output_prefix}.sumstats.gz", sep="\t", index=False, compression="gzip"
    )
    np.savez_compressed(
        f"{output_prefix}.ld.npz", ld=sample_locus.ld.r.astype(np.float16)
    )
    sample_locus.ld.map.to_csv(
        f"{output_prefix}.ldmap.gz", sep="\t", index=False, compression="gzip"
    )

    # Load locus
    loaded_locus = load_locus(
        prefix=output_prefix,
        popu="EUR",
        cohort="test",
        sample_size=1000,
        locus_start=1000,
        locus_end=3000,
    )

    assert isinstance(loaded_locus, Locus)
    assert len(loaded_locus.sumstats) == len(sample_locus.sumstats)
    assert loaded_locus.ld.r.shape == sample_locus.ld.r.shape


def _make_test_locus(
    popu: str,
    cohort: str,
    beta_scale: float,
    p_values: Optional[List[float]] = None,
) -> Locus:
    if p_values is None:
        p_values = [1e-8, 5e-8, 1e-7]
    sumstats = pd.DataFrame(
        {
            ColName.SNPID: ["1-100-A-G", "1-200-A-G", "1-300-A-G"],
            ColName.CHR: [1, 1, 1],
            ColName.BP: [100, 200, 300],
            ColName.RSID: ["rs1", "rs2", "rs3"],
            ColName.EA: ["A", "A", "A"],
            ColName.NEA: ["G", "G", "G"],
            ColName.EAF: [0.2, 0.25, 0.3],
            ColName.MAF: [0.2, 0.25, 0.3],
            ColName.A1: ["A", "A", "A"],
            ColName.A2: ["G", "G", "G"],
            ColName.BETA: [0.2 * beta_scale, 0.15 * beta_scale, 0.1 * beta_scale],
            ColName.SE: [0.05, 0.05, 0.05],
            ColName.P: p_values,
        }
    )
    ld = LDMatrix(sumstats, np.eye(len(sumstats), dtype=float))
    return Locus(
        popu=popu,
        cohort=cohort,
        sample_size=1000,
        sumstats=sumstats,
        locus_start=0,
        locus_end=400,
        ld=ld,
    )


def test_fine_map_embeds_per_dataset_columns():
    locus_primary = _make_test_locus("EUR", "cohort1", beta_scale=1.0)
    locus_secondary = _make_test_locus("AFR", "cohort2", beta_scale=0.8)
    locus_set = LocusSet([locus_primary, locus_secondary])

    combined = fine_map(
        locus_set,
        tool="abf",
        max_causal=1,
        set_L_by_cojo=False,
        combine_cred="union",
        combine_pip="max",
    )

    combined_df = combined.create_enhanced_pips_df(locus_set)
    assert not combined_df.empty
    assert ColName.PIP in combined_df.columns

    expected_columns = {
        "EUR_cohort1_PIP",
        "EUR_cohort1_CRED",
        "AFR_cohort2_PIP",
        "AFR_cohort2_CRED",
    }
    assert expected_columns.issubset(set(combined_df.columns))

    for col in expected_columns:
        if col.endswith("_PIP"):
            assert combined_df[col].ge(0).all()
        else:
            assert combined_df[col].dtype.kind in {"i", "f"}


def test_single_input_returns_zero_without_significant_snp():
    high_pvals = [1e-4, 2e-4, 1e-3]
    locus_primary = _make_test_locus(
        "EUR", "cohort1", beta_scale=1.0, p_values=high_pvals
    )
    locus_secondary = _make_test_locus(
        "AFR", "cohort2", beta_scale=0.8, p_values=high_pvals
    )
    locus_set = LocusSet([locus_primary, locus_secondary])

    prefixes = [
        f"{locus_primary.popu}_{locus_primary.cohort}_",
        f"{locus_secondary.popu}_{locus_secondary.cohort}_",
    ]

    for tool in ("finemap", "susie", "rsparsepro"):
        result = fine_map(
            locus_set,
            tool=tool,
            max_causal=1,
            set_L_by_cojo=False,
            significant_threshold=5e-8,
        )
        assert result.n_cs == 0
        assert (result.pips == 0).all()
        df = result.create_enhanced_pips_df(locus_set)
        assert (df["PIP"] == 0).all()
        assert (df["CRED"] == 0).all()
        for prefix in prefixes:
            pip_col = f"{prefix}PIP"
            cred_col = f"{prefix}CRED"
            assert pip_col in df.columns
            assert cred_col in df.columns
            assert (df[pip_col] == 0).all()
            assert (df[cred_col] == 0).all()


# ---------------------------------------------------------------------------
# TestIsSuccess
# ---------------------------------------------------------------------------
class TestIsSuccess:
    """Tests for _is_success function."""

    @pytest.mark.parametrize(
        "n_cs, max_causal, expected",
        [
            (0, 5, False),  # 0 CS is failure
            (3, 5, True),  # 3 < 5 is success
            (5, 5, False),  # n_cs == max_causal is failure (saturated)
            (6, 5, False),  # n_cs > max_causal is failure
            (1, 2, True),  # 1 < 2 is success
        ],
    )
    def test_parametrized(self, n_cs, max_causal, expected):
        cs = CredibleSet(
            tool="test",
            parameters={},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[1] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}"] for i in range(n_cs)],
            pips=pd.Series(dtype=float),
        )
        assert _is_success(cs, max_causal) == expected


# ---------------------------------------------------------------------------
# TestEmptyCredibleSet
# ---------------------------------------------------------------------------
class TestEmptyCredibleSet:
    """Tests for _empty_credible_set function."""

    def test_returns_zero_cs(self):
        result = _empty_credible_set("susie")
        assert result.n_cs == 0

    def test_tool_name_correct(self):
        result = _empty_credible_set("finemap")
        assert result.tool == "finemap"

    def test_adaptive_failed_flag(self):
        result = _empty_credible_set("susie")
        assert result.parameters.get("adaptive_failed") is True

    def test_empty_pips(self):
        result = _empty_credible_set("susie")
        assert len(result.pips) == 0


# ---------------------------------------------------------------------------
# TestAdaptiveFinemap
# ---------------------------------------------------------------------------
class TestAdaptiveFinemap:
    """Tests for _adaptive_fine_map function."""

    def _make_cs(self, n_cs, tool="susie"):
        return CredibleSet(
            tool=tool,
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}", f"s{i}_b"] for i in range(n_cs)],
            pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
        )

    def test_success_on_first_try(self):
        """Tool returns n_cs < max_causal on first try."""
        tool_func = MagicMock(return_value=self._make_cs(3))
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 5, tool_func, {})
        assert result.n_cs == 3
        assert tool_func.call_count == 1

    def test_saturated_then_increase(self):
        """n_cs == max_causal triggers increase."""
        # First call: n_cs=5 (saturated at max_causal=5)
        # Second call with max_causal=10: n_cs=7 (success)
        call_count = 0

        def side_effect(locus, max_causal=5, **kwargs):
            nonlocal call_count
            call_count += 1
            if max_causal == 5:
                return self._make_cs(5)
            else:
                return self._make_cs(7)

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 5, tool_func, {})
        assert result.n_cs == 7

    def test_initial_failure_then_decrease(self):
        """Initial max_causal fails, then decreasing works."""

        def side_effect(locus, max_causal=5, **kwargs):
            if max_causal >= 4:
                raise RuntimeError("convergence failed")
            return self._make_cs(2)

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 5, tool_func, {})
        assert result.n_cs == 2

    def test_all_attempts_fail(self):
        """All attempts fail, return empty."""
        tool_func = MagicMock(side_effect=RuntimeError("fail"))
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 3, tool_func, {})
        assert result.n_cs == 0
        assert result.parameters.get("adaptive_failed") is True


# ---------------------------------------------------------------------------
# TestAdaptiveFineMapMulti
# ---------------------------------------------------------------------------
class TestAdaptiveFineMapMulti:
    """Tests for _adaptive_fine_map_multi function."""

    def _make_cs(self, n_cs, tool="multisusie"):
        return CredibleSet(
            tool=tool,
            parameters={},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}"] for i in range(n_cs)],
            pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
        )

    def test_success_on_first_try(self):
        tool_func = MagicMock(return_value=self._make_cs(3))
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 5, tool_func, {})
        assert result.n_cs == 3

    def test_all_fail_returns_empty(self):
        tool_func = MagicMock(side_effect=RuntimeError("fail"))
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 3, tool_func, {})
        assert result.n_cs == 0


# ---------------------------------------------------------------------------
# TestFineMapBranches
# ---------------------------------------------------------------------------
class TestFineMapBranches:
    """Tests for fine_map function branches."""

    def test_unknown_tool_raises(self):
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        with pytest.raises(ValueError, match="not recognized"):
            fine_map(locus_set, tool="unknown_tool", set_L_by_cojo=False)

    def test_negative_timeout_raises(self):
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        with pytest.raises(ValueError, match="timeout_minutes must be a positive"):
            fine_map(
                locus_set,
                tool="finemap",
                timeout_minutes=-5,
                set_L_by_cojo=False,
            )

    def test_strategy_deprecation_warning(self):
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fine_map(
                locus_set,
                tool="abf",
                strategy="independent",
                max_causal=1,
                set_L_by_cojo=False,
            )
            assert len(w) >= 1
            assert any(issubclass(x.category, DeprecationWarning) for x in w)

    def test_single_locus_abf(self):
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        result = fine_map(locus_set, tool="abf", max_causal=1, set_L_by_cojo=False)
        assert isinstance(result, CredibleSet)

    def test_purity_filtering(self):
        """Purity filter parameter should be passed through."""
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        result = fine_map(
            locus_set, tool="abf", max_causal=1, set_L_by_cojo=False, purity=0.99
        )
        # With purity=0.99 and identity LD, all CS should survive
        assert isinstance(result, CredibleSet)


# ---------------------------------------------------------------------------
# TestGenerateRunSummary
# ---------------------------------------------------------------------------
class TestGenerateRunSummary:
    """Tests for _generate_run_summary function."""

    def test_writes_log_file(self, tmp_path):
        summary = {
            "start_time": "2024-01-01T00:00:00",
            "end_time": "2024-01-01T00:01:00",
            "total_loci": 10,
            "successful_loci": 8,
            "failed_loci": 2,
            "errors": ["Error 1", "Error 2"],
            "tool": "susie",
            "meta_method": "meta_all",
            "parameters": {"max_causal": 5},
        }
        log_file = str(tmp_path / "run_summary.log")
        _generate_run_summary(summary, log_file)
        assert os.path.exists(log_file)

    def test_content_correct(self, tmp_path):
        summary = {
            "start_time": "2024-01-01T00:00:00",
            "end_time": "2024-01-01T00:01:00",
            "total_loci": 5,
            "successful_loci": 5,
            "failed_loci": 0,
            "errors": [],
            "tool": "finemap",
            "meta_method": "no_meta",
            "parameters": {"coverage": 0.95},
        }
        log_file = str(tmp_path / "run_summary.log")
        _generate_run_summary(summary, log_file)
        with open(log_file) as f:
            content = f.read()
        assert "CREDTOOLS FINE-MAPPING RUN SUMMARY" in content
        assert "Total Loci: 5" in content
        assert "Successful: 5" in content
        assert "finemap" in content

    def test_errors_section(self, tmp_path):
        summary = {
            "start_time": "t1",
            "end_time": "t2",
            "total_loci": 1,
            "successful_loci": 0,
            "failed_loci": 1,
            "errors": ["Something broke"],
            "tool": "susie",
            "meta_method": "meta_all",
            "parameters": {},
        }
        log_file = str(tmp_path / "run_summary.log")
        _generate_run_summary(summary, log_file)
        with open(log_file) as f:
            content = f.read()
        assert "Error Details:" in content
        assert "Something broke" in content
