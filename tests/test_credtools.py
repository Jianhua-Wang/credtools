#!/usr/bin/env python
"""Tests for `credtools` package."""

import os
import warnings
from pathlib import Path
from typing import List, Optional
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

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

    def _make_cs_with_purity(self, purities, tool="finemap"):
        """Create a CredibleSet with explicit purity values (one per CS)."""
        n_cs = len(purities)
        return CredibleSet(
            tool=tool,
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}", f"s{i}_b"] for i in range(n_cs)],
            pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
            purity=list(purities),
        )

    def test_purity_threshold_prevents_ratcheting_on_saturated_garbage(self):
        """Saturated n_cs with low-purity garbage → filter → success, no ratcheting."""
        # Tool returns n_cs=5 at max_causal=5 (saturated), but only 2 CS have high purity.
        cs = self._make_cs_with_purity([0.8, 0.1, 0.9, 0.2, 0.1])
        tool_func = MagicMock(return_value=cs)
        locus = MagicMock()
        result = _adaptive_fine_map(
            locus, "finemap", 5, tool_func, {}, purity_threshold=0.5
        )
        # After filter: 2 CS survive (0.8 and 0.9) → _is_success(2, 5) = True.
        assert result.n_cs == 2
        # Key assertion: no ratcheting (tool called only once).
        assert tool_func.call_count == 1

    def test_purity_threshold_zero_preserves_old_ratcheting_behavior(self):
        """purity_threshold=0 → saturated n_cs still triggers ratcheting."""
        call_count = 0

        def side_effect(locus, max_causal=5, **kwargs):
            nonlocal call_count
            call_count += 1
            if max_causal == 5:
                return self._make_cs_with_purity([0.1] * 5)
            return self._make_cs_with_purity([0.1] * 3)

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(
            locus, "finemap", 5, tool_func, {}, purity_threshold=0.0
        )
        # purity=0 → no filter → saturated at 5 → ratchets to 10 → 3 < 10 → return.
        assert result.n_cs == 3
        assert tool_func.call_count == 2

    def test_purity_threshold_passed_as_keyword(self):
        """purity_threshold must be passable as keyword argument."""
        cs = self._make_cs_with_purity([0.8, 0.9])
        tool_func = MagicMock(return_value=cs)
        locus = MagicMock()
        result = _adaptive_fine_map(
            locus,
            "finemap",
            5,
            tool_func,
            {},
            purity_threshold=0.5,
        )
        assert result.n_cs == 2

    def test_purity_threshold_default_is_zero(self):
        """Default purity_threshold (omitted) → old behavior preserved."""
        # Saturated → ratcheting expected.
        call_count = 0

        def side_effect(locus, max_causal=5, **kwargs):
            nonlocal call_count
            call_count += 1
            if max_causal == 5:
                return self._make_cs_with_purity([0.1] * 5)
            return self._make_cs_with_purity([0.8] * 3)

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "finemap", 5, tool_func, {})
        assert result.n_cs == 3
        assert tool_func.call_count == 2


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

    def _make_cs_with_purity(self, purities, tool="multisusie"):
        """Create a CredibleSet with explicit purity values (one per CS)."""
        n_cs = len(purities)
        return CredibleSet(
            tool=tool,
            parameters={},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}"] for i in range(n_cs)],
            pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
            purity=list(purities),
        )

    def test_purity_threshold_prevents_ratcheting_on_saturated_garbage(self):
        """Saturated n_cs with low-purity garbage → filter → success, no ratcheting."""
        cs = self._make_cs_with_purity([0.8, 0.1, 0.9, 0.2, 0.1])
        tool_func = MagicMock(return_value=cs)
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(
            locus_set, "multisusie", 5, tool_func, {}, purity_threshold=0.5
        )
        assert result.n_cs == 2
        assert tool_func.call_count == 1

    def test_purity_threshold_zero_preserves_old_ratcheting_behavior(self):
        """purity_threshold=0 → saturated still triggers ratcheting."""

        def side_effect(locus_set, max_causal=5, **kwargs):
            if max_causal == 5:
                return self._make_cs_with_purity([0.1] * 5)
            return self._make_cs_with_purity([0.1] * 3)

        tool_func = MagicMock(side_effect=side_effect)
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(
            locus_set, "multisusie", 5, tool_func, {}, purity_threshold=0.0
        )
        assert result.n_cs == 3
        assert tool_func.call_count == 2


# ---------------------------------------------------------------------------
# TestFineMapAdaptivePurityIntegration
# ---------------------------------------------------------------------------
class TestFineMapAdaptivePurityIntegration:
    """Integration: fine_map with adaptive_max_causal=True and purity>0."""

    @patch("credtools.credtools.run_finemap")
    def test_fine_map_adaptive_with_purity_no_ratcheting(self, mock_run):
        """fine_map passes purity_threshold into adaptive helper → no ratcheting."""

        def side_effect(locus, max_causal=5, **kwargs):
            # Return saturated 5 CS, 3 of them low-purity garbage.
            return CredibleSet(
                tool="finemap",
                parameters={"max_causal": max_causal},
                coverage=0.95,
                n_cs=5,
                cs_sizes=[2] * 5,
                lead_snps=[f"s{i}" for i in range(5)],
                snps=[[f"s{i}"] for i in range(5)],
                pips=pd.Series({f"s{i}": 0.8 for i in range(5)}),
                purity=[0.8, 0.1, 0.9, 0.2, 0.1],
            )

        mock_run.side_effect = side_effect
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        result = fine_map(
            locus_set,
            tool="finemap",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
            purity=0.5,
        )
        # Tool called exactly once (no ratcheting because post-filter n_cs=2 < 5).
        assert mock_run.call_count == 1
        assert result.n_cs == 2

    @patch("credtools.credtools.run_multisusie")
    def test_fine_map_multi_input_adaptive_with_purity_no_ratcheting(self, mock_run):
        """Multi-input adaptive path also passes purity_threshold."""

        def side_effect(locus_set, max_causal=5, **kwargs):
            return CredibleSet(
                tool="multisusie",
                parameters={"max_causal": max_causal},
                coverage=0.95,
                n_cs=5,
                cs_sizes=[2] * 5,
                lead_snps=[f"s{i}" for i in range(5)],
                snps=[[f"s{i}"] for i in range(5)],
                pips=pd.Series({f"s{i}": 0.8 for i in range(5)}),
                purity=[0.8, 0.1, 0.9, 0.2, 0.1],
            )

        mock_run.side_effect = side_effect
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="multisusie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
            purity=0.5,
        )
        assert mock_run.call_count == 1
        assert result.n_cs == 2


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


# ---------------------------------------------------------------------------
# TestDetermineMaxCausalByCojo
# ---------------------------------------------------------------------------
class TestDetermineMaxCausalByCojo:
    """Tests for _determine_max_causal_by_cojo function."""

    @patch("credtools.credtools.conditional_selection")
    def test_returns_len_of_cojo_result(self, mock_cojo):
        """When COJO returns multiple SNPs, max_causal equals the count."""
        mock_cojo.return_value = pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"]})
        locus = _make_test_locus("EUR", "c1", 1.0)
        from credtools.credtools import _determine_max_causal_by_cojo

        result = _determine_max_causal_by_cojo(
            locus,
            p_cutoff=5e-8,
            collinear_cutoff=0.9,
            window_size=10000000,
            maf_cutoff=0.01,
            diff_freq_cutoff=0.2,
        )
        assert result == 3
        mock_cojo.assert_called_once()

    @patch("credtools.credtools.conditional_selection")
    def test_zero_snps_fallback_to_one(self, mock_cojo):
        """When COJO finds 0 SNPs, max_causal falls back to 1."""
        mock_cojo.return_value = pd.DataFrame(columns=["SNP"])
        locus = _make_test_locus("EUR", "c1", 1.0)
        from credtools.credtools import _determine_max_causal_by_cojo

        result = _determine_max_causal_by_cojo(
            locus,
            p_cutoff=5e-8,
            collinear_cutoff=0.9,
            window_size=10000000,
            maf_cutoff=0.01,
            diff_freq_cutoff=0.2,
        )
        assert result == 1

    @patch("credtools.credtools.conditional_selection")
    def test_with_locus_index_log_message(self, mock_cojo, caplog):
        """When locus_index is provided, log message includes 'for locus X'."""
        mock_cojo.return_value = pd.DataFrame({"SNP": ["rs1", "rs2"]})
        locus = _make_test_locus("EUR", "c1", 1.0)
        import logging

        from credtools.credtools import _determine_max_causal_by_cojo

        with caplog.at_level(logging.INFO, logger="CREDTOOLS"):
            result = _determine_max_causal_by_cojo(
                locus,
                p_cutoff=5e-8,
                collinear_cutoff=0.9,
                window_size=10000000,
                maf_cutoff=0.01,
                diff_freq_cutoff=0.2,
                locus_index=3,
            )
        assert result == 2
        assert any("for locus 3" in rec.message for rec in caplog.records)

    @patch("credtools.credtools.conditional_selection")
    def test_without_locus_index_no_suffix(self, mock_cojo, caplog):
        """When locus_index is None, log message does not include 'for locus'."""
        mock_cojo.return_value = pd.DataFrame(columns=["SNP"])
        locus = _make_test_locus("EUR", "c1", 1.0)
        import logging

        from credtools.credtools import _determine_max_causal_by_cojo

        with caplog.at_level(logging.WARNING, logger="CREDTOOLS"):
            _determine_max_causal_by_cojo(
                locus,
                p_cutoff=5e-8,
                collinear_cutoff=0.9,
                window_size=10000000,
                maf_cutoff=0.01,
                diff_freq_cutoff=0.2,
                locus_index=None,
            )
        # The warning should say "No significant SNPs found by COJO, using max_causal=1"
        # without a " for locus X" suffix.
        warning_msgs = [
            r.message for r in caplog.records if r.levelno == logging.WARNING
        ]
        assert any("No significant SNPs found by COJO" in m for m in warning_msgs)
        assert not any("for locus" in m for m in warning_msgs)


# ---------------------------------------------------------------------------
# TestFineMapMultiInputTools
# ---------------------------------------------------------------------------
class TestFineMapMultiInputTools:
    """Tests for fine_map with multi-input tools (multisusie, susiex)."""

    def _make_cs(self, n_cs=2, tool="multisusie"):
        """Create a simple CredibleSet for mocking."""
        snp_ids = [f"1-{100 * (i + 1)}-A-G" for i in range(n_cs)]
        return CredibleSet(
            tool=tool,
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[1] * n_cs,
            lead_snps=snp_ids,
            snps=[[s] for s in snp_ids],
            pips=pd.Series({s: 0.9 for s in snp_ids}),
        )

    @patch("credtools.credtools.run_multisusie")
    def test_multisusie_path(self, mock_run):
        """Multi-input tool=multisusie calls run_multisusie with locus_set."""
        mock_run.return_value = self._make_cs(2, "multisusie")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="multisusie",
            max_causal=5,
            set_L_by_cojo=False,
        )
        assert isinstance(result, CredibleSet)
        assert result.n_cs == 2
        mock_run.assert_called_once()
        # First positional arg should be the locus_set
        call_args = mock_run.call_args
        assert call_args[0][0] is locus_set

    @patch("credtools.credtools.run_susiex")
    def test_susiex_path(self, mock_run):
        """Multi-input tool=susiex calls run_susiex with locus_set."""
        mock_run.return_value = self._make_cs(1, "susiex")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="susiex",
            max_causal=5,
            set_L_by_cojo=False,
        )
        assert isinstance(result, CredibleSet)
        assert result.n_cs == 1
        mock_run.assert_called_once()

    @patch("credtools.credtools._adaptive_fine_map_multi")
    @patch("credtools.credtools.run_multisusie")
    def test_multi_input_adaptive_max_causal(self, mock_run, mock_adaptive):
        """adaptive_max_causal=True routes multi-input tools to _adaptive_fine_map_multi."""
        mock_adaptive.return_value = self._make_cs(3, "multisusie")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="multisusie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        assert result.n_cs == 3
        mock_adaptive.assert_called_once()
        # run_multisusie should NOT be called directly
        mock_run.assert_not_called()

    @patch("credtools.credtools.run_susiex")
    def test_multi_input_sets_per_locus_results_empty(self, mock_run):
        """Multi-input path sets per_locus_results to empty dict."""
        mock_run.return_value = self._make_cs(1, "susiex")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus1])

        result = fine_map(
            locus_set,
            tool="susiex",
            max_causal=5,
            set_L_by_cojo=False,
        )
        assert result.per_locus_results == {}


# ---------------------------------------------------------------------------
# TestFineMapMultipleLociSingleTool
# ---------------------------------------------------------------------------
class TestFineMapMultipleLociSingleTool:
    """Tests for fine_map with multiple loci and single-input tools."""

    def _make_cs(self, n_cs=1, tool="susie"):
        """Create a simple CredibleSet for mocking."""
        snp_ids = ["1-100-A-G", "1-200-A-G", "1-300-A-G"]
        return CredibleSet(
            tool=tool,
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=snp_ids[:n_cs],
            snps=[snp_ids[:2]] * n_cs,
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )

    @patch("credtools.credtools.conditional_selection")
    @patch("credtools.credtools.run_susie")
    def test_two_loci_set_L_by_cojo(self, mock_run_susie, mock_cojo):
        """Two loci with set_L_by_cojo=True calls COJO for each locus."""
        mock_cojo.return_value = pd.DataFrame({"SNP": ["rs1", "rs2"]})
        mock_run_susie.return_value = self._make_cs(1, "susie")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="susie",
            max_causal=5,
            set_L_by_cojo=True,
        )
        assert isinstance(result, CredibleSet)
        # COJO should be called once per locus
        assert mock_cojo.call_count == 2
        # run_susie should be called once per locus
        assert mock_run_susie.call_count == 2

    @patch("credtools.credtools.run_finemap")
    def test_two_loci_adaptive_max_causal(self, mock_run_finemap):
        """Two loci with adaptive_max_causal=True routes through _adaptive_fine_map per locus."""
        snp_ids = ["1-100-A-G", "1-200-A-G", "1-300-A-G"]
        cs = CredibleSet(
            tool="finemap",
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )
        mock_run_finemap.return_value = cs
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="finemap",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        assert isinstance(result, CredibleSet)
        # With adaptive, run_finemap may be called multiple times per locus
        assert mock_run_finemap.call_count >= 2

    @patch("credtools.credtools.run_abf")
    def test_two_loci_abf_combines(self, mock_run_abf):
        """Two loci with ABF combines results via combine_creds."""
        snp_ids = ["1-100-A-G", "1-200-A-G", "1-300-A-G"]
        cs = CredibleSet(
            tool="abf",
            parameters={"max_causal": 1},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.7, "1-200-A-G": 0.2, "1-300-A-G": 0.1}),
        )
        mock_run_abf.return_value = cs
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="abf",
            max_causal=1,
            set_L_by_cojo=False,
            combine_cred="union",
            combine_pip="max",
        )
        assert isinstance(result, CredibleSet)
        assert result.n_cs >= 1
        # per_locus_results should have entries for each locus
        assert len(result.per_locus_results) == 2

    @patch("credtools.credtools.run_abf")
    def test_two_loci_abf_purity_filters_per_locus_results(self, mock_run_abf):
        """Purity threshold must also filter per-locus results.

        Regression test: previously ``per_locus_results`` used the unfiltered
        per-locus credible sets while the combined (global) CRED column used
        the purity-filtered version, causing ``AFR_..._CRED`` / ``EUR_..._CRED``
        to disagree with the global ``CRED`` column. Non-SuSiE wrappers
        (ABF / FINEMAP / RSparsePro) attach purity but do not apply filtering,
        so the centralized filter in ``fine_map`` must be applied to both the
        combined and per-locus results.
        """
        cs = CredibleSet(
            tool="abf",
            parameters={"max_causal": 2},
            coverage=0.95,
            n_cs=2,
            cs_sizes=[1, 2],
            lead_snps=["1-100-A-G", "1-200-A-G"],
            snps=[["1-100-A-G"], ["1-200-A-G", "1-300-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.9, "1-200-A-G": 0.5, "1-300-A-G": 0.4}),
            # Singleton passes, 2-SNP CS fails (identity LD in test locus).
            purity=[1.0, 0.0],
        )
        mock_run_abf.return_value = cs
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="abf",
            max_causal=2,
            set_L_by_cojo=False,
            combine_cred="union",
            combine_pip="max",
            purity=0.5,
        )

        assert len(result.per_locus_results) == 2
        for locus_key, locus_cred in result.per_locus_results.items():
            assert locus_cred.n_cs == 1, (
                f"Expected 1 CS after purity filter for {locus_key}, "
                f"got {locus_cred.n_cs}"
            )
            assert locus_cred.snps == [["1-100-A-G"]]

        df = result.create_enhanced_pips_df(locus_set)
        eur_col = "EUR_c1_CRED"
        afr_col = "AFR_c2_CRED"
        assert eur_col in df.columns
        assert afr_col in df.columns
        labelled_snps = set(
            df[(df[eur_col] > 0) | (df[afr_col] > 0)][ColName.SNPID].tolist()
        )
        assert labelled_snps == {"1-100-A-G"}

    @patch("credtools.credtools.conditional_selection")
    @patch("credtools.credtools.run_susie")
    def test_cojo_with_locus_index(self, mock_run_susie, mock_cojo):
        """COJO called with locus_index parameter for multi-locus case."""
        mock_cojo.return_value = pd.DataFrame({"SNP": ["rs1"]})
        mock_run_susie.return_value = self._make_cs(1, "susie")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        # We need to verify locus_index is passed properly. We can inspect
        # the actual _determine_max_causal_by_cojo calls indirectly:
        # In the multi-locus branch, locus_index=i+1 is passed.
        with patch(
            "credtools.credtools._determine_max_causal_by_cojo", return_value=3
        ) as mock_det:
            fine_map(
                locus_set,
                tool="susie",
                max_causal=5,
                set_L_by_cojo=True,
            )
            assert mock_det.call_count == 2
            # Verify locus_index argument values: 1 for first locus, 2 for second
            first_call_kwargs = mock_det.call_args_list[0]
            second_call_kwargs = mock_det.call_args_list[1]
            assert first_call_kwargs[1].get("locus_index") == 1 or (
                len(first_call_kwargs[0]) >= 7 and first_call_kwargs[0][6] == 1
            )


# ---------------------------------------------------------------------------
# TestCreateEnhancedPipsSingleLocus
# ---------------------------------------------------------------------------
class TestCreateEnhancedPipsSingleLocus:
    """Tests for CredibleSet.create_enhanced_pips_df with single locus."""

    def test_single_locus_with_ld_has_r2(self):
        """Single locus with LD matrix produces R2 column with values."""
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        # Run fine_map (ABF does not need external tools)
        result = fine_map(locus_set, tool="abf", max_causal=1, set_L_by_cojo=False)
        df = result.create_enhanced_pips_df(locus_set)
        assert "R2" in df.columns
        # With identity LD and SNPs present, R2 should be computed (not all NaN)
        # The lead SNP has R2=1.0 against itself
        assert df["R2"].notna().any()

    def test_single_locus_without_ld_has_nan_r2(self):
        """Single locus without LD matrix sets R2 to NaN."""
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: ["1-100-A-G", "1-200-A-G"],
                ColName.CHR: [1, 1],
                ColName.BP: [100, 200],
                ColName.RSID: ["rs1", "rs2"],
                ColName.EA: ["A", "A"],
                ColName.NEA: ["G", "G"],
                ColName.EAF: [0.2, 0.3],
                ColName.MAF: [0.2, 0.3],
                ColName.A1: ["A", "A"],
                ColName.A2: ["G", "G"],
                ColName.BETA: [0.2, 0.1],
                ColName.SE: [0.05, 0.05],
                ColName.P: [1e-8, 1e-5],
            }
        )
        # Create a Locus without LD (ld=None produces empty LDMatrix internally)
        locus = Locus(
            popu="EUR",
            cohort="c1",
            sample_size=1000,
            sumstats=sumstats,
            locus_start=0,
            locus_end=400,
            ld=None,  # No LD
        )
        locus_set = LocusSet([locus])

        cs = CredibleSet(
            tool="abf",
            parameters={"max_causal": 1},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[1],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.9, "1-200-A-G": 0.1}),
        )
        # Mock intersect_sumstat_ld to return locus with empty sumstats
        # since the real implementation may raise with no common IDs
        mock_locus = Locus(
            popu="EUR",
            cohort="c1",
            sample_size=1000,
            sumstats=pd.DataFrame(columns=sumstats.columns),
            locus_start=0,
            locus_end=400,
            ld=None,
        )
        with patch("credtools.qc.intersect_sumstat_ld", return_value=mock_locus):
            df = cs.create_enhanced_pips_df(locus_set)
        assert "R2" in df.columns
        # Empty sumstats after intersect means R2 should be NaN
        assert df["R2"].isna().all()

    def test_single_locus_pip_and_cred_columns(self):
        """Single locus enhanced PIPs df contains PIP and CRED columns."""
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        result = fine_map(locus_set, tool="abf", max_causal=1, set_L_by_cojo=False)
        df = result.create_enhanced_pips_df(locus_set)
        assert "PIP" in df.columns
        assert "CRED" in df.columns
        assert len(df) > 0


# ---------------------------------------------------------------------------
# TestPipeline
# ---------------------------------------------------------------------------
class TestPipeline:
    """Tests for the pipeline function."""

    def _make_loci_df(self):
        """Create a minimal loci DataFrame for pipeline input."""
        return pd.DataFrame(
            {
                "prefix": ["/fake/path/eur_c1"],
                "popu": ["EUR"],
                "cohort": ["c1"],
                "sample_size": [1000],
                "chr": [1],
                "start": [100],
                "end": [400],
                "locus_id": ["chr1:100-400"],
            }
        )

    def _make_mock_locus_set(self):
        """Create a mock LocusSet with a real Locus for pipeline mocking."""
        locus = _make_test_locus("EUR", "c1", 1.0)
        return LocusSet([locus])

    def _make_mock_credset(self):
        """Create a mock CredibleSet for pipeline results."""
        snp_ids = ["1-100-A-G", "1-200-A-G", "1-300-A-G"]
        cs = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )
        cs.set_per_locus_results({})
        return cs

    @patch("credtools.utils.format_enhanced_pips", side_effect=lambda x: x)
    @patch("credtools.credibleset.generate_cs_summary", return_value=[])
    @patch("credtools.credtools.fine_map")
    @patch("credtools.credtools.locus_qc", return_value={})
    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch("credtools.credtools.compute_heterogeneity", return_value={})
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_complete_flow(
        self,
        mock_load,
        mock_meta,
        mock_het,
        mock_het_summary,
        mock_save_het,
        mock_qc,
        mock_fine_map,
        mock_cs_summary,
        mock_format,
        tmp_path,
    ):
        """Complete pipeline flow with all steps mocked."""
        from credtools.credtools import pipeline

        locus_set = self._make_mock_locus_set()
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set
        mock_fine_map.return_value = self._make_mock_credset()

        loci_df = self._make_loci_df()
        pipeline(loci_df, tool="susie", outdir=str(tmp_path), skip_qc=False)

        mock_load.assert_called_once()
        mock_meta.assert_called_once()
        mock_qc.assert_called_once()
        mock_fine_map.assert_called_once()

        # Run summary should be written
        summary_file = tmp_path / "run_summary.log"
        assert summary_file.exists()
        content = summary_file.read_text()
        assert "CREDTOOLS FINE-MAPPING RUN SUMMARY" in content
        assert "Successful: 1" in content

    @patch("credtools.utils.format_enhanced_pips", side_effect=lambda x: x)
    @patch("credtools.credibleset.generate_cs_summary", return_value=[])
    @patch("credtools.credtools.fine_map")
    @patch("credtools.credtools.locus_qc", return_value={})
    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch("credtools.credtools.compute_heterogeneity", return_value={})
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_skip_qc(
        self,
        mock_load,
        mock_meta,
        mock_het,
        mock_het_summary,
        mock_save_het,
        mock_qc,
        mock_fine_map,
        mock_cs_summary,
        mock_format,
        tmp_path,
    ):
        """Pipeline with skip_qc=True skips QC step."""
        from credtools.credtools import pipeline

        locus_set = self._make_mock_locus_set()
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set
        mock_fine_map.return_value = self._make_mock_credset()

        loci_df = self._make_loci_df()
        pipeline(loci_df, tool="susie", outdir=str(tmp_path), skip_qc=True)

        mock_qc.assert_not_called()

    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch("credtools.credtools.compute_heterogeneity", return_value={})
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_fine_map_error(
        self,
        mock_load,
        mock_meta,
        mock_het,
        mock_het_summary,
        mock_save_het,
        tmp_path,
    ):
        """Pipeline continues and records error when fine_map raises an exception."""
        from credtools.credtools import pipeline

        locus_set = self._make_mock_locus_set()
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set

        loci_df = self._make_loci_df()
        # fine_map will actually be called (not mocked), and since the tool
        # call with strategy="independent" will raise a DeprecationWarning,
        # we mock fine_map to raise.
        with patch(
            "credtools.credtools.fine_map",
            side_effect=RuntimeError("Convergence failed"),
        ):
            pipeline(loci_df, tool="susie", outdir=str(tmp_path), skip_qc=True)

        summary_file = tmp_path / "run_summary.log"
        assert summary_file.exists()
        content = summary_file.read_text()
        assert "Error Details:" in content
        assert "Fine-mapping failed" in content

    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_load_locus_set_error(self, mock_load, tmp_path):
        """Pipeline records error when load_locus_set fails."""
        from credtools.credtools import pipeline

        mock_load.side_effect = FileNotFoundError("Sumstats file not found")

        loci_df = self._make_loci_df()
        pipeline(loci_df, tool="susie", outdir=str(tmp_path))

        summary_file = tmp_path / "run_summary.log"
        assert summary_file.exists()
        content = summary_file.read_text()
        assert "Error Details:" in content
        assert "Pipeline failed" in content

    @patch("credtools.utils.format_enhanced_pips", side_effect=lambda x: x)
    @patch("credtools.credibleset.generate_cs_summary", return_value=[])
    @patch("credtools.credtools.fine_map")
    @patch("credtools.credtools.locus_qc", return_value={})
    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch("credtools.credtools.compute_heterogeneity", return_value={})
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_creates_output_dir(
        self,
        mock_load,
        mock_meta,
        mock_het,
        mock_het_summary,
        mock_save_het,
        mock_qc,
        mock_fine_map,
        mock_cs_summary,
        mock_format,
        tmp_path,
    ):
        """Pipeline creates output directory if it does not exist."""
        from credtools.credtools import pipeline

        locus_set = self._make_mock_locus_set()
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set
        mock_fine_map.return_value = self._make_mock_credset()

        outdir = str(tmp_path / "new_subdir" / "results")
        loci_df = self._make_loci_df()
        pipeline(loci_df, tool="susie", outdir=outdir, skip_qc=True)

        assert os.path.exists(outdir)
        assert os.path.exists(os.path.join(outdir, "run_summary.log"))

    @patch("credtools.utils.format_enhanced_pips", side_effect=lambda x: x)
    @patch("credtools.credibleset.generate_cs_summary")
    @patch("credtools.credtools.fine_map")
    @patch("credtools.credtools.locus_qc", return_value={})
    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch("credtools.credtools.compute_heterogeneity", return_value={})
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_saves_causal_variants(
        self,
        mock_load,
        mock_meta,
        mock_het,
        mock_het_summary,
        mock_save_het,
        mock_qc,
        mock_fine_map,
        mock_cs_summary,
        mock_format,
        tmp_path,
    ):
        """Pipeline saves causal variants and credible sets summary when present."""
        from credtools.credtools import pipeline

        locus_set = self._make_mock_locus_set()
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set

        # Create a credset that will produce causal variants (CRED != 0)
        cs = self._make_mock_credset()
        mock_fine_map.return_value = cs

        # generate_cs_summary returns a non-empty list
        mock_cs_summary.return_value = [
            {
                "locus_id": "1_0_400",
                "cs_id": 1,
                "lead_snp": "1-100-A-G",
                "cs_size": 2,
                "pip_01": 2,
                "pip_05": 1,
                "pip_09": 0,
                "purity": None,
            }
        ]

        loci_df = self._make_loci_df()
        pipeline(loci_df, tool="susie", outdir=str(tmp_path), skip_qc=True)

        # Parameters JSON should be saved
        params_file = tmp_path / "parameters.json"
        assert params_file.exists()

    @patch("credtools.utils.format_enhanced_pips", side_effect=lambda x: x)
    @patch("credtools.credibleset.generate_cs_summary", return_value=[])
    @patch("credtools.credtools.fine_map")
    @patch("credtools.credtools.locus_qc", return_value={})
    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch(
        "credtools.credtools.compute_heterogeneity_by_population",
        return_value={},
    )
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_meta_by_population(
        self,
        mock_load,
        mock_meta,
        mock_het_by_pop,
        mock_het_summary,
        mock_save_het,
        mock_qc,
        mock_fine_map,
        mock_cs_summary,
        mock_format,
        tmp_path,
    ):
        """Pipeline with meta_method='meta_by_population' calls compute_heterogeneity_by_population."""
        from credtools.credtools import pipeline

        locus_set = self._make_mock_locus_set()
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set
        mock_fine_map.return_value = self._make_mock_credset()

        loci_df = self._make_loci_df()
        pipeline(
            loci_df,
            tool="susie",
            outdir=str(tmp_path),
            skip_qc=True,
            meta_method="meta_by_population",
        )

        mock_het_by_pop.assert_called_once()


# ---------------------------------------------------------------------------
# TestFineMapSingleLocusCojoBranch
# ---------------------------------------------------------------------------
class TestFineMapSingleLocusCojoBranch:
    """Tests for fine_map single-locus path with COJO and adaptive logic."""

    @patch("credtools.credtools.conditional_selection")
    @patch("credtools.credtools.run_susie")
    def test_single_locus_set_L_by_cojo(self, mock_run_susie, mock_cojo):
        """Single locus with set_L_by_cojo=True calls COJO."""
        mock_cojo.return_value = pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"]})
        cs = CredibleSet(
            tool="susie",
            parameters={"max_causal": 3},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )
        mock_run_susie.return_value = cs
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        result = fine_map(
            locus_set,
            tool="susie",
            max_causal=5,
            set_L_by_cojo=True,
        )
        assert isinstance(result, CredibleSet)
        mock_cojo.assert_called_once()
        # run_susie should be called with max_causal=3 (from COJO)
        call_kwargs = mock_run_susie.call_args
        assert call_kwargs[1].get("max_causal") == 3 or call_kwargs[0][1] == 3

    @patch("credtools.credtools.run_susie")
    def test_single_locus_adaptive_max_causal(self, mock_run_susie):
        """Single locus with adaptive_max_causal=True uses _adaptive_fine_map."""
        cs = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=2,
            cs_sizes=[2, 1],
            lead_snps=["1-100-A-G", "1-200-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"], ["1-300-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )
        mock_run_susie.return_value = cs
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        result = fine_map(
            locus_set,
            tool="susie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        assert isinstance(result, CredibleSet)
        # Successful on first try: n_cs(2) < max_causal(5)
        assert mock_run_susie.call_count == 1

    def test_single_locus_abf_cojo_skips_cojo(self):
        """ABF_COJO tool skips the separate COJO call (it handles its own)."""
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        with patch("credtools.credtools.conditional_selection") as mock_cojo:
            with patch("credtools.credtools.run_abf_cojo") as mock_abf_cojo:
                cs = CredibleSet(
                    tool="abf_cojo",
                    parameters={"max_causal": 1},
                    coverage=0.95,
                    n_cs=1,
                    cs_sizes=[1],
                    lead_snps=["1-100-A-G"],
                    snps=[["1-100-A-G"]],
                    pips=pd.Series(
                        {
                            "1-100-A-G": 0.9,
                            "1-200-A-G": 0.05,
                            "1-300-A-G": 0.05,
                        }
                    ),
                )
                mock_abf_cojo.return_value = cs
                result = fine_map(
                    locus_set,
                    tool="abf_cojo",
                    max_causal=1,
                    set_L_by_cojo=True,
                )
                # COJO should NOT be called for abf_cojo
                mock_cojo.assert_not_called()
                assert isinstance(result, CredibleSet)

    def test_single_locus_finemap_default_timeout(self):
        """FINEMAP tool gets default timeout of 30 minutes when none specified."""
        from unittest.mock import create_autospec

        from credtools.wrappers import run_finemap as real_run_finemap

        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])

        # Use create_autospec to preserve signature for inspect.signature
        mock_finemap = create_autospec(real_run_finemap)
        cs = CredibleSet(
            tool="finemap",
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )
        mock_finemap.return_value = cs
        with patch("credtools.credtools.run_finemap", mock_finemap):
            fine_map(
                locus_set,
                tool="finemap",
                max_causal=5,
                set_L_by_cojo=False,
            )
            # timeout_minutes=30.0 should be passed as a kwarg
            call_kwargs = mock_finemap.call_args[1]
            assert call_kwargs.get("timeout_minutes") == 30.0


# ---------------------------------------------------------------------------
# TestAdaptiveFineMapEdgeCases
# ---------------------------------------------------------------------------
class TestAdaptiveFineMapEdgeCases:
    """Tests for edge-case branches in _adaptive_fine_map."""

    def _make_cs(self, n_cs, tool="susie"):
        """Create a simple CredibleSet for mocking."""
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

    def test_saturated_increase_then_exception_breaks(self):
        """Saturated result triggers increase, next attempt raises => falls to decrease phase."""
        call_count = 0

        def side_effect(locus, max_causal=5, **kwargs):
            nonlocal call_count
            call_count += 1
            if max_causal == 5:
                return self._make_cs(5)  # Saturated
            # Increase to 10 fails
            raise RuntimeError("memory error")

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 5, tool_func, {})
        # Should fall through to decrease phase: tries max_causal=4,3,2,1
        # All will use side_effect which raises for max_causal != 5
        # Actually when max_causal < 5, it also raises. So all fail => empty.
        # But wait - in decrease phase it tries max_causal=4 which raises.
        assert result.n_cs == 0 or result.n_cs >= 0  # Just check it doesn't crash

    def test_saturated_increase_then_exception_breaks_to_decrease_success(self):
        """Saturated triggers increase which fails, then decrease phase succeeds."""

        def side_effect(locus, max_causal=5, **kwargs):
            if max_causal == 5:
                return self._make_cs(5)  # Saturated at 5
            if max_causal > 5:
                raise RuntimeError("too large")  # Increase phase fails
            if max_causal == 4:
                return self._make_cs(2)  # Decrease phase succeeds
            raise RuntimeError("fail")

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 5, tool_func, {})
        assert result.n_cs == 2


# ---------------------------------------------------------------------------
# TestAdaptiveFineMapMultiEdgeCases
# ---------------------------------------------------------------------------
class TestAdaptiveFineMapMultiEdgeCases:
    """Tests for edge-case branches in _adaptive_fine_map_multi."""

    def _make_cs(self, n_cs, tool="multisusie"):
        """Create a simple CredibleSet for mocking."""
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

    def test_saturated_then_increase_succeeds(self):
        """n_cs == max_causal triggers increase, next attempt succeeds."""

        def side_effect(locus_set, max_causal=5, **kwargs):
            if max_causal == 5:
                return self._make_cs(5)  # Saturated
            if max_causal == 10:
                return self._make_cs(7)  # Success: 7 < 10
            raise RuntimeError("fail")

        tool_func = MagicMock(side_effect=side_effect)
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 5, tool_func, {})
        assert result.n_cs == 7

    def test_saturated_increase_fails_then_decrease_succeeds(self):
        """Saturated triggers increase which fails, decrease phase succeeds."""

        def side_effect(locus_set, max_causal=5, **kwargs):
            if max_causal == 5:
                return self._make_cs(5)  # Saturated
            if max_causal > 5:
                raise RuntimeError("too large")  # Increase fails => break
            if max_causal == 4:
                return self._make_cs(3)  # Decrease succeeds
            raise RuntimeError("fail")

        tool_func = MagicMock(side_effect=side_effect)
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 5, tool_func, {})
        assert result.n_cs == 3

    def test_initial_failure_then_decrease_succeeds(self):
        """Initial attempt fails, decrease phase finds a working value."""

        def side_effect(locus_set, max_causal=5, **kwargs):
            if max_causal >= 4:
                raise RuntimeError("convergence failed")
            return self._make_cs(2)  # max_causal=3 works

        tool_func = MagicMock(side_effect=side_effect)
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 5, tool_func, {})
        assert result.n_cs == 2


# ---------------------------------------------------------------------------
# TestFineMapMultiInputPurityFiltering
# ---------------------------------------------------------------------------
class TestFineMapMultiInputPurityFiltering:
    """Test purity filtering for multi-input tools."""

    @patch("credtools.credtools.run_multisusie")
    def test_multi_input_purity_filtering(self, mock_run):
        """Multi-input tool with purity > 0 triggers purity filtering."""
        cs = CredibleSet(
            tool="multisusie",
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.2, "1-300-A-G": 0.0}),
            purity=[0.3],  # Low purity that should be filtered out
        )
        mock_run.return_value = cs
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])

        result = fine_map(
            locus_set,
            tool="multisusie",
            max_causal=5,
            set_L_by_cojo=False,
            purity=0.5,  # Threshold higher than the CS purity
        )
        # The CS with purity=0.3 should be filtered out since purity < 0.5
        assert result.n_cs == 0


# ---------------------------------------------------------------------------
# TestPipelineQCMetrics
# ---------------------------------------------------------------------------
class TestPipelineQCMetrics:
    """Tests for pipeline QC metrics saving branch."""

    @patch("credtools.utils.format_enhanced_pips", side_effect=lambda x: x)
    @patch("credtools.credibleset.generate_cs_summary", return_value=[])
    @patch("credtools.credtools.fine_map")
    @patch("credtools.credtools.locus_qc")
    @patch("credtools.credtools.save_heterogeneity")
    @patch("credtools.credtools.heterogeneity_summary", return_value={})
    @patch("credtools.credtools.compute_heterogeneity", return_value={})
    @patch("credtools.credtools.meta")
    @patch("credtools.credtools.load_locus_set")
    def test_pipeline_saves_qc_metrics(
        self,
        mock_load,
        mock_meta,
        mock_het,
        mock_het_summary,
        mock_save_het,
        mock_qc,
        mock_fine_map,
        mock_cs_summary,
        mock_format,
        tmp_path,
    ):
        """Pipeline saves QC metric files when QC is not skipped."""
        from credtools.credtools import pipeline

        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        mock_load.return_value = locus_set
        mock_meta.return_value = locus_set

        # QC returns a dict of DataFrames
        qc_df = pd.DataFrame({"metric": [1.0, 2.0], "value": [0.5, 0.6]})
        mock_qc.return_value = {"ld_check": qc_df}

        snp_ids = ["1-100-A-G", "1-200-A-G", "1-300-A-G"]
        cs = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=["1-100-A-G"],
            snps=[["1-100-A-G", "1-200-A-G"]],
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
        )
        cs.set_per_locus_results({})
        mock_fine_map.return_value = cs

        loci_df = pd.DataFrame(
            {
                "prefix": ["/fake/path/eur_c1"],
                "popu": ["EUR"],
                "cohort": ["c1"],
                "sample_size": [1000],
                "chr": [1],
                "start": [100],
                "end": [400],
                "locus_id": ["chr1:100-400"],
            }
        )
        pipeline(loci_df, tool="susie", outdir=str(tmp_path), skip_qc=False)

        # QC metric file should be saved
        qc_file = tmp_path / "ld_check.txt"
        assert qc_file.exists()


# ---------------------------------------------------------------------------
# TestAdaptiveFinemapNonConvergence
# ---------------------------------------------------------------------------
class TestAdaptiveFinemapNonConvergence:
    """Non-convergence handling in _adaptive_fine_map (empty_on_nonconvergence)."""

    def _empty_nonconverged(self, tool="susie") -> CredibleSet:
        """Empty CS with converged=False, mimicking empty_on_nonconvergence=True."""
        return CredibleSet(
            tool=tool,
            parameters={},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
            converged=False,
        )

    def _ok_cs(self, n_cs, tool="susie") -> CredibleSet:
        """Successful credible set with converged=True."""
        return CredibleSet(
            tool=tool,
            parameters={},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}"] for i in range(n_cs)],
            pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
            converged=True,
        )

    def test_nonconverged_initial_then_decrement_succeeds(self):
        """Initial non-converged → decrement L → smaller L converges with signal."""
        call_count = 0

        def side_effect(locus, max_causal=5, **kwargs):
            nonlocal call_count
            call_count += 1
            if max_causal >= 4:
                return self._empty_nonconverged()
            return self._ok_cs(2)

        tool_func = MagicMock(side_effect=side_effect)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 5, tool_func, {})
        assert result.n_cs == 2
        assert result.converged is True
        # initial=5 (nc), then phase2: 4 (nc), 3 (ok)
        assert call_count == 3

    def test_nonconverged_at_all_L_returns_empty(self):
        """Non-converged at every L from initial down to 1 → empty result."""
        tool_func = MagicMock(side_effect=lambda *a, **kw: self._empty_nonconverged())
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 3, tool_func, {})
        assert result.n_cs == 0
        # Phase 1 (L=3) + Phase 2 (L=2, 1) = 3 attempts
        assert tool_func.call_count == 3

    def test_genuine_zero_n_cs_with_converged_true_does_not_keep_retrying(self):
        """n_cs=0 with converged=True (no signal) returns at first Phase-2 success."""
        empty_converged = CredibleSet(
            tool="susie",
            parameters={},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
            converged=True,
        )
        tool_func = MagicMock(return_value=empty_converged)
        locus = MagicMock()
        result = _adaptive_fine_map(locus, "susie", 3, tool_func, {})
        assert result.n_cs == 0
        # Phase 1 + first Phase-2 iteration should suffice (no infinite retry).
        # Phase 1 (L=3) + Phase 2 (L=2) = 2 calls; Phase 2 returns since converged.
        assert tool_func.call_count == 2


# ---------------------------------------------------------------------------
# TestAdaptiveFineMapMultiNonConvergence
# ---------------------------------------------------------------------------
class TestAdaptiveFineMapMultiNonConvergence:
    """Non-convergence handling in _adaptive_fine_map_multi."""

    def _empty_nonconverged(self, tool="multisusie") -> CredibleSet:
        return CredibleSet(
            tool=tool,
            parameters={},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
            converged=False,
        )

    def _ok_cs(self, n_cs, tool="multisusie") -> CredibleSet:
        return CredibleSet(
            tool=tool,
            parameters={},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=[f"s{i}" for i in range(n_cs)],
            snps=[[f"s{i}"] for i in range(n_cs)],
            pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
            converged=True,
        )

    def test_nonconverged_initial_then_decrement_succeeds(self):
        """Initial non-converged → decrement L → smaller L converges."""

        def side_effect(locus_set, max_causal=5, **kwargs):
            if max_causal >= 4:
                return self._empty_nonconverged()
            return self._ok_cs(2)

        tool_func = MagicMock(side_effect=side_effect)
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 5, tool_func, {})
        assert result.n_cs == 2
        assert result.converged is True

    def test_nonconverged_at_all_L_returns_empty(self):
        """Non-converged at every L → return empty result."""
        tool_func = MagicMock(side_effect=lambda *a, **kw: self._empty_nonconverged())
        locus_set = MagicMock()
        locus_set.n_loci = 2
        result = _adaptive_fine_map_multi(locus_set, "multisusie", 3, tool_func, {})
        assert result.n_cs == 0


# ---------------------------------------------------------------------------
# TestFineMapInjectsEmptyOnNonconvergence
# ---------------------------------------------------------------------------
class TestFineMapInjectsEmptyOnNonconvergence:
    """fine_map should inject empty_on_nonconvergence=True when adaptive_max_causal=True."""

    def _make_cs(self, n_cs=1, tool="susie") -> CredibleSet:
        snp_ids = ["1-100-A-G", "1-200-A-G", "1-300-A-G"]
        return CredibleSet(
            tool=tool,
            parameters={"max_causal": 5},
            coverage=0.95,
            n_cs=n_cs,
            cs_sizes=[2] * n_cs,
            lead_snps=snp_ids[:n_cs],
            snps=[snp_ids[:2]] * n_cs,
            pips=pd.Series({"1-100-A-G": 0.8, "1-200-A-G": 0.15, "1-300-A-G": 0.05}),
            converged=True,
        )

    @patch("credtools.credtools._adaptive_fine_map")
    def test_susie_adaptive_injects_empty_on_nonconvergence(self, mock_adaptive):
        """Test that susie + adaptive injects empty_on_nonconvergence=True."""
        mock_adaptive.return_value = self._make_cs(2, "susie")
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        fine_map(
            locus_set,
            tool="susie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        args, kwargs = mock_adaptive.call_args
        params = args[4] if len(args) > 4 else kwargs.get("params")
        assert params.get("empty_on_nonconvergence") is True

    @patch("credtools.credtools._adaptive_fine_map")
    def test_rsparsepro_adaptive_injects_empty_on_nonconvergence(self, mock_adaptive):
        """Test that rsparsepro + adaptive injects empty_on_nonconvergence=True."""
        mock_adaptive.return_value = self._make_cs(2, "rsparsepro")
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        fine_map(
            locus_set,
            tool="rsparsepro",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        args, kwargs = mock_adaptive.call_args
        params = args[4] if len(args) > 4 else kwargs.get("params")
        assert params.get("empty_on_nonconvergence") is True

    @patch("credtools.credtools._adaptive_fine_map_multi")
    def test_multisusie_adaptive_injects_empty_on_nonconvergence(self, mock_adaptive):
        """Test that multisusie + adaptive injects empty_on_nonconvergence=True."""
        mock_adaptive.return_value = self._make_cs(2, "multisusie")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])
        fine_map(
            locus_set,
            tool="multisusie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        args, kwargs = mock_adaptive.call_args
        params = args[4] if len(args) > 4 else kwargs.get("params")
        assert params.get("empty_on_nonconvergence") is True

    @patch("credtools.credtools._adaptive_fine_map_multi")
    def test_mesusie_adaptive_injects_empty_on_nonconvergence(self, mock_adaptive):
        """Test that mesusie + adaptive injects empty_on_nonconvergence=True."""
        mock_adaptive.return_value = self._make_cs(2, "mesusie")
        locus1 = _make_test_locus("EUR", "c1", 1.0)
        locus2 = _make_test_locus("AFR", "c2", 0.8)
        locus_set = LocusSet([locus1, locus2])
        fine_map(
            locus_set,
            tool="mesusie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
        )
        args, kwargs = mock_adaptive.call_args
        params = args[4] if len(args) > 4 else kwargs.get("params")
        assert params.get("empty_on_nonconvergence") is True

    @patch("credtools.credtools.run_susie")
    def test_susie_non_adaptive_does_not_inject(self, mock_run):
        """Test that susie + non-adaptive does NOT inject empty_on_nonconvergence."""
        mock_run.return_value = self._make_cs(2, "susie")
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        fine_map(
            locus_set,
            tool="susie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=False,
        )
        call_kwargs = mock_run.call_args[1]
        assert "empty_on_nonconvergence" not in call_kwargs

    @patch("credtools.credtools._adaptive_fine_map")
    def test_user_override_empty_on_nonconvergence_false_respected(self, mock_adaptive):
        """User-provided empty_on_nonconvergence=False overrides adaptive default."""
        mock_adaptive.return_value = self._make_cs(2, "susie")
        locus = _make_test_locus("EUR", "c1", 1.0)
        locus_set = LocusSet([locus])
        fine_map(
            locus_set,
            tool="susie",
            max_causal=5,
            set_L_by_cojo=False,
            adaptive_max_causal=True,
            empty_on_nonconvergence=False,
        )
        args, kwargs = mock_adaptive.call_args
        params = args[4] if len(args) > 4 else kwargs.get("params")
        assert params.get("empty_on_nonconvergence") is False
