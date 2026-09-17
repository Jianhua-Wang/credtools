"""Regression tests for the adaptive max_causal (L) control flow.

These tests exercise ``_adaptive_fine_map`` and ``_adaptive_fine_map_multi``
with a scripted ``tool_func`` so that the exact sequence of attempted L values
can be asserted. Both adaptive entry points are covered through the same
parametrized scenarios to guarantee they share identical control logic.
"""

from typing import Callable, Dict, List, Optional
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from credtools.credibleset import CredibleSet
from credtools.credtools import (_adaptive_fine_map, _adaptive_fine_map_multi,
                                 _is_valid_result, fine_map)


def _make_cs(
    n_cs: int,
    max_causal: int,
    converged: Optional[bool] = True,
    tool: str = "susie",
    purities: Optional[List[float]] = None,
) -> CredibleSet:
    return CredibleSet(
        tool=tool,
        parameters={"max_causal": max_causal},
        coverage=0.95,
        n_cs=n_cs,
        cs_sizes=[2] * n_cs,
        lead_snps=[f"s{i}" for i in range(n_cs)],
        snps=[[f"s{i}", f"s{i}_b"] for i in range(n_cs)],
        pips=pd.Series({f"s{i}": 0.8 for i in range(n_cs)}),
        purity=purities,
        converged=converged,
    )


class ScriptedTool:
    """A tool_func stand-in driven by a ``{L: outcome}`` script.

    Outcome values:
      - int ``k``: return a converged result with ``n_cs=k``
      - ``"raise"``: raise RuntimeError
      - ``"nonconverged"``: return ``n_cs=0, converged=False``
      - ``"empty"``: return ``n_cs=0, converged=None`` (e.g. significance gate)
      - ``("cs", k, converged)``: return ``n_cs=k`` with explicit converged flag
    Any L not in the script raises, so unexpected attempts are loud.
    """

    def __init__(self, script: Dict[int, object], tool: str = "susie"):
        self.script = script
        self.tool = tool
        self.calls: List[int] = []

    def __call__(self, data, max_causal: int, **kwargs) -> CredibleSet:
        self.calls.append(max_causal)
        if max_causal not in self.script:
            raise AssertionError(f"unexpected attempt with max_causal={max_causal}")
        outcome = self.script[max_causal]
        if outcome == "raise":
            raise RuntimeError(f"tool failed at L={max_causal}")
        if outcome == "nonconverged":
            return _make_cs(0, max_causal, converged=False, tool=self.tool)
        if outcome == "empty":
            return _make_cs(0, max_causal, converged=None, tool=self.tool)
        if isinstance(outcome, tuple):
            _, k, conv = outcome
            return _make_cs(k, max_causal, converged=conv, tool=self.tool)
        return _make_cs(int(outcome), max_causal, converged=True, tool=self.tool)


def _run_single(tool: ScriptedTool, initial: int, **kw) -> CredibleSet:
    return _adaptive_fine_map(MagicMock(), "susie", initial, tool, {}, **kw)


def _run_multi(tool: ScriptedTool, initial: int, **kw) -> CredibleSet:
    locus_set = MagicMock()
    locus_set.n_loci = 2
    return _adaptive_fine_map_multi(locus_set, "multisusie", initial, tool, {}, **kw)


RUNNERS = [
    pytest.param(_run_single, id="single"),
    pytest.param(_run_multi, id="multi"),
]


@pytest.mark.parametrize("runner", RUNNERS)
class TestAdaptiveControlFlow:
    """Scenario tests shared by the single- and multi-input adaptive paths."""

    def test_fallback_starts_from_failed_L_minus_one(self, runner: Callable):
        """5 sat -> 10 sat -> 15 fails -> 14 succeeds: return the L=14 result."""
        tool = ScriptedTool({5: 5, 10: 10, 15: "raise", 14: 6})
        result = runner(tool, 5)
        assert tool.calls == [5, 10, 15, 14]
        assert result.n_cs == 6
        assert result.parameters["max_causal"] == 14

    def test_fallback_reuses_previous_successful_L(self, runner: Callable):
        """5 sat -> 10 sat -> 15 fails -> 14..11 fail: reuse the cached L=10 result."""
        tool = ScriptedTool(
            {
                5: 5,
                10: 10,
                15: "raise",
                14: "raise",
                13: "nonconverged",
                12: "raise",
                11: "nonconverged",
            }
        )
        result = runner(tool, 5)
        assert tool.calls == [5, 10, 15, 14, 13, 12, 11]
        assert result.n_cs == 10
        assert result.parameters["max_causal"] == 10
        assert result.parameters.get("adaptive_failed") is None

    def test_initial_failure_then_decrease(self, runner: Callable):
        """Initial L=5 fails -> L=4 succeeds (allmeta-style behavior preserved)."""
        tool = ScriptedTool({5: "raise", 4: 4})
        result = runner(tool, 5)
        assert tool.calls == [5, 4]
        assert result.n_cs == 4
        assert result.parameters["max_causal"] == 4

    def test_initial_nonconverged_then_decrease(self, runner: Callable):
        """Initial L=5 non-converged empty -> L=4 non-converged -> L=3 succeeds."""
        tool = ScriptedTool({5: "nonconverged", 4: "nonconverged", 3: 2})
        result = runner(tool, 5)
        assert tool.calls == [5, 4, 3]
        assert result.n_cs == 2

    def test_saturated_up_to_cap_never_exceeds_20(self, runner: Callable):
        """Saturated at 5/10/15/20: return the L=20 result, never try 25 or fall back."""
        tool = ScriptedTool({5: 5, 10: 10, 15: 15, 20: 20})
        result = runner(tool, 5)
        assert tool.calls == [5, 10, 15, 20]
        assert result.n_cs == 20
        assert result.parameters["max_causal"] == 20

    def test_cap_is_clamped_for_non_multiple_initial(self, runner: Callable):
        """Initial L=3 steps 8 -> 13 -> 18 -> 20 (clamped), never 23."""
        tool = ScriptedTool({3: 3, 8: 8, 13: 13, 18: 18, 20: 12})
        result = runner(tool, 3)
        assert tool.calls == [3, 8, 13, 18, 20]
        assert result.n_cs == 12

    def test_expansion_nonconverged_empty_is_failure(self, runner: Callable):
        """5 sat -> 10 returns non-converged empty: treat as failure, fall back from 9."""
        tool = ScriptedTool({5: 5, 10: "nonconverged", 9: 7})
        result = runner(tool, 5)
        assert tool.calls == [5, 10, 9]
        assert result.n_cs == 7
        assert result.parameters["max_causal"] == 9

    def test_expansion_nonconverged_then_reuse_cached(self, runner: Callable):
        """5 sat -> 10 non-converged -> 9..6 fail: reuse L=5 rather than jump to 4."""
        tool = ScriptedTool(
            {
                5: 5,
                10: "nonconverged",
                9: "raise",
                8: "raise",
                7: "nonconverged",
                6: "raise",
            }
        )
        result = runner(tool, 5)
        assert tool.calls == [5, 10, 9, 8, 7, 6]
        assert result.n_cs == 5
        assert result.parameters["max_causal"] == 5

    def test_nonconverged_with_credible_sets_is_not_success(self, runner: Callable):
        """converged=False with n_cs>0 (empty_on_nonconvergence=False) is not accepted."""
        tool = ScriptedTool({5: ("cs", 3, False), 4: 2})
        result = runner(tool, 5)
        assert tool.calls == [5, 4]
        assert result.n_cs == 2

    def test_genuine_empty_result_is_accepted(self, runner: Callable):
        """n_cs=0 with converged=None (significance gate) is a valid final result."""
        tool = ScriptedTool({5: "empty"})
        result = runner(tool, 5)
        assert tool.calls == [5]
        assert result.n_cs == 0
        assert result.converged is None
        assert result.parameters.get("adaptive_failed") is None

    def test_converged_empty_result_is_accepted(self, runner: Callable):
        """n_cs=0 with converged=True is a genuine no-signal result, not a failure."""
        tool = ScriptedTool({5: ("cs", 0, True)})
        result = runner(tool, 5)
        assert tool.calls == [5]
        assert result.n_cs == 0
        assert result.parameters.get("adaptive_failed") is None

    def test_all_attempts_fail(self, runner: Callable):
        """Every L from initial down to 1 fails -> empty result flagged adaptive_failed."""
        tool = ScriptedTool({3: "raise", 2: "nonconverged", 1: "raise"})
        result = runner(tool, 3)
        assert tool.calls == [3, 2, 1]
        assert result.n_cs == 0
        assert result.parameters.get("adaptive_failed") is True

    def test_success_on_first_try(self, runner: Callable):
        tool = ScriptedTool({5: 3})
        result = runner(tool, 5)
        assert tool.calls == [5]
        assert result.n_cs == 3

    def test_purity_filter_applied_before_decisions(self, runner: Callable):
        """Purity filtering still gates saturation decisions inside the loop."""

        def tool(data, max_causal, **kwargs):
            # 5 CS at L=5 but only 2 pass purity -> not saturated -> accept.
            return _make_cs(5, max_causal, purities=[0.9, 0.1, 0.8, 0.2, 0.1])

        mock = MagicMock(side_effect=tool)
        result = runner(mock, 5, purity_threshold=0.5)
        assert result.n_cs == 2
        assert mock.call_count == 1

    def test_log_records_attempts_and_final_choice(self, runner: Callable, caplog):
        tool = ScriptedTool({5: 5, 10: "raise", 9: 4})
        with caplog.at_level("INFO", logger="CREDTOOLS"):
            runner(tool, 5)
        text = caplog.text
        assert "max_causal=5" in text and "n_cs=5" in text
        assert "max_causal=10" in text and "failed" in text.lower()
        assert "max_causal=9" in text
        assert "selected max_causal=9" in text


class TestIsValidResult:
    """Tests for _is_valid_result."""

    def test_converged_true(self):
        assert _is_valid_result(_make_cs(2, 5, converged=True))

    def test_converged_none(self):
        assert _is_valid_result(_make_cs(0, 5, converged=None))

    def test_converged_false(self):
        assert not _is_valid_result(_make_cs(0, 5, converged=False))
        assert not _is_valid_result(_make_cs(3, 5, converged=False))


class TestNonAdaptiveUnaffected:
    """Non-adaptive fine_map must call the tool exactly once with the given L."""

    @patch("credtools.credtools.run_susie")
    def test_single_call_with_requested_L(self, mock_run):
        mock_run.return_value = _make_cs(0, 7, converged=False)
        locus = MagicMock()
        locus_set = MagicMock()
        locus_set.n_loci = 1
        locus_set.loci = [locus]
        result = fine_map(
            locus_set,
            tool="susie",
            max_causal=7,
            adaptive_max_causal=False,
            set_L_by_cojo=False,
        )
        assert mock_run.call_count == 1
        assert mock_run.call_args.kwargs["max_causal"] == 7
        # Non-adaptive returns whatever the tool produced, even if non-converged.
        assert result.converged is False
        # Non-adaptive does not inject the empty_on_nonconvergence default.
        assert "empty_on_nonconvergence" not in mock_run.call_args.kwargs
