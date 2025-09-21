"""Tests for timeout handling in credtools."""

import sys
from types import SimpleNamespace

import pytest

from credtools.credtools import fine_map
from credtools.utils import ExternalTool


def test_external_tool_run_timeout(tmp_path):
    """ExternalTool.run should raise TimeoutError and log when timeout is reached."""
    log_file = tmp_path / "timeout.log"
    tool = ExternalTool("python")
    tool.set_custom_path(sys.executable)

    with pytest.raises(TimeoutError) as exc:
        tool.run(
            ["-c", "import time; time.sleep(0.5)"],
            str(log_file),
            timeout=0.1,
        )

    assert "timed out" in str(exc.value)
    assert "[timeout]" in log_file.read_text()


def test_fine_map_default_timeout_for_finemap(monkeypatch):
    """fine_map should default to 30 minutes per locus for FINEMAP when timeout is not provided."""
    captured = {}

    def fake_run_finemap(locus, max_causal, timeout_minutes=None, **kwargs):
        captured["timeout_minutes"] = timeout_minutes
        captured["extra_kwargs"] = kwargs
        return "mock-result"

    monkeypatch.setattr("credtools.credtools.run_finemap", fake_run_finemap)

    dummy_locus_set = SimpleNamespace(n_loci=1, loci=[object()])

    result = fine_map(
        dummy_locus_set,
        tool="finemap",
        set_L_by_cojo=False,
        timeout_minutes=None,
    )

    assert result == "mock-result"
    assert captured["timeout_minutes"] == 30.0
    assert captured.get("extra_kwargs", {}) == {}

    captured.clear()
    result = fine_map(
        dummy_locus_set,
        tool="finemap",
        set_L_by_cojo=False,
        timeout_minutes=15,
    )

    assert result == "mock-result"
    assert captured["timeout_minutes"] == 15.0
    assert captured.get("extra_kwargs", {}) == {}
