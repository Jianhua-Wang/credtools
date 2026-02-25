"""Tests for credtools.utils module."""

import logging
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pytest

from credtools.utils import (
    ExternalTool,
    ToolManager,
    create_float_format_dict,
    format_enhanced_pips,
    format_float,
    format_pvalue,
    get_float_format,
    io_in_tempdir,
    tool_manager,
)


# ---------------------------------------------------------------------------
# TestFormatFloat
# ---------------------------------------------------------------------------
class TestFormatFloat:
    """Tests for format_float function."""

    def test_nan_returns_empty(self):
        assert format_float(float("nan")) == ""

    def test_none_returns_empty(self):
        assert format_float(None) == ""

    def test_pd_na_returns_empty(self):
        assert format_float(pd.NA) == ""

    def test_non_numeric_returns_str(self):
        assert format_float("hello") == "hello"

    def test_small_number_scientific(self):
        result = format_float(1e-6)
        assert "e" in result.lower()

    def test_large_number_scientific(self):
        result = format_float(99999)
        assert "e" in result.lower()

    def test_regular_number_decimal(self):
        result = format_float(0.1234)
        assert result == "0.1234"

    def test_custom_decimals(self):
        result = format_float(0.123456, decimals=2)
        assert result == "0.12"

    def test_zero(self):
        result = format_float(0.0)
        assert "e" in result.lower() or result == "0.0000"

    def test_negative_number(self):
        result = format_float(-0.5)
        assert result == "-0.5000"

    def test_string_number(self):
        result = format_float("0.5")
        assert result == "0.5000"


# ---------------------------------------------------------------------------
# TestFormatPvalue
# ---------------------------------------------------------------------------
class TestFormatPvalue:
    """Tests for format_pvalue function."""

    def test_nan_returns_empty(self):
        assert format_pvalue(float("nan")) == ""

    def test_normal_pvalue(self):
        result = format_pvalue(5e-8)
        assert "e" in result.lower()
        assert result == "5.000e-08"

    def test_large_pvalue(self):
        result = format_pvalue(0.5)
        assert result == "5.000e-01"


# ---------------------------------------------------------------------------
# TestGetFloatFormat
# ---------------------------------------------------------------------------
class TestGetFloatFormat:
    """Tests for get_float_format function."""

    @pytest.mark.parametrize(
        "col_name, expected",
        [
            ("EUR_p", "%.3e"),
            ("p", "%.3e"),
            ("P", "%.3e"),
            ("cohort_P", "%.3e"),
        ],
    )
    def test_pvalue_columns(self, col_name, expected):
        assert get_float_format(col_name) == expected

    @pytest.mark.parametrize(
        "col_name, expected",
        [
            ("EUR_eaf", "%.4f"),
            ("eaf", "%.4f"),
            ("cohort_maf", "%.4f"),
            ("maf", "%.4f"),
            ("EUR_pip", "%.4f"),
            ("pip", "%.4f"),
            ("cohort_r2", "%.4f"),
            ("r2", "%.4f"),
        ],
    )
    def test_frequency_columns(self, col_name, expected):
        assert get_float_format(col_name) == expected

    @pytest.mark.parametrize(
        "col_name, expected",
        [
            ("EUR_beta", "%.4f"),
            ("beta", "%.4f"),
            ("cohort_se", "%.4f"),
            ("se", "%.4f"),
        ],
    )
    def test_effect_columns(self, col_name, expected):
        assert get_float_format(col_name) == expected

    def test_unknown_column_returns_none(self):
        assert get_float_format("snpid") is None
        assert get_float_format("chr") is None


# ---------------------------------------------------------------------------
# TestCreateFloatFormatDict
# ---------------------------------------------------------------------------
class TestCreateFloatFormatDict:
    """Tests for create_float_format_dict function."""

    def test_mixed_type_dataframe(self):
        df = pd.DataFrame(
            {
                "snpid": ["rs1", "rs2"],
                "EUR_p": [1e-8, 0.05],
                "EUR_beta": [0.1, 0.2],
                "chr": [1, 2],
            }
        )
        result = create_float_format_dict(df)
        assert "EUR_p" in result
        assert result["EUR_p"] == "%.3e"
        assert "EUR_beta" in result
        assert result["EUR_beta"] == "%.4f"
        assert "snpid" not in result

    def test_empty_dataframe(self):
        df = pd.DataFrame()
        result = create_float_format_dict(df)
        assert result == {}

    def test_no_numeric_columns(self):
        df = pd.DataFrame({"a": ["x", "y"], "b": ["z", "w"]})
        result = create_float_format_dict(df)
        assert result == {}


# ---------------------------------------------------------------------------
# TestFormatEnhancedPips
# ---------------------------------------------------------------------------
class TestFormatEnhancedPips:
    """Tests for format_enhanced_pips function."""

    def test_pvalue_scientific_notation(self):
        df = pd.DataFrame({"EUR_p": [1e-8, 0.05], "snpid": ["rs1", "rs2"]})
        result = format_enhanced_pips(df)
        assert "e" in str(result["EUR_p"].iloc[0]).lower()

    def test_pip_4_decimals(self):
        df = pd.DataFrame({"pip": [0.123456, 0.789], "snpid": ["rs1", "rs2"]})
        result = format_enhanced_pips(df)
        assert result["pip"].iloc[0] == "0.1235"

    def test_nan_values_become_empty(self):
        df = pd.DataFrame({"EUR_p": [1e-8, float("nan")], "snpid": ["rs1", "rs2"]})
        result = format_enhanced_pips(df)
        assert result["EUR_p"].iloc[1] == ""

    def test_non_matching_columns_unchanged(self):
        df = pd.DataFrame({"chr": [1, 2], "snpid": ["rs1", "rs2"]})
        result = format_enhanced_pips(df)
        assert result["chr"].iloc[0] == 1


# ---------------------------------------------------------------------------
# TestExternalTool
# ---------------------------------------------------------------------------
class TestExternalTool:
    """Tests for ExternalTool class."""

    def test_set_custom_path_nonexistent_raises(self):
        tool = ExternalTool("fake_tool")
        with pytest.raises(FileNotFoundError, match="does not exist"):
            tool.set_custom_path("/nonexistent/path/fake_tool")

    def test_set_custom_path_existing(self, tmp_path):
        tool = ExternalTool("test_tool")
        fake_bin = tmp_path / "test_tool"
        fake_bin.touch()
        tool.set_custom_path(str(fake_bin))
        assert tool.custom_path == str(fake_bin)

    def test_get_path_custom(self, tmp_path):
        tool = ExternalTool("test_tool")
        fake_bin = tmp_path / "test_tool"
        fake_bin.touch()
        tool.set_custom_path(str(fake_bin))
        assert tool.get_path() == str(fake_bin)

    def test_get_path_system(self):
        # python is always in system PATH
        tool = ExternalTool("python3")
        path = tool.get_path()
        assert "python" in path.lower()

    def test_get_path_not_found_raises(self):
        tool = ExternalTool("totally_nonexistent_tool_xyz_12345")
        with pytest.raises(FileNotFoundError, match="Could not find"):
            tool.get_path()

    def test_run_negative_timeout_raises(self, tmp_path):
        tool = ExternalTool("python3")
        log = str(tmp_path / "test.log")
        with pytest.raises(ValueError, match="timeout must be a positive"):
            tool.run(["--version"], log, timeout=-1)

    def test_run_success(self, tmp_path):
        tool = ExternalTool("python3")
        log = str(tmp_path / "test.log")
        tool.run(["--version"], log)
        assert os.path.exists(log)

    def test_run_output_file_check(self, tmp_path):
        tool = ExternalTool("python3")
        log = str(tmp_path / "test.log")
        with pytest.raises(FileNotFoundError, match="Expected output file"):
            tool.run(["--version"], log, output_file_path="/nonexistent/output.txt")

    def test_run_output_file_list(self, tmp_path):
        tool = ExternalTool("python3")
        log = str(tmp_path / "test.log")
        with pytest.raises(FileNotFoundError, match="Expected output file"):
            tool.run(
                ["--version"],
                log,
                output_file_path=["/nonexistent/a.txt", "/nonexistent/b.txt"],
            )

    def test_run_timeout(self, tmp_path):
        tool = ExternalTool("python3")
        log = str(tmp_path / "test.log")
        with pytest.raises(TimeoutError, match="timed out"):
            tool.run(["-c", "import time; time.sleep(30)"], log, timeout=0.1)

    def test_run_command_failure(self, tmp_path):
        tool = ExternalTool("python3")
        log = str(tmp_path / "test.log")
        with pytest.raises(subprocess.CalledProcessError):
            tool.run(["-c", "raise SystemExit(1)"], log)


# ---------------------------------------------------------------------------
# TestToolManager
# ---------------------------------------------------------------------------
class TestToolManager:
    """Tests for ToolManager class."""

    def test_set_unregistered_tool_raises(self):
        mgr = ToolManager()
        with pytest.raises(KeyError, match="not registered"):
            mgr.set_tool_path("nonexistent", "/fake")

    def test_get_unregistered_tool_raises(self):
        mgr = ToolManager()
        with pytest.raises(KeyError, match="not registered"):
            mgr.get_tool("nonexistent")

    def test_run_unregistered_tool_raises(self):
        mgr = ToolManager()
        with pytest.raises(KeyError, match="not registered"):
            mgr.run_tool("nonexistent", [], "log.txt")

    def test_register_and_get(self):
        mgr = ToolManager()
        mgr.register_tool("my_tool")
        tool = mgr.get_tool("my_tool")
        assert isinstance(tool, ExternalTool)
        assert tool.name == "my_tool"

    def test_register_with_default_path(self):
        mgr = ToolManager()
        mgr.register_tool("my_tool", "bin/my_tool")
        tool = mgr.get_tool("my_tool")
        assert tool.default_path == "bin/my_tool"

    def test_global_tool_manager_has_finemap(self):
        assert "finemap" in tool_manager.tools
        assert "SuSiEx" in tool_manager.tools


# ---------------------------------------------------------------------------
# TestIoInTempdir
# ---------------------------------------------------------------------------
class TestIoInTempdir:
    """Tests for io_in_tempdir decorator."""

    def test_tempdir_created_and_cleaned(self, tmp_path):
        parent = str(tmp_path / "tmp")

        @io_in_tempdir(dir=parent)
        def dummy(temp_dir=None):
            assert os.path.isdir(temp_dir)
            return temp_dir

        # Set logger level to INFO so temp dir is cleaned
        logging.getLogger("Utils").setLevel(logging.INFO)
        used_dir = dummy()
        assert not os.path.exists(used_dir)

    def test_tempdir_retained_on_debug(self, tmp_path):
        parent = str(tmp_path / "tmp")

        @io_in_tempdir(dir=parent)
        def dummy(temp_dir=None):
            return temp_dir

        logging.getLogger("Utils").setLevel(logging.DEBUG)
        used_dir = dummy()
        assert os.path.exists(used_dir)
        # Cleanup
        logging.getLogger("Utils").setLevel(logging.WARNING)

    def test_exception_propagated(self, tmp_path):
        parent = str(tmp_path / "tmp")

        @io_in_tempdir(dir=parent)
        def fail(temp_dir=None):
            raise RuntimeError("test error")

        with pytest.raises(RuntimeError, match="test error"):
            fail()
