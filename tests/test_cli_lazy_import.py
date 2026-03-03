"""Tests for lazy import optimization in cli.py.

Verify that importing credtools.cli does NOT trigger loading heavy modules
like numpy, pandas, scipy, or sklearn at module level.
"""

import subprocess
import sys


def test_cli_module_does_not_import_pandas_at_top_level():
    """Importing credtools.cli should not trigger pandas import."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; "
            "sys.modules.pop('pandas', None); "
            "import credtools.cli; "
            "assert 'pandas' not in sys.modules, "
            "'pandas was imported at module level by credtools.cli'",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"pandas imported at module level:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )


def test_cli_module_does_not_import_numpy_at_top_level():
    """Importing credtools.cli should not trigger numpy import."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; "
            "sys.modules.pop('numpy', None); "
            "import credtools.cli; "
            "assert 'numpy' not in sys.modules, "
            "'numpy was imported at module level by credtools.cli'",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"numpy imported at module level:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )


def test_cli_module_does_not_import_scipy_at_top_level():
    """Importing credtools.cli should not trigger scipy import."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; "
            "import credtools.cli; "
            "assert 'scipy' not in sys.modules, "
            "'scipy was imported at module level by credtools.cli'",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"scipy imported at module level:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )


def test_cli_module_does_not_import_sklearn_at_top_level():
    """Importing credtools.cli should not trigger sklearn import."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; "
            "import credtools.cli; "
            "assert 'sklearn' not in sys.modules, "
            "'sklearn was imported at module level by credtools.cli'",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"sklearn imported at module level:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
