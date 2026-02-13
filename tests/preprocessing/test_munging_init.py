"""Unit tests for credtools/preprocessing/munging/__init__.py."""

import gzip
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from credtools.preprocessing.munging import (
    create_example_config,
    load_and_munge,
    read_config,
    validate_config,
)
from credtools.preprocessing.munging.constants import OUTPUT_COLS, ColName


# ==================== Tests for read_config() ====================


def test_read_config_valid_json(sample_config_file):
    """Test read_config() with valid JSON file returns dict."""
    config = read_config(sample_config_file)

    assert isinstance(config, dict)
    assert "column_mapping" in config
    assert isinstance(config["column_mapping"], dict)
    assert config["column_mapping"]["CHROM"] == "CHR"


def test_read_config_file_not_found():
    """Test read_config() raises FileNotFoundError when file doesn't exist."""
    non_existent_path = "/tmp/non_existent_config_12345.json"

    with pytest.raises(FileNotFoundError, match="Configuration file not found"):
        read_config(non_existent_path)


def test_read_config_invalid_json(tmp_path):
    """Test read_config() raises ValueError for invalid JSON."""
    invalid_json_file = tmp_path / "invalid.json"
    invalid_json_file.write_text("{ invalid json content }")

    with pytest.raises(ValueError, match="Invalid JSON in configuration file"):
        read_config(str(invalid_json_file))


# ==================== Tests for load_and_munge() ====================


def test_load_and_munge_tab_file_auto_mapping(tab_sumstats_file):
    """Test load_and_munge() with tab file, auto separator and auto mapping."""
    result = load_and_munge(tab_sumstats_file)

    # Should return a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Should have rows (may be filtered but should have some)
    assert len(result) > 0

    # Should have all output columns
    for col in OUTPUT_COLS:
        assert col in result.columns, f"Missing expected output column: {col}"

    # CHR and BP should be integer types (munge converts CHR to int8, BP to int32)
    assert result["CHR"].dtype in [np.int8, np.int16, np.int32, np.int64]
    assert result["BP"].dtype in [np.int8, np.int16, np.int32, np.int64]


def test_load_and_munge_with_config_mapping(tmp_path, sample_config_dict):
    """Test load_and_munge() applies config column_mapping correctly."""
    # Create a file with non-standard column names
    df = pd.DataFrame(
        {
            "CHROM": [1, 2],
            "POS": [1000, 2000],
            "A1": ["A", "C"],
            "A2": ["G", "T"],
            "BETA": [0.1, -0.2],
            "SE": [0.05, 0.06],
            "PVAL": [0.01, 0.001],
        }
    )
    filepath = tmp_path / "custom_headers.tsv"
    df.to_csv(filepath, sep="\t", index=False)

    # Load with config mapping
    result = load_and_munge(str(filepath), config=sample_config_dict)

    # Should have standard column names after mapping
    assert "CHR" in result.columns
    assert "BP" in result.columns
    assert "P" in result.columns
    # Original column names should not be in result
    assert "CHROM" not in result.columns
    assert "POS" not in result.columns
    assert "PVAL" not in result.columns


def test_load_and_munge_gzipped_file(gzipped_sumstats_file):
    """Test load_and_munge() successfully reads gzipped files."""
    result = load_and_munge(gzipped_sumstats_file)

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0

    # Should have all mandatory columns
    for col in ColName.mandatory_cols:
        assert col in result.columns


def test_load_and_munge_explicit_sep_parameter(tmp_path, minimal_gwas_df):
    """Test load_and_munge() uses explicitly provided sep parameter."""
    # Create comma-separated file
    filepath = tmp_path / "comma_sep.csv"
    minimal_gwas_df.to_csv(filepath, sep=",", index=False)

    # Load with explicit sep
    result = load_and_munge(str(filepath), sep=",")

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0


def test_load_and_munge_missing_mandatory_cols_raises_error(tmp_path):
    """Test load_and_munge() raises ValueError when mandatory columns missing."""
    # Create file with unmappable columns
    df = pd.DataFrame(
        {
            "foo": [1, 2, 3],
            "bar": [100, 200, 300],
            "baz": ["A", "B", "C"],
        }
    )
    filepath = tmp_path / "bad_headers.tsv"
    df.to_csv(filepath, sep="\t", index=False)

    with pytest.raises(ValueError, match="missing required columns"):
        load_and_munge(str(filepath))


def test_load_and_munge_output_has_expected_format(tab_sumstats_file):
    """Test load_and_munge() output has all 11 output_cols columns."""
    result = load_and_munge(tab_sumstats_file)

    # Should have exactly the output columns (11 columns)
    expected_cols = set(OUTPUT_COLS)
    actual_cols = set(result.columns)

    # Check all output columns are present
    assert expected_cols.issubset(
        actual_cols
    ), f"Missing columns: {expected_cols - actual_cols}"

    # Verify it's 11 columns as specified
    assert len(OUTPUT_COLS) == 11


# ==================== Tests for validate_config() ====================


def test_validate_config_valid_complete_config(sample_config_dict):
    """Test validate_config() returns True for complete valid config."""
    result = validate_config(sample_config_dict)

    assert result is True


def test_validate_config_missing_column_mapping_key():
    """Test validate_config() returns False when column_mapping key missing."""
    invalid_config = {"some_other_key": "value"}

    result = validate_config(invalid_config)

    assert result is False


def test_validate_config_column_mapping_not_dict():
    """Test validate_config() returns False when column_mapping is not a dict."""
    invalid_config = {"column_mapping": "not a dictionary"}

    result = validate_config(invalid_config)

    assert result is False


def test_validate_config_missing_mandatory_columns():
    """Test validate_config() returns False when mapping missing mandatory columns."""
    # Only map some columns, missing others
    incomplete_config = {
        "column_mapping": {
            "CHROM": "CHR",
            "POS": "BP",
            # Missing EA, NEA, BETA, SE, P
        }
    }

    result = validate_config(incomplete_config)

    assert result is False


# ==================== Tests for create_example_config() ====================


def test_create_example_config_creates_file(tmp_path):
    """Test create_example_config() creates JSON file at specified path."""
    headers = ["CHROM", "POS", "A1", "A2", "BETA", "SE", "PVAL"]
    output_path = tmp_path / "test_config.json"

    result_path = create_example_config(headers, output_path=str(output_path))

    # Check file was created
    assert os.path.exists(output_path)

    # Check it's valid JSON
    with open(output_path, "r") as f:
        config = json.load(f)

    assert isinstance(config, dict)
    assert "column_mapping" in config


def test_create_example_config_returns_file_path(tmp_path):
    """Test create_example_config() returns the file path."""
    headers = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
    output_path = tmp_path / "example_config.json"

    result_path = create_example_config(headers, output_path=str(output_path))

    # Should return the path
    assert result_path == str(output_path)

    # Path should be valid
    assert os.path.exists(result_path)
