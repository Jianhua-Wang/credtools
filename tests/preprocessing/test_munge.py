"""Unit tests for credtools.preprocessing.munge module."""

import gzip
import json
import os
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from credtools.preprocessing.munge import (
    create_munge_config,
    munge_sumstats,
    validate_munged_files,
)


@pytest.fixture
def valid_sumstats_file(tmp_path):
    """Create a valid temporary sumstats file."""
    df = pd.DataFrame(
        {
            "CHR": [1, 1, 2],
            "BP": [1000, 2000, 3000],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "C"],
            "BETA": [0.1, -0.2, 0.3],
            "SE": [0.05, 0.06, 0.07],
            "P": [0.01, 0.001, 0.0001],
        }
    )
    filepath = tmp_path / "test_sumstats.tsv"
    df.to_csv(filepath, sep="\t", index=False)
    return str(filepath)


@pytest.fixture
def munged_sumstats_file(tmp_path):
    """Create a munged sumstats file with all required columns."""
    df = pd.DataFrame(
        {
            "CHR": [1, 1, 2],
            "BP": [1000, 2000, 3000],
            "SNPID": ["1-1000-A-G", "1-2000-C-T", "2-3000-G-C"],
            "EA": ["A", "C", "G"],
            "NEA": ["G", "T", "C"],
            "EAF": [0.3, 0.45, 0.1],
            "BETA": [0.1, -0.2, 0.3],
            "SE": [0.05, 0.06, 0.07],
            "P": [0.01, 0.001, 0.0001],
            "RSID": ["rs1", "rs2", "rs3"],
        }
    )
    filepath = tmp_path / "test_munged.txt.gz"
    df.to_csv(filepath, sep="\t", index=False, compression="gzip")
    return str(filepath)


# The module path where load_and_munge/read_config/etc live (imported locally
# inside munge_sumstats via ``from .munging import ...``).
_MUNGING = "credtools.preprocessing.munging"


class TestMungeSumstatsInputNormalization:
    """Tests for munge_sumstats input normalization to dict format."""

    @patch(f"{_MUNGING}.load_and_munge")
    def test_string_input_normalized_to_dict(self, mock_load, tmp_path, valid_sumstats_file):
        """String input should be normalized to dict with filename as key."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        output_dir = str(tmp_path / "output")

        result = munge_sumstats(valid_sumstats_file, output_dir)

        expected_key = "test_sumstats"
        assert expected_key in result
        assert result[expected_key].endswith(f"{expected_key}.munged.txt.gz")

    @patch(f"{_MUNGING}.load_and_munge")
    def test_list_input_normalized_to_dict(self, mock_load, tmp_path, valid_sumstats_file):
        """List input should be normalized to dict with filenames as keys."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        file2 = tmp_path / "file2.tsv"
        pd.DataFrame({"CHR": [2]}).to_csv(file2, sep="\t", index=False)

        output_dir = str(tmp_path / "output")
        result = munge_sumstats([valid_sumstats_file, str(file2)], output_dir)

        assert "test_sumstats" in result
        assert "file2" in result
        assert len(result) == 2

    @patch(f"{_MUNGING}.load_and_munge")
    def test_dict_input_used_directly(self, mock_load, tmp_path, valid_sumstats_file):
        """Dict input should be used directly without modification."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        output_dir = str(tmp_path / "output")
        input_dict = {"EUR": valid_sumstats_file}

        result = munge_sumstats(input_dict, output_dir)

        assert "EUR" in result
        assert result["EUR"].endswith("EUR.munged.txt.gz")

    def test_invalid_type_raises_value_error(self, tmp_path):
        """Invalid input type (e.g., int) should raise ValueError."""
        output_dir = str(tmp_path / "output")

        with pytest.raises(ValueError, match="must be a string, list of strings, or dictionary"):
            munge_sumstats(12345, output_dir)


class TestMungeSumstatsFileHandling:
    """Tests for munge_sumstats file existence and error handling."""

    def test_file_not_found_raises_error(self, tmp_path):
        """Non-existent file should raise FileNotFoundError."""
        output_dir = str(tmp_path / "output")
        non_existent_file = str(tmp_path / "nonexistent.tsv")

        with pytest.raises(FileNotFoundError, match="Input file not found"):
            munge_sumstats(non_existent_file, output_dir)

    @patch(f"{_MUNGING}.load_and_munge")
    def test_output_dir_auto_created(self, mock_load, tmp_path, valid_sumstats_file):
        """Output directory should be created if it doesn't exist."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        output_dir = str(tmp_path / "nested" / "output" / "dir")

        munge_sumstats(valid_sumstats_file, output_dir)

        assert os.path.exists(output_dir)

    @patch(f"{_MUNGING}.load_and_munge")
    def test_output_filename_format(self, mock_load, tmp_path, valid_sumstats_file):
        """Output filename should be {identifier}.munged.txt.gz."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        output_dir = str(tmp_path / "output")

        result = munge_sumstats({"MyData": valid_sumstats_file}, output_dir)

        assert "MyData" in result
        assert result["MyData"].endswith("MyData.munged.txt.gz")
        assert os.path.dirname(result["MyData"]) == output_dir


class TestMungeSumstatsOutputFormat:
    """Tests for munge_sumstats output format and content."""

    @patch(f"{_MUNGING}.load_and_munge")
    def test_output_is_gzipped_tsv(self, mock_load, tmp_path, valid_sumstats_file):
        """Output file should be gzip-compressed TSV."""
        munged_data = pd.DataFrame(
            {
                "CHR": [1, 2],
                "BP": [1000, 2000],
                "EA": ["A", "C"],
                "NEA": ["G", "T"],
                "BETA": [0.1, -0.2],
                "SE": [0.05, 0.06],
                "P": [0.01, 0.001],
            }
        )
        mock_load.return_value = munged_data
        output_dir = str(tmp_path / "output")

        result = munge_sumstats(valid_sumstats_file, output_dir)

        output_file = result["test_sumstats"]
        assert output_file.endswith(".gz")

        # Verify it's readable as gzipped TSV
        df_read = pd.read_csv(output_file, sep="\t", compression="gzip")
        pd.testing.assert_frame_equal(df_read, munged_data)


class TestMungeSumstatsOverwriteBehavior:
    """Tests for munge_sumstats force_overwrite parameter."""

    @patch(f"{_MUNGING}.load_and_munge")
    def test_force_false_skips_existing_files(self, mock_load, tmp_path, valid_sumstats_file):
        """force_overwrite=False should skip existing output files without modification."""
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)

        # Create existing output file with specific content
        existing_output = os.path.join(output_dir, "test_sumstats.munged.txt.gz")
        original_data = pd.DataFrame({"ORIGINAL": [1, 2, 3]})
        original_data.to_csv(existing_output, sep="\t", index=False, compression="gzip")

        # Mock should return different data
        mock_load.return_value = pd.DataFrame({"NEW": [4, 5, 6]})

        # Run with force_overwrite=False
        result = munge_sumstats(valid_sumstats_file, output_dir, force_overwrite=False)

        # Verify file was not modified (still has original content)
        df_check = pd.read_csv(existing_output, sep="\t", compression="gzip")
        pd.testing.assert_frame_equal(df_check, original_data)

        # Mock should not have been called
        mock_load.assert_not_called()

    @patch(f"{_MUNGING}.load_and_munge")
    def test_force_true_overwrites_existing_files(self, mock_load, tmp_path, valid_sumstats_file):
        """force_overwrite=True should overwrite existing output files."""
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)

        # Create existing output file with specific content
        existing_output = os.path.join(output_dir, "test_sumstats.munged.txt.gz")
        original_data = pd.DataFrame({"ORIGINAL": [1, 2, 3]})
        original_data.to_csv(existing_output, sep="\t", index=False, compression="gzip")

        # Mock should return different data
        new_data = pd.DataFrame({"NEW": [4, 5, 6]})
        mock_load.return_value = new_data

        # Run with force_overwrite=True
        result = munge_sumstats(valid_sumstats_file, output_dir, force_overwrite=True)

        # Verify file was modified (now has new content)
        df_check = pd.read_csv(existing_output, sep="\t", compression="gzip")
        pd.testing.assert_frame_equal(df_check, new_data)

        # Mock should have been called
        mock_load.assert_called_once()


class TestMungeSumstatsConfigHandling:
    """Tests for munge_sumstats config_file parameter."""

    @patch(f"{_MUNGING}.read_config")
    @patch(f"{_MUNGING}.load_and_munge")
    def test_config_file_loaded_when_provided(
        self, mock_load, mock_read_config, tmp_path, valid_sumstats_file
    ):
        """Config file should be loaded when provided and exists."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        mock_config = {"test_sumstats": {"column_mapping": {"CHROM": "CHR"}}}
        mock_read_config.return_value = mock_config

        config_file = tmp_path / "config.json"
        with open(config_file, "w") as f:
            json.dump(mock_config, f)

        output_dir = str(tmp_path / "output")
        result = munge_sumstats(valid_sumstats_file, output_dir, config_file=str(config_file))

        mock_read_config.assert_called_once_with(str(config_file))


class TestMungeSumstatsReturnValue:
    """Tests for munge_sumstats return value."""

    @patch(f"{_MUNGING}.load_and_munge")
    def test_returns_dict_mapping_identifier_to_path(
        self, mock_load, tmp_path, valid_sumstats_file
    ):
        """Should return dict mapping identifier to output path."""
        mock_load.return_value = pd.DataFrame({"CHR": [1], "BP": [1000]})
        output_dir = str(tmp_path / "output")
        input_dict = {"EUR": valid_sumstats_file, "ASN": valid_sumstats_file}

        result = munge_sumstats(input_dict, output_dir)

        assert isinstance(result, dict)
        assert "EUR" in result
        assert "ASN" in result
        assert result["EUR"].endswith("EUR.munged.txt.gz")
        assert result["ASN"].endswith("ASN.munged.txt.gz")


class TestMungeSumstatsErrorHandling:
    """Tests for munge_sumstats error handling during processing."""

    @patch(f"{_MUNGING}.load_and_munge")
    def test_processing_failure_cleans_up_output_file(
        self, mock_load, tmp_path, valid_sumstats_file
    ):
        """Processing failure should clean up partially created output file."""
        mock_load.side_effect = ValueError("Munging failed")
        output_dir = str(tmp_path / "output")
        os.makedirs(output_dir, exist_ok=True)

        with pytest.raises(ValueError, match="Munging failed"):
            munge_sumstats(valid_sumstats_file, output_dir)

        # Check that no output file remains
        output_file = os.path.join(output_dir, "test_sumstats.munged.txt.gz")
        assert not os.path.exists(output_file)


class TestCreateMungeConfig:
    """Tests for create_munge_config function."""

    @patch(f"{_MUNGING}.create_config_template")
    @patch(f"{_MUNGING}.inspect_headers")
    def test_non_interactive_creates_json_config(
        self, mock_inspect, mock_create_template, tmp_path, valid_sumstats_file
    ):
        """Non-interactive mode should create JSON config file."""
        mock_inspect.return_value = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
        mock_create_template.return_value = {
            "column_mapping": {"CHR": "CHR", "BP": "BP"}
        }

        output_config = str(tmp_path / "config.json")
        sample_files = {"test": valid_sumstats_file}

        create_munge_config(sample_files, output_config, interactive=False)

        assert os.path.exists(output_config)
        with open(output_config, "r") as f:
            config = json.load(f)
        assert "test" in config
        assert "column_mapping" in config["test"]

    @patch(f"{_MUNGING}.create_config_template")
    @patch(f"{_MUNGING}.inspect_headers")
    def test_config_contains_entries_for_each_input(
        self, mock_inspect, mock_create_template, tmp_path, valid_sumstats_file
    ):
        """Config should contain entries for each input file."""
        mock_inspect.return_value = ["CHR", "BP", "EA"]
        mock_create_template.return_value = {"column_mapping": {}}

        file2 = tmp_path / "file2.tsv"
        pd.DataFrame({"CHR": [1]}).to_csv(file2, sep="\t", index=False)

        output_config = str(tmp_path / "config.json")
        sample_files = {"EUR": valid_sumstats_file, "ASN": str(file2)}

        create_munge_config(sample_files, output_config, interactive=False)

        with open(output_config, "r") as f:
            config = json.load(f)
        assert "EUR" in config
        assert "ASN" in config


class TestValidateMungedFiles:
    """Tests for validate_munged_files function."""

    def test_valid_munged_file_passes_validation(self, munged_sumstats_file):
        """Valid munged file should pass validation."""
        munged_files = {"test": munged_sumstats_file}

        result = validate_munged_files(munged_files)

        assert result["test"]["validation_passed"] is True
        assert result["test"]["file_exists"] is True
        assert len(result["test"]["missing_columns"]) == 0
        assert result["test"]["n_variants"] == 3

    def test_file_not_found_returns_false(self, tmp_path):
        """Non-existent file should return file_exists=False."""
        non_existent = str(tmp_path / "nonexistent.txt.gz")
        munged_files = {"test": non_existent}

        result = validate_munged_files(munged_files)

        assert result["test"]["file_exists"] is False
        assert result["test"]["validation_passed"] is False

    def test_missing_columns_populated_in_result(self, tmp_path):
        """Missing required columns should be listed in missing_columns."""
        # Create file with only some required columns
        df = pd.DataFrame(
            {
                "CHR": [1, 2],
                "BP": [1000, 2000],
                "EA": ["A", "C"],
                "NEA": ["G", "T"],
                # Missing SNPID, EAF, BETA, SE, P, RSID
            }
        )
        filepath = tmp_path / "incomplete.txt.gz"
        df.to_csv(filepath, sep="\t", index=False, compression="gzip")

        munged_files = {"test": str(filepath)}
        result = validate_munged_files(munged_files)

        assert result["test"]["validation_passed"] is False
        assert "SNPID" in result["test"]["missing_columns"]
        assert "BETA" in result["test"]["missing_columns"]
        assert "SE" in result["test"]["missing_columns"]
        assert "P" in result["test"]["missing_columns"]

    def test_n_variants_correct_count(self, munged_sumstats_file):
        """n_variants should contain correct variant count."""
        munged_files = {"test": munged_sumstats_file}

        result = validate_munged_files(munged_files)

        assert result["test"]["n_variants"] == 3

    def test_corrupted_file_has_error_field(self, tmp_path):
        """Corrupted/invalid file should populate error field."""
        # Create a file that is NOT valid gzip (raw bytes with .gz extension)
        filepath = tmp_path / "corrupted.txt.gz"
        filepath.write_bytes(b"\x00\x01\x02\x03\x04\x05")

        munged_files = {"test": str(filepath)}
        result = validate_munged_files(munged_files)

        assert "error" in result["test"]
        assert result["test"]["validation_passed"] is False
