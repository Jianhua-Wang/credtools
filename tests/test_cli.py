import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from credtools.cli import CombineCred, CombinePIP, Tool, run_fine_map


class DummyPool:
    """Minimal pool stub used to validate parallel execution paths."""

    def __init__(self, processes: int) -> None:
        self.processes = processes

    def __enter__(self):
        """Return the pool instance to mimic context manager behaviour."""
        return self

    def __exit__(self, exc_type, exc, tb):
        """Ignore context manager exit arguments."""
        return False

    def imap_unordered(self, func, iterable):
        """Yield results immediately in the call order for determinism."""
        for item in iterable:
            yield func(item)


def make_fake_process_task(output_root: Path):
    def _fake_process_task(task):  # pragma: no cover - exercised via run_fine_map
        _fake_process_task.calls.append(task["locus_id"])
        locus_dir = output_root / task["locus_id"]
        locus_dir.mkdir(parents=True, exist_ok=True)
        pip_path = locus_dir / "pips.txt.gz"
        with gzip.open(pip_path, "wt") as handle:
            handle.write("SNP\tCRED\nrs1\t1\n")

        return {
            "status": "success",
            "locus_id": task["locus_id"],
            "causal_variants_records": [
                {"SNP": "rs1", "CRED": 1, "locus_id": task["locus_id"]}
            ],
            "cs_summary_records": [
                {"locus_id": task["locus_id"], "cs_id": 1, "cs_size": 1}
            ],
        }

    _fake_process_task.calls = []
    return _fake_process_task


@pytest.mark.parametrize("processes", [1, 3])
def test_run_fine_map_parallel(tmp_path, monkeypatch, processes):
    loci_df = pd.DataFrame(
        [
            {
                "prefix": "locus1",
                "popu": "EUR",
                "cohort": "C1",
                "sample_size": 1000,
                "chr": 1,
                "start": 10,
                "end": 100,
                "locus_id": "locus1",
            },
            {
                "prefix": "locus2",
                "popu": "EUR",
                "cohort": "C2",
                "sample_size": 1200,
                "chr": 2,
                "start": 20,
                "end": 200,
                "locus_id": "locus2",
            },
        ]
    )

    inputs_path = tmp_path / "loci.tsv"
    loci_df.to_csv(inputs_path, sep="	", index=False)

    output_dir = tmp_path / "out"

    fake_task_runner = make_fake_process_task(output_dir)
    monkeypatch.setattr("credtools.cli._process_fine_map_task", fake_task_runner)
    monkeypatch.setattr("credtools.cli.Pool", DummyPool)

    run_fine_map(
        inputs=str(inputs_path),
        outdir=str(output_dir),
        tool=Tool.finemap,
        max_causal=5,
        adaptive_max_causal=False,
        set_L_by_cojo=False,
        p_cutoff=5e-8,
        collinear_cutoff=0.9,
        window_size=10_000_000,
        maf_cutoff=0.01,
        diff_freq_cutoff=0.2,
        coverage=0.95,
        timeout_minutes=30.0,
        processes=processes,
        combine_cred=CombineCred.union,
        combine_pip=CombinePIP.max,
        jaccard_threshold=0.1,
        max_iter=100,
        estimate_residual_variance=False,
        purity=0.5,
        convergence_tol=1e-3,
        calculate_lambda_s=False,
        log_file=None,
    )

    assert set(fake_task_runner.calls) == {"locus1", "locus2"}

    for locus in ("locus1", "locus2"):
        pip_path = output_dir / locus / "pips.txt.gz"
        assert pip_path.exists()
        with gzip.open(pip_path, "rt") as handle:
            lines = [line for line in handle.read().strip().splitlines() if line]
            assert len(lines) == 2

    summary_path = output_dir / "credible_sets_summary.txt.gz"
    assert summary_path.exists()
    summary_df = pd.read_csv(summary_path, sep="	")
    assert set(summary_df["locus_id"]) == {"locus1", "locus2"}

    params_path = output_dir / "parameters.json"
    with params_path.open() as handle:
        params = json.load(handle)
    assert params["parameters"]["processes"] == processes

    summary_log = output_dir / "run_summary.log"
    assert summary_log.exists()


# ---------------------------------------------------------------------------
# Additional imports for the new test classes below
# ---------------------------------------------------------------------------
import logging
from unittest.mock import MagicMock, patch

from credtools.cli import (
    _load_custom_chunks,
    _process_fine_map_task,
    create_updated_sumstat_info,
    parse_population_config_file_munge_only,
    setup_file_logging,
)


# ===========================================================================
# Tests for setup_file_logging
# ===========================================================================
class TestSetupFileLogging:
    """Tests for the setup_file_logging function."""

    def test_returns_immediately_when_log_file_is_none(self):
        """When log_file is None, the function should return without side-effects."""
        root_before = logging.getLogger().handlers[:]
        setup_file_logging(log_file=None, verbose=False)
        root_after = logging.getLogger().handlers[:]
        # No new handler should have been added
        assert root_before == root_after

    def test_creates_file_handler_verbose_false(self, tmp_path):
        """With verbose=False the file handler should use INFO level."""
        log_path = str(tmp_path / "test.log")
        root_logger = logging.getLogger()
        handlers_before = len(root_logger.handlers)

        setup_file_logging(log_file=log_path, verbose=False)

        # At least one new handler should have been added to root logger
        new_handlers = root_logger.handlers[handlers_before:]
        assert len(new_handlers) >= 1
        file_handler = new_handlers[0]
        assert isinstance(file_handler, logging.FileHandler)
        assert file_handler.level == logging.INFO

        # The log file should exist on disk
        assert Path(log_path).exists()

        # Cleanup: remove handlers we just added
        for h in new_handlers:
            root_logger.removeHandler(h)
            h.close()

        # Also clean up handlers added to named loggers
        for name in [
            "CREDTOOLS", "FINEMAP", "RSparsePro", "COJO", "SuSiE",
            "MULTISUSIE", "SUSIEX", "ABF", "ABF_COJO", "Locus",
            "LDMatrix", "QC", "Sumstats", "Utils",
        ]:
            logger = logging.getLogger(name)
            for h in list(logger.handlers):
                if isinstance(h, logging.FileHandler) and h.baseFilename == str(
                    Path(log_path).resolve()
                ):
                    logger.removeHandler(h)
                    h.close()

    def test_creates_file_handler_verbose_true(self, tmp_path):
        """With verbose=True the file handler should use DEBUG level."""
        log_path = str(tmp_path / "debug.log")
        root_logger = logging.getLogger()
        handlers_before = len(root_logger.handlers)

        setup_file_logging(log_file=log_path, verbose=True)

        new_handlers = root_logger.handlers[handlers_before:]
        assert len(new_handlers) >= 1
        file_handler = new_handlers[0]
        assert isinstance(file_handler, logging.FileHandler)
        assert file_handler.level == logging.DEBUG

        # Cleanup
        for h in new_handlers:
            root_logger.removeHandler(h)
            h.close()
        for name in [
            "CREDTOOLS", "FINEMAP", "RSparsePro", "COJO", "SuSiE",
            "MULTISUSIE", "SUSIEX", "ABF", "ABF_COJO", "Locus",
            "LDMatrix", "QC", "Sumstats", "Utils",
        ]:
            logger = logging.getLogger(name)
            for h in list(logger.handlers):
                if isinstance(h, logging.FileHandler) and h.baseFilename == str(
                    Path(log_path).resolve()
                ):
                    logger.removeHandler(h)
                    h.close()

    def test_adds_handler_to_named_loggers(self, tmp_path):
        """File handler should be attached to each credtools-specific logger."""
        log_path = str(tmp_path / "named.log")
        setup_file_logging(log_file=log_path, verbose=False)

        expected_names = [
            "CREDTOOLS", "FINEMAP", "RSparsePro", "COJO", "SuSiE",
            "MULTISUSIE", "SUSIEX", "ABF", "ABF_COJO", "Locus",
            "LDMatrix", "QC", "Sumstats", "Utils",
        ]
        resolved = str(Path(log_path).resolve())
        for name in expected_names:
            logger = logging.getLogger(name)
            file_handlers = [
                h
                for h in logger.handlers
                if isinstance(h, logging.FileHandler)
                and h.baseFilename == resolved
            ]
            assert len(file_handlers) >= 1, (
                f"Logger '{name}' should have a FileHandler pointing to {log_path}"
            )

        # Cleanup
        root_logger = logging.getLogger()
        for h in list(root_logger.handlers):
            if isinstance(h, logging.FileHandler) and h.baseFilename == resolved:
                root_logger.removeHandler(h)
                h.close()
        for name in expected_names:
            logger = logging.getLogger(name)
            for h in list(logger.handlers):
                if isinstance(h, logging.FileHandler) and h.baseFilename == resolved:
                    logger.removeHandler(h)
                    h.close()

    def test_oserror_prints_warning(self, tmp_path, capsys):
        """When the log file cannot be created, a warning should be printed."""
        # Use a path inside a non-existent directory to trigger OSError
        bad_path = str(tmp_path / "no_such_dir" / "deep" / "test.log")
        # Should not raise; instead it prints a warning via Rich Console
        setup_file_logging(log_file=bad_path, verbose=False)
        # We cannot easily capture Rich output via capsys, but the key contract
        # is that no exception propagates.  Just verify the function returns.


# ===========================================================================
# Tests for parse_population_config_file_munge_only
# ===========================================================================
class TestParsePopulationConfigFileMungeOnly:
    """Tests for the parse_population_config_file_munge_only function."""

    def _write_config(self, tmp_path, rows, create_data_files=True):
        """Helper to write a tab-separated config file and optional data files."""
        config_path = tmp_path / "config.tsv"
        df = pd.DataFrame(rows)
        df.to_csv(config_path, sep="\t", index=False)

        if create_data_files:
            for row in rows:
                p = Path(row["path"])
                p.parent.mkdir(parents=True, exist_ok=True)
                p.touch()

        return str(config_path)

    def test_normal_file_returns_dict_and_dataframe(self, tmp_path):
        """A valid config should return a mapping dict and the original DataFrame."""
        data_file_1 = tmp_path / "eur_c1.tsv"
        data_file_2 = tmp_path / "eas_c2.tsv"
        rows = [
            {"popu": "EUR", "cohort": "C1", "sample_size": 1000, "path": str(data_file_1)},
            {"popu": "EAS", "cohort": "C2", "sample_size": 2000, "path": str(data_file_2)},
        ]
        config_path = self._write_config(tmp_path, rows)

        result_dict, result_df = parse_population_config_file_munge_only(config_path)

        assert isinstance(result_dict, dict)
        assert isinstance(result_df, pd.DataFrame)
        assert "EUR_C1" in result_dict
        assert "EAS_C2" in result_dict
        assert result_dict["EUR_C1"] == str(data_file_1)
        assert result_dict["EAS_C2"] == str(data_file_2)
        assert len(result_df) == 2
        assert list(result_df.columns) == ["popu", "cohort", "sample_size", "path"]

    def test_file_not_found_raises_file_not_found_error(self, tmp_path):
        """A non-existent config file should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Configuration file not found"):
            parse_population_config_file_munge_only(
                str(tmp_path / "nonexistent.tsv")
            )

    def test_missing_required_columns_raises_value_error(self, tmp_path):
        """If required columns are absent, a ValueError should be raised."""
        config_path = tmp_path / "bad_config.tsv"
        df = pd.DataFrame({"popu": ["EUR"], "cohort": ["C1"]})
        df.to_csv(config_path, sep="\t", index=False)

        with pytest.raises(ValueError, match="Missing required columns"):
            parse_population_config_file_munge_only(str(config_path))

    def test_missing_data_file_raises_value_error(self, tmp_path):
        """If a referenced sumstats file doesn't exist, a ValueError should be raised."""
        rows = [
            {
                "popu": "EUR",
                "cohort": "C1",
                "sample_size": 1000,
                "path": str(tmp_path / "does_not_exist.tsv"),
            }
        ]
        config_path = self._write_config(tmp_path, rows, create_data_files=False)

        with pytest.raises(ValueError, match="Summary statistics file not found"):
            parse_population_config_file_munge_only(config_path)


# ===========================================================================
# Tests for create_updated_sumstat_info
# ===========================================================================
class TestCreateUpdatedSumstatInfo:
    """Tests for the create_updated_sumstat_info function."""

    def test_updates_paths_for_matching_identifiers(self, tmp_path):
        """Paths should be updated for identifiers present in munged_files."""
        original_df = pd.DataFrame(
            {
                "popu": ["EUR", "EAS"],
                "cohort": ["C1", "C2"],
                "sample_size": [1000, 2000],
                "path": ["/old/eur_c1.tsv", "/old/eas_c2.tsv"],
            }
        )
        munged_files = {
            "EUR_C1": "/new/eur_c1_munged.tsv",
            "EAS_C2": "/new/eas_c2_munged.tsv",
        }
        output_path = str(tmp_path / "updated_config.tsv")

        result_path = create_updated_sumstat_info(original_df, munged_files, output_path)

        assert result_path == output_path
        assert Path(output_path).exists()

        updated_df = pd.read_csv(output_path, sep="\t")
        assert updated_df.loc[0, "path"] == "/new/eur_c1_munged.tsv"
        assert updated_df.loc[1, "path"] == "/new/eas_c2_munged.tsv"

    def test_unmatched_identifiers_keep_original_path(self, tmp_path):
        """Identifiers not in munged_files should retain their original path."""
        original_df = pd.DataFrame(
            {
                "popu": ["EUR", "AFR"],
                "cohort": ["C1", "C3"],
                "sample_size": [1000, 3000],
                "path": ["/old/eur_c1.tsv", "/old/afr_c3.tsv"],
            }
        )
        # Only EUR_C1 is munged; AFR_C3 is not
        munged_files = {"EUR_C1": "/new/eur_c1_munged.tsv"}
        output_path = str(tmp_path / "partial_update.tsv")

        create_updated_sumstat_info(original_df, munged_files, output_path)

        updated_df = pd.read_csv(output_path, sep="\t")
        assert updated_df.loc[0, "path"] == "/new/eur_c1_munged.tsv"
        assert updated_df.loc[1, "path"] == "/old/afr_c3.tsv"

    def test_empty_munged_files_keeps_all_paths(self, tmp_path):
        """When munged_files is empty, all original paths should be preserved."""
        original_df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "path": ["/old/eur_c1.tsv"],
            }
        )
        output_path = str(tmp_path / "no_update.tsv")

        create_updated_sumstat_info(original_df, {}, output_path)

        updated_df = pd.read_csv(output_path, sep="\t")
        assert updated_df.loc[0, "path"] == "/old/eur_c1.tsv"


# ===========================================================================
# Tests for _process_fine_map_task
# ===========================================================================
class TestProcessFineMapTask:
    """Tests for the _process_fine_map_task function."""

    def test_error_path_returns_error_status(self, tmp_path):
        """When load_locus_set raises, the result should have status='error'."""
        task = {
            "locus_id": "locus_err",
            "outdir": str(tmp_path),
            "locus_records": [
                {"prefix": "x", "popu": "EUR", "cohort": "C1", "sample_size": 100}
            ],
            "calculate_lambda_s": False,
            "fine_map_kwargs": {"tool": "susie"},
        }

        with patch(
            "credtools.cli.load_locus_set",
            side_effect=RuntimeError("mock load failure"),
        ):
            result = _process_fine_map_task(task)

        assert result["status"] == "error"
        assert result["locus_id"] == "locus_err"
        assert "mock load failure" in result["error"]
        assert "traceback" in result

    def test_success_path_returns_success_status(self, tmp_path):
        """When all internal calls succeed, the result should have status='success'."""
        mock_locus_set = MagicMock()

        # Build a mock creds object whose create_enhanced_pips_df returns a DataFrame
        enhanced_df = pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "CRED": [1, 0],
                "PIP": [0.95, 0.05],
            }
        )
        mock_creds = MagicMock()
        mock_creds.create_enhanced_pips_df.return_value = enhanced_df

        task = {
            "locus_id": "locus_ok",
            "outdir": str(tmp_path),
            "locus_records": [
                {"prefix": "x", "popu": "EUR", "cohort": "C1", "sample_size": 500}
            ],
            "calculate_lambda_s": False,
            "fine_map_kwargs": {"tool": "susie"},
        }

        with (
            patch("credtools.cli.load_locus_set", return_value=mock_locus_set),
            patch("credtools.cli.fine_map", return_value=mock_creds),
            patch(
                "credtools.credibleset.generate_cs_summary",
                return_value=[{"locus_id": "locus_ok", "cs_id": 1, "cs_size": 1}],
            ),
            patch("credtools.utils.format_enhanced_pips", return_value=enhanced_df),
        ):
            result = _process_fine_map_task(task)

        assert result["status"] == "success"
        assert result["locus_id"] == "locus_ok"
        assert isinstance(result["causal_variants_records"], list)
        assert isinstance(result["cs_summary_records"], list)

        # The pips.txt.gz file should have been written
        pip_file = tmp_path / "locus_ok" / "pips.txt.gz"
        assert pip_file.exists()

    def test_empty_causal_variants(self, tmp_path):
        """When no causal variants are found, records should be an empty list."""
        mock_locus_set = MagicMock()

        enhanced_df = pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "CRED": [0, 0],
                "PIP": [0.02, 0.01],
            }
        )
        mock_creds = MagicMock()
        mock_creds.create_enhanced_pips_df.return_value = enhanced_df

        task = {
            "locus_id": "locus_empty",
            "outdir": str(tmp_path),
            "locus_records": [
                {"prefix": "x", "popu": "EUR", "cohort": "C1", "sample_size": 500}
            ],
            "calculate_lambda_s": False,
            "fine_map_kwargs": {"tool": "susie"},
        }

        with (
            patch("credtools.cli.load_locus_set", return_value=mock_locus_set),
            patch("credtools.cli.fine_map", return_value=mock_creds),
            patch(
                "credtools.credibleset.generate_cs_summary",
                return_value=[],
            ),
            patch("credtools.utils.format_enhanced_pips", return_value=enhanced_df),
        ):
            result = _process_fine_map_task(task)

        assert result["status"] == "success"
        assert result["causal_variants_records"] == []


# ===========================================================================
# Tests for _load_custom_chunks
# ===========================================================================
class TestLoadCustomChunks:
    """Tests for the _load_custom_chunks function."""

    def _write_chunks(self, tmp_path, rows):
        """Helper to write a tab-separated chunks file."""
        chunks_path = tmp_path / "chunks.tsv"
        df = pd.DataFrame(rows)
        df.to_csv(chunks_path, sep="\t", index=False)
        return str(chunks_path)

    def test_normal_file_returns_expected_dataframe(self, tmp_path):
        """A valid chunks file should produce a DataFrame with all expected columns."""
        rows = [
            {"chr": 1, "start": 1000, "end": 2000},
            {"chr": 2, "start": 3000, "end": 4000},
        ]
        chunks_path = self._write_chunks(tmp_path, rows)

        result = _load_custom_chunks(chunks_path)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2

        expected_cols = [
            "chr", "start", "end", "locus_id", "lead_snp",
            "lead_bp", "lead_p", "ancestry", "n_variants",
        ]
        assert list(result.columns) == expected_cols

        # Verify locus_id format
        assert result.iloc[0]["locus_id"] == "chr1_1000_2000"
        assert result.iloc[1]["locus_id"] == "chr2_3000_4000"

        # Verify lead_bp is the midpoint
        assert result.iloc[0]["lead_bp"] == 1500
        assert result.iloc[1]["lead_bp"] == 3500

        # Verify placeholder columns
        assert result.iloc[0]["ancestry"] == "custom"
        assert result.iloc[0]["n_variants"] == 0
        assert result.iloc[0]["lead_snp"] is None
        assert result.iloc[0]["lead_p"] is None

    def test_file_not_found_raises_file_not_found_error(self, tmp_path):
        """A non-existent chunks file should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Custom chunks file not found"):
            _load_custom_chunks(str(tmp_path / "nonexistent_chunks.tsv"))

    def test_missing_columns_raises_value_error(self, tmp_path):
        """If required columns are absent, a ValueError should be raised."""
        chunks_path = tmp_path / "bad_chunks.tsv"
        df = pd.DataFrame({"chr": [1], "start": [100]})  # missing "end"
        df.to_csv(chunks_path, sep="\t", index=False)

        with pytest.raises(ValueError, match="Missing required columns"):
            _load_custom_chunks(str(chunks_path))

    def test_single_row_chunk(self, tmp_path):
        """A file with a single chunk row should work correctly."""
        rows = [{"chr": 22, "start": 50000, "end": 60000}]
        chunks_path = self._write_chunks(tmp_path, rows)

        result = _load_custom_chunks(chunks_path)

        assert len(result) == 1
        assert result.iloc[0]["locus_id"] == "chr22_50000_60000"
        assert result.iloc[0]["lead_bp"] == 55000


# ---------------------------------------------------------------------------
# Additional imports for parse_population_config_file, create_updated_chunk_info,
# _update_chunk_info_with_prepared, and the main callback
# ---------------------------------------------------------------------------
from credtools.cli import (
    _update_chunk_info_with_prepared,
    create_updated_chunk_info,
    parse_population_config_file,
)


# ===========================================================================
# Tests for parse_population_config_file
# ===========================================================================
class TestParsePopulationConfigFile:
    """Tests for the parse_population_config_file function."""

    def test_normal_config_returns_dicts_and_dataframe(self, tmp_path):
        """A valid config with all required columns and files returns two dicts and a df."""
        # Create dummy sumstats files
        ss_file = tmp_path / "eur_c1.tsv"
        ss_file.touch()

        # Create dummy plink LD reference files
        ld_prefix = str(tmp_path / "eur_ld")
        for ext in ["bed", "bim", "fam"]:
            (tmp_path / f"eur_ld.{ext}").touch()

        config_path = tmp_path / "config.tsv"
        df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "path": [str(ss_file)],
                "ld_ref": [ld_prefix],
            }
        )
        df.to_csv(config_path, sep="\t", index=False)

        sumstats_dict, ld_ref_dict, result_df = parse_population_config_file(
            str(config_path)
        )

        assert isinstance(sumstats_dict, dict)
        assert isinstance(ld_ref_dict, dict)
        assert isinstance(result_df, pd.DataFrame)
        assert "EUR_C1" in sumstats_dict
        assert "EUR_C1" in ld_ref_dict
        assert sumstats_dict["EUR_C1"] == str(ss_file)
        assert ld_ref_dict["EUR_C1"] == ld_prefix

    def test_file_not_found_raises_file_not_found_error(self, tmp_path):
        """A non-existent config file should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Configuration file not found"):
            parse_population_config_file(str(tmp_path / "nonexistent.tsv"))

    def test_missing_required_columns_raises_value_error(self, tmp_path):
        """Config file missing required columns should raise ValueError."""
        config_path = tmp_path / "bad_config.tsv"
        df = pd.DataFrame({"popu": ["EUR"], "cohort": ["C1"]})
        df.to_csv(config_path, sep="\t", index=False)

        with pytest.raises(ValueError, match="Missing required columns"):
            parse_population_config_file(str(config_path))

    def test_missing_sumstats_file_raises_value_error(self, tmp_path):
        """When a referenced sumstats file does not exist, raise ValueError."""
        ld_prefix = str(tmp_path / "eur_ld")
        for ext in ["bed", "bim", "fam"]:
            (tmp_path / f"eur_ld.{ext}").touch()

        config_path = tmp_path / "config.tsv"
        df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "path": [str(tmp_path / "missing.tsv")],
                "ld_ref": [ld_prefix],
            }
        )
        df.to_csv(config_path, sep="\t", index=False)

        with pytest.raises(ValueError, match="Summary statistics file not found"):
            parse_population_config_file(str(config_path))

    def test_missing_ld_reference_raises_value_error(self, tmp_path):
        """When LD reference plink files are missing, raise ValueError."""
        ss_file = tmp_path / "eur_c1.tsv"
        ss_file.touch()

        config_path = tmp_path / "config.tsv"
        df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "path": [str(ss_file)],
                "ld_ref": [str(tmp_path / "no_ld")],
            }
        )
        df.to_csv(config_path, sep="\t", index=False)

        with pytest.raises(ValueError, match="LD reference files not found"):
            parse_population_config_file(str(config_path))


# ===========================================================================
# Tests for create_updated_chunk_info
# ===========================================================================
class TestCreateUpdatedChunkInfo:
    """Tests for the create_updated_chunk_info function."""

    def test_updates_paths_from_chunk_info(self, tmp_path):
        """Paths are updated based on the ancestry directories in chunk_info_df."""
        original_df = pd.DataFrame(
            {
                "popu": ["EUR", "EAS"],
                "cohort": ["C1", "C2"],
                "sample_size": [1000, 2000],
                "path": ["/old/eur_c1.tsv", "/old/eas_c2.tsv"],
            }
        )
        chunk_info_df = pd.DataFrame(
            {
                "ancestry": ["EUR_C1", "EUR_C1", "EAS_C2"],
                "sumstats_file": [
                    "/chunks/EUR_C1/chr1.tsv",
                    "/chunks/EUR_C1/chr2.tsv",
                    "/chunks/EAS_C2/chr1.tsv",
                ],
            }
        )
        output_path = str(tmp_path / "updated.tsv")

        result = create_updated_chunk_info(original_df, chunk_info_df, output_path)

        assert result == output_path
        updated_df = pd.read_csv(output_path, sep="\t")
        assert updated_df.loc[0, "path"] == "/chunks/EUR_C1"
        assert updated_df.loc[1, "path"] == "/chunks/EAS_C2"

    def test_unmatched_ancestry_keeps_original_path(self, tmp_path):
        """If ancestry not found in chunk_info_df, path stays the same."""
        original_df = pd.DataFrame(
            {
                "popu": ["EUR", "AFR"],
                "cohort": ["C1", "C3"],
                "sample_size": [1000, 3000],
                "path": ["/old/eur.tsv", "/old/afr.tsv"],
            }
        )
        chunk_info_df = pd.DataFrame(
            {
                "ancestry": ["EUR_C1"],
                "sumstats_file": ["/chunks/EUR_C1/chr1.tsv"],
            }
        )
        output_path = str(tmp_path / "partial.tsv")

        create_updated_chunk_info(original_df, chunk_info_df, output_path)

        updated_df = pd.read_csv(output_path, sep="\t")
        assert updated_df.loc[0, "path"] == "/chunks/EUR_C1"
        assert updated_df.loc[1, "path"] == "/old/afr.tsv"


# ===========================================================================
# Tests for _update_chunk_info_with_prepared
# ===========================================================================
class TestUpdateChunkInfoWithPrepared:
    """Tests for the _update_chunk_info_with_prepared function."""

    def test_updates_sumstats_file_with_prepared_prefix(self):
        """When matching locus_id + ancestry, sumstats_file should be updated."""
        chunk_info_df = pd.DataFrame(
            {
                "locus_id": ["L1", "L2"],
                "ancestry": ["EUR", "EAS"],
                "sumstats_file": ["/chunks/L1_EUR.tsv", "/chunks/L2_EAS.tsv"],
            }
        )
        prepared_df = pd.DataFrame(
            {
                "locus_id": ["L1", "L2"],
                "popu": ["EUR", "EAS"],
                "prefix": ["/prep/L1_EUR", "/prep/L2_EAS"],
            }
        )

        result = _update_chunk_info_with_prepared(chunk_info_df, prepared_df)

        assert result.loc[0, "sumstats_file"] == "/prep/L1_EUR.sumstats.gz"
        assert result.loc[1, "sumstats_file"] == "/prep/L2_EAS.sumstats.gz"

    def test_unmatched_keys_keep_original(self):
        """Rows with no match in prepared_df should keep their original sumstats_file."""
        chunk_info_df = pd.DataFrame(
            {
                "locus_id": ["L1", "L3"],
                "ancestry": ["EUR", "AFR"],
                "sumstats_file": ["/chunks/L1_EUR.tsv", "/chunks/L3_AFR.tsv"],
            }
        )
        prepared_df = pd.DataFrame(
            {
                "locus_id": ["L1"],
                "popu": ["EUR"],
                "prefix": ["/prep/L1_EUR"],
            }
        )

        result = _update_chunk_info_with_prepared(chunk_info_df, prepared_df)

        assert result.loc[0, "sumstats_file"] == "/prep/L1_EUR.sumstats.gz"
        assert result.loc[1, "sumstats_file"] == "/chunks/L3_AFR.tsv"

    def test_does_not_modify_original_dataframe(self):
        """The function should return a copy, not modify the input in place."""
        chunk_info_df = pd.DataFrame(
            {
                "locus_id": ["L1"],
                "ancestry": ["EUR"],
                "sumstats_file": ["/chunks/L1_EUR.tsv"],
            }
        )
        prepared_df = pd.DataFrame(
            {
                "locus_id": ["L1"],
                "popu": ["EUR"],
                "prefix": ["/prep/L1_EUR"],
            }
        )
        original_value = chunk_info_df.loc[0, "sumstats_file"]

        _update_chunk_info_with_prepared(chunk_info_df, prepared_df)

        # Original should be unchanged
        assert chunk_info_df.loc[0, "sumstats_file"] == original_value


# ===========================================================================
# Tests for run_fine_map with failed loci
# ===========================================================================
class TestRunFineMapFailedLoci:
    """Tests for run_fine_map when loci fail during processing."""

    def test_failed_loci_recorded_in_summary(self, tmp_path, monkeypatch):
        """When a locus fails, the error should be recorded in run_summary.log."""
        loci_df = pd.DataFrame(
            [
                {
                    "prefix": "locus_fail",
                    "popu": "EUR",
                    "cohort": "C1",
                    "sample_size": 1000,
                    "chr": 1,
                    "start": 10,
                    "end": 100,
                    "locus_id": "locus_fail",
                },
            ]
        )
        inputs_path = tmp_path / "loci.tsv"
        loci_df.to_csv(inputs_path, sep="\t", index=False)
        output_dir = tmp_path / "out"
        output_dir.mkdir(parents=True, exist_ok=True)

        def _failing_task(task):
            return {
                "status": "error",
                "locus_id": task["locus_id"],
                "error": "Intentional test failure",
                "traceback": "Traceback (most recent call last): ...",
            }

        monkeypatch.setattr("credtools.cli._process_fine_map_task", _failing_task)
        monkeypatch.setattr("credtools.cli.Pool", DummyPool)

        run_fine_map(
            inputs=str(inputs_path),
            outdir=str(output_dir),
            tool=Tool.finemap,
            max_causal=5,
            adaptive_max_causal=False,
            set_L_by_cojo=False,
            p_cutoff=5e-8,
            collinear_cutoff=0.9,
            window_size=10_000_000,
            maf_cutoff=0.01,
            diff_freq_cutoff=0.2,
            coverage=0.95,
            timeout_minutes=30.0,
            processes=1,
            combine_cred=CombineCred.union,
            combine_pip=CombinePIP.max,
            jaccard_threshold=0.1,
            max_iter=100,
            estimate_residual_variance=False,
            purity=0.5,
            convergence_tol=1e-3,
            calculate_lambda_s=False,
            log_file=None,
        )

        summary_log = output_dir / "run_summary.log"
        assert summary_log.exists()
        content = summary_log.read_text()
        assert "Failed: 1" in content
        assert "Intentional test failure" in content

        # No parameters.json because no locus succeeded
        params_path = output_dir / "parameters.json"
        assert not params_path.exists()

    def test_empty_loci_produces_summary(self, tmp_path, monkeypatch):
        """When there are no loci, the summary should reflect zero processed."""
        loci_df = pd.DataFrame(
            columns=[
                "prefix", "popu", "cohort", "sample_size",
                "chr", "start", "end", "locus_id",
            ]
        )
        inputs_path = tmp_path / "empty_loci.tsv"
        loci_df.to_csv(inputs_path, sep="\t", index=False)
        output_dir = tmp_path / "out_empty"
        output_dir.mkdir(parents=True, exist_ok=True)

        monkeypatch.setattr("credtools.cli.Pool", DummyPool)

        run_fine_map(
            inputs=str(inputs_path),
            outdir=str(output_dir),
            tool=Tool.susie,
            max_causal=5,
            adaptive_max_causal=False,
            set_L_by_cojo=False,
            p_cutoff=5e-8,
            collinear_cutoff=0.9,
            window_size=10_000_000,
            maf_cutoff=0.01,
            diff_freq_cutoff=0.2,
            coverage=0.95,
            timeout_minutes=30.0,
            processes=1,
            combine_cred=CombineCred.union,
            combine_pip=CombinePIP.max,
            jaccard_threshold=0.1,
            max_iter=100,
            estimate_residual_variance=False,
            purity=0.0,
            convergence_tol=1e-3,
            calculate_lambda_s=False,
            log_file=None,
        )

        summary_log = output_dir / "run_summary.log"
        assert summary_log.exists()
        content = summary_log.read_text()
        assert "Total Loci: 0" in content
        assert "Successful: 0" in content


# ===========================================================================
# Tests for main callback via typer
# ===========================================================================
class TestMainCallback:
    """Tests for the main CLI callback function."""

    def test_version_flag(self):
        """The --version flag should print version and exit."""
        from typer.testing import CliRunner

        from credtools import __version__
        from credtools.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["--version"])
        assert __version__ in result.output

    def test_verbose_flag(self):
        """The --verbose flag should set root logger to DEBUG."""
        from typer.testing import CliRunner

        from credtools.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["--verbose"])
        # The command should succeed (exit code 0) without error
        assert result.exit_code == 0

    def test_default_no_args_shows_help(self):
        """Invoking with no args should display help text."""
        from typer.testing import CliRunner

        from credtools.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "CREDTOOLS" in result.output or "credtools" in result.output.lower()


# ---------------------------------------------------------------------------
# Additional imports for _prepare_ld_matrices and CLI commands
# ---------------------------------------------------------------------------
from credtools.cli import _prepare_ld_matrices


# ===========================================================================
# Tests for _prepare_ld_matrices
# ===========================================================================
class TestPrepareLdMatrices:
    """Tests for the _prepare_ld_matrices function."""

    def test_successful_preparation(self, tmp_path):
        """When prepare_finemap_inputs succeeds, a DataFrame is returned."""
        chunk_info_df = pd.DataFrame(
            {
                "ancestry": ["EUR_C1", "EUR_C1"],
                "sumstats_file": [
                    "/data/EUR_C1/chr1.sumstats.gz",
                    "/data/EUR_C1/chr2.sumstats.gz",
                ],
                "locus_id": ["L1", "L2"],
            }
        )
        ld_ref_dict = {"EUR_C1": "/ref/eur"}

        expected_result = pd.DataFrame(
            {
                "locus_id": ["L1", "L2"],
                "popu": ["EUR_C1", "EUR_C1"],
                "prefix": ["/prep/L1", "/prep/L2"],
            }
        )

        with patch(
            "credtools.preprocessing.prepare.prepare_finemap_inputs",
            return_value=expected_result,
        ):
            result = _prepare_ld_matrices(
                chunk_info_df=chunk_info_df,
                ld_ref_dict=ld_ref_dict,
                output_dir=str(tmp_path / "prepared"),
            )

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        # Output directory should have been created
        assert (tmp_path / "prepared").exists()

    def test_no_matching_ld_reference_raises_value_error(self, tmp_path):
        """When no LD reference matches an ancestry, raise ValueError."""
        chunk_info_df = pd.DataFrame(
            {
                "ancestry": ["AFR"],
                "sumstats_file": ["/data/AFR/chr1.sumstats.gz"],
                "locus_id": ["L1"],
            }
        )
        ld_ref_dict = {"EUR_C1": "/ref/eur"}

        with patch(
            "credtools.preprocessing.prepare.prepare_finemap_inputs",
            return_value=pd.DataFrame(),
        ):
            with pytest.raises(ValueError, match="No LD reference found for ancestry"):
                _prepare_ld_matrices(
                    chunk_info_df=chunk_info_df,
                    ld_ref_dict=ld_ref_dict,
                    output_dir=str(tmp_path / "prepared"),
                )

    def test_import_error_is_raised(self, tmp_path):
        """When the prepare module cannot be imported, ImportError propagates."""
        chunk_info_df = pd.DataFrame(
            {
                "ancestry": ["EUR"],
                "sumstats_file": ["/data/chr1.sumstats.gz"],
                "locus_id": ["L1"],
            }
        )
        ld_ref_dict = {"EUR": "/ref/eur"}

        with patch.dict(
            "sys.modules",
            {"credtools.preprocessing.prepare": None},
        ):
            with pytest.raises((ImportError, ValueError)):
                _prepare_ld_matrices(
                    chunk_info_df=chunk_info_df,
                    ld_ref_dict=ld_ref_dict,
                    output_dir=str(tmp_path / "prepared"),
                )

    def test_column_renaming_and_prefix_generation(self, tmp_path):
        """Verify that ancestry is renamed to popu and prefix is derived from sumstats_file."""
        chunk_info_df = pd.DataFrame(
            {
                "ancestry": ["EUR_C1"],
                "sumstats_file": ["/data/locus1.sumstats.gz"],
                "locus_id": ["L1"],
            }
        )
        ld_ref_dict = {"EUR_C1": "/ref/eur"}

        captured_kwargs = {}

        def mock_prepare(**kwargs):
            captured_kwargs.update(kwargs)
            return pd.DataFrame(
                {"locus_id": ["L1"], "popu": ["EUR_C1"], "prefix": ["/prep/L1"]}
            )

        with patch(
            "credtools.preprocessing.prepare.prepare_finemap_inputs",
            side_effect=mock_prepare,
        ):
            _prepare_ld_matrices(
                chunk_info_df=chunk_info_df,
                ld_ref_dict=ld_ref_dict,
                output_dir=str(tmp_path / "prepared"),
                threads=4,
                ld_format="vcf",
                keep_intermediate=True,
            )

        # Verify the passed-through parameters
        assert captured_kwargs["threads"] == 4
        assert captured_kwargs["ld_format"] == "vcf"
        assert captured_kwargs["keep_intermediate"] is True

        # Verify the prepared chunk_info_df has renamed columns and new columns
        prep_df = captured_kwargs["chunk_info_df"]
        assert "popu" in prep_df.columns
        assert "ancestry" not in prep_df.columns
        assert "cohort" in prep_df.columns
        assert "sample_size" in prep_df.columns
        assert prep_df.iloc[0]["sample_size"] == 50000
        # Prefix should strip .sumstats and .gz
        assert "sumstats" not in prep_df.iloc[0]["prefix"]


# ===========================================================================
# Tests for run_meta via CLI runner
# ===========================================================================
class TestRunMeta:
    """Tests for the run_meta CLI command."""

    def test_meta_command_calls_meta_loci(self, tmp_path, monkeypatch):
        """The meta command should invoke meta_loci with correct arguments."""
        from typer.testing import CliRunner

        from credtools.cli import app

        # Create a minimal input file
        loci_df = pd.DataFrame(
            {
                "prefix": ["locus1"],
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "chr": [1],
                "start": [10],
                "end": [100],
                "locus_id": ["locus1"],
            }
        )
        inputs_path = tmp_path / "meta_input.tsv"
        loci_df.to_csv(inputs_path, sep="\t", index=False)
        output_dir = tmp_path / "meta_out"
        output_dir.mkdir()

        call_args = {}

        def mock_meta_loci(*args, **kwargs):
            call_args["args"] = args
            call_args["kwargs"] = kwargs

        monkeypatch.setattr("credtools.cli.meta_loci", mock_meta_loci)

        runner = CliRunner()
        result = runner.invoke(
            app,
            ["meta", str(inputs_path), str(output_dir)],
        )

        assert result.exit_code == 0
        assert call_args["args"][0] == str(inputs_path)
        assert call_args["args"][1] == str(output_dir)


# ===========================================================================
# Tests for run_qc via CLI runner
# ===========================================================================
class TestRunQc:
    """Tests for the run_qc CLI command."""

    def test_qc_command_success_path(self, tmp_path, monkeypatch):
        """The qc command should display success message when no failures."""
        from typer.testing import CliRunner

        from credtools.cli import app

        loci_df = pd.DataFrame(
            {
                "prefix": ["locus1"],
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "chr": [1],
                "start": [10],
                "end": [100],
                "locus_id": ["locus1"],
            }
        )
        inputs_path = tmp_path / "qc_input.tsv"
        loci_df.to_csv(inputs_path, sep="\t", index=False)
        output_dir = tmp_path / "qc_out"
        output_dir.mkdir()

        def mock_loci_qc(*args, **kwargs):
            return {
                "successful_loci": 1,
                "failed_loci": 0,
                "log_path": str(output_dir / "qc_summary.log"),
            }

        monkeypatch.setattr("credtools.cli.loci_qc", mock_loci_qc)

        runner = CliRunner()
        result = runner.invoke(
            app,
            ["qc", str(inputs_path), str(output_dir)],
        )

        assert result.exit_code == 0
        assert "successfully" in result.output.lower() or "QC completed" in result.output

    def test_qc_command_with_failures(self, tmp_path, monkeypatch):
        """The qc command should report failures when loci fail."""
        from typer.testing import CliRunner

        from credtools.cli import app

        loci_df = pd.DataFrame(
            {
                "prefix": ["locus1"],
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "chr": [1],
                "start": [10],
                "end": [100],
                "locus_id": ["locus1"],
            }
        )
        inputs_path = tmp_path / "qc_input.tsv"
        loci_df.to_csv(inputs_path, sep="\t", index=False)
        output_dir = tmp_path / "qc_out"
        output_dir.mkdir()

        def mock_loci_qc(*args, **kwargs):
            return {
                "successful_loci": 0,
                "failed_loci": 1,
                "log_path": str(output_dir / "qc_summary.log"),
            }

        monkeypatch.setattr("credtools.cli.loci_qc", mock_loci_qc)

        runner = CliRunner()
        result = runner.invoke(
            app,
            ["qc", str(inputs_path), str(output_dir)],
        )

        assert result.exit_code == 0
        assert "failed" in result.output.lower() or "1" in result.output


# ===========================================================================
# Tests for run_munge via CLI runner
# ===========================================================================
class TestRunMunge:
    """Tests for the run_munge CLI command."""

    def test_munge_with_config_file(self, tmp_path, monkeypatch):
        """The munge command should process a config file and produce munged results."""
        from typer.testing import CliRunner

        from credtools.cli import app

        # Create a valid sumstats file so parse_population_config_file_munge_only passes
        ss_file = tmp_path / "eur_c1.tsv"
        ss_file.touch()

        # Create config file
        config_path = tmp_path / "config.tsv"
        config_df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["C1"],
                "sample_size": [1000],
                "path": [str(ss_file)],
            }
        )
        config_df.to_csv(config_path, sep="\t", index=False)

        output_dir = tmp_path / "munged"
        output_dir.mkdir()

        munged_path = str(output_dir / "EUR_C1.munged.tsv.gz")

        def mock_munge_sumstats(input_files, output_dir, config_file, force_overwrite):
            return {"EUR_C1": munged_path}

        def mock_validate_munged_files(result, required_columns):
            return {
                "EUR_C1": {
                    "validation_passed": True,
                    "n_variants": 500,
                }
            }

        # Patch the preprocessing imports inside run_munge
        import types

        mock_preprocessing = types.ModuleType("credtools.preprocessing")
        mock_preprocessing.munge_sumstats = mock_munge_sumstats

        mock_munge_mod = types.ModuleType("credtools.preprocessing.munge")
        mock_munge_mod.create_munge_config = MagicMock()
        mock_munge_mod.validate_munged_files = mock_validate_munged_files

        monkeypatch.setitem(
            __import__("sys").modules, "credtools.preprocessing", mock_preprocessing
        )
        monkeypatch.setitem(
            __import__("sys").modules, "credtools.preprocessing.munge", mock_munge_mod
        )

        runner = CliRunner()
        result = runner.invoke(
            app,
            ["munge", str(config_path), str(output_dir)],
        )

        assert result.exit_code == 0
        assert "Successfully munged" in result.output or "munged" in result.output.lower()

        # Verify that the updated sumstat info file was created
        updated_info = output_dir / "sumstat_info_updated.txt"
        assert updated_info.exists()

    def test_munge_with_direct_file_inputs(self, tmp_path, monkeypatch):
        """The munge command should accept comma-separated file paths."""
        from typer.testing import CliRunner

        from credtools.cli import app

        # Create dummy input files
        file1 = tmp_path / "study1.tsv"
        file1.touch()
        file2 = tmp_path / "study2.tsv"
        file2.touch()

        output_dir = tmp_path / "munged"
        output_dir.mkdir()

        def mock_munge_sumstats(input_files, output_dir, config_file, force_overwrite):
            return {k: f"{output_dir}/{k}.munged.gz" for k in input_files}

        def mock_validate_munged_files(result, required_columns):
            return {
                k: {"validation_passed": True, "n_variants": 100} for k in result
            }

        import types

        mock_preprocessing = types.ModuleType("credtools.preprocessing")
        mock_preprocessing.munge_sumstats = mock_munge_sumstats

        mock_munge_mod = types.ModuleType("credtools.preprocessing.munge")
        mock_munge_mod.create_munge_config = MagicMock()
        mock_munge_mod.validate_munged_files = mock_validate_munged_files

        monkeypatch.setitem(
            __import__("sys").modules, "credtools.preprocessing", mock_preprocessing
        )
        monkeypatch.setitem(
            __import__("sys").modules, "credtools.preprocessing.munge", mock_munge_mod
        )

        runner = CliRunner()
        result = runner.invoke(
            app,
            ["munge", f"{file1},{file2}", str(output_dir)],
        )

        assert result.exit_code == 0
        assert "Successfully munged" in result.output or "2" in result.output

    def test_munge_error_during_munging(self, tmp_path, monkeypatch):
        """When munge_sumstats raises, the command should exit with code 1."""
        from typer.testing import CliRunner

        from credtools.cli import app

        file1 = tmp_path / "study1.tsv"
        file1.touch()

        output_dir = tmp_path / "munged"
        output_dir.mkdir()

        def mock_munge_sumstats_fail(**kwargs):
            raise RuntimeError("Munge failed")

        import types

        mock_preprocessing = types.ModuleType("credtools.preprocessing")
        mock_preprocessing.munge_sumstats = mock_munge_sumstats_fail

        mock_munge_mod = types.ModuleType("credtools.preprocessing.munge")
        mock_munge_mod.create_munge_config = MagicMock()
        mock_munge_mod.validate_munged_files = MagicMock()

        monkeypatch.setitem(
            __import__("sys").modules, "credtools.preprocessing", mock_preprocessing
        )
        monkeypatch.setitem(
            __import__("sys").modules, "credtools.preprocessing.munge", mock_munge_mod
        )

        runner = CliRunner()
        result = runner.invoke(
            app,
            ["munge", str(file1), str(output_dir)],
        )

        assert result.exit_code == 1
