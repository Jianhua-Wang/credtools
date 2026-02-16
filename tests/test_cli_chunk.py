"""Unit tests for chunk-related helper functions in credtools.cli."""

import os

import numpy as np
import pandas as pd
import pytest

from credtools.cli import (
    _load_custom_chunks,
    _update_chunk_info_with_prepared,
    create_updated_chunk_info,
    parse_population_config_file,
)


# ---------------------------------------------------------------------------
# TestLoadCustomChunks
# ---------------------------------------------------------------------------


class TestLoadCustomChunks:
    """Tests for cli._load_custom_chunks."""

    def _write_chunks_file(self, tmp_path, rows, filename="chunks.tsv"):
        """Helper to write a TSV file and return its path."""
        df = pd.DataFrame(rows)
        fpath = tmp_path / filename
        df.to_csv(fpath, sep="\t", index=False)
        return str(fpath)

    def test_valid_file_returns_dataframe(self, tmp_path):
        """Valid file with chr/start/end returns a DataFrame."""
        fpath = self._write_chunks_file(
            tmp_path,
            [
                {"chr": 1, "start": 100, "end": 500},
                {"chr": 2, "start": 200, "end": 600},
            ],
        )
        result = _load_custom_chunks(fpath)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2

    def test_locus_id_format(self, tmp_path):
        """locus_id should follow chr{N}_{start}_{end}."""
        fpath = self._write_chunks_file(
            tmp_path, [{"chr": 3, "start": 1000, "end": 5000}]
        )
        result = _load_custom_chunks(fpath)
        assert result.iloc[0]["locus_id"] == "chr3_1000_5000"

    def test_placeholder_columns(self, tmp_path):
        """Placeholder columns: ancestry='custom', lead_snp=None, etc."""
        fpath = self._write_chunks_file(
            tmp_path, [{"chr": 1, "start": 100, "end": 500}]
        )
        result = _load_custom_chunks(fpath)
        row = result.iloc[0]
        assert row["ancestry"] == "custom"
        assert row["lead_snp"] is None
        assert row["lead_p"] is None
        assert row["n_variants"] == 0

    def test_lead_bp_midpoint(self, tmp_path):
        """lead_bp should be (start + end) // 2."""
        fpath = self._write_chunks_file(
            tmp_path, [{"chr": 1, "start": 100, "end": 500}]
        )
        result = _load_custom_chunks(fpath)
        assert result.iloc[0]["lead_bp"] == (100 + 500) // 2

    def test_output_columns_count_and_order(self, tmp_path):
        """Result should have exactly 9 columns in the expected order."""
        fpath = self._write_chunks_file(
            tmp_path, [{"chr": 1, "start": 100, "end": 500}]
        )
        result = _load_custom_chunks(fpath)
        expected_cols = [
            "chr",
            "start",
            "end",
            "locus_id",
            "lead_snp",
            "lead_bp",
            "lead_p",
            "ancestry",
            "n_variants",
        ]
        assert list(result.columns) == expected_cols

    def test_file_not_found_raises(self, tmp_path):
        """Non-existent file → FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            _load_custom_chunks(str(tmp_path / "nonexistent.tsv"))

    def test_missing_required_columns_raises(self, tmp_path):
        """File without 'end' column → ValueError."""
        fpath = self._write_chunks_file(
            tmp_path,
            [{"chr": 1, "start": 100}],  # missing 'end'
            filename="bad.tsv",
        )
        with pytest.raises(ValueError, match="Missing required columns"):
            _load_custom_chunks(fpath)

    def test_extra_columns_ignored(self, tmp_path):
        """Extra columns are not in the output."""
        fpath = self._write_chunks_file(
            tmp_path,
            [{"chr": 1, "start": 100, "end": 500, "gene": "BRCA1"}],
        )
        result = _load_custom_chunks(fpath)
        assert "gene" not in result.columns


# ---------------------------------------------------------------------------
# TestUpdateChunkInfoWithPrepared
# ---------------------------------------------------------------------------


class TestUpdateChunkInfoWithPrepared:
    """Tests for cli._update_chunk_info_with_prepared."""

    def _make_chunk_info(self, rows):
        return pd.DataFrame(rows)

    def _make_prepared(self, rows):
        return pd.DataFrame(rows)

    def test_basic_update(self):
        """Matching locus_id+ancestry → sumstats_file updated."""
        chunk_df = self._make_chunk_info(
            [
                {
                    "locus_id": "chr1_1_500",
                    "ancestry": "EUR",
                    "sumstats_file": "/old/EUR.chr1_1_500.sumstats.gz",
                }
            ]
        )
        prepared_df = self._make_prepared(
            [
                {
                    "locus_id": "chr1_1_500",
                    "popu": "EUR",
                    "prefix": "/new/EUR.chr1_1_500",
                }
            ]
        )
        result = _update_chunk_info_with_prepared(chunk_df, prepared_df)
        assert result.iloc[0]["sumstats_file"] == "/new/EUR.chr1_1_500.sumstats.gz"

    def test_no_match_unchanged(self):
        """No matching key → sumstats_file stays the same."""
        chunk_df = self._make_chunk_info(
            [
                {
                    "locus_id": "chr1_1_500",
                    "ancestry": "EUR",
                    "sumstats_file": "/old/file.gz",
                }
            ]
        )
        prepared_df = self._make_prepared(
            [
                {
                    "locus_id": "chr2_1_500",
                    "popu": "ASN",
                    "prefix": "/new/something",
                }
            ]
        )
        result = _update_chunk_info_with_prepared(chunk_df, prepared_df)
        assert result.iloc[0]["sumstats_file"] == "/old/file.gz"

    def test_partial_match(self):
        """Only matching rows are updated; others stay unchanged."""
        chunk_df = self._make_chunk_info(
            [
                {
                    "locus_id": "chr1_1_500",
                    "ancestry": "EUR",
                    "sumstats_file": "/old/eur.gz",
                },
                {
                    "locus_id": "chr2_1_500",
                    "ancestry": "ASN",
                    "sumstats_file": "/old/asn.gz",
                },
            ]
        )
        prepared_df = self._make_prepared(
            [{"locus_id": "chr1_1_500", "popu": "EUR", "prefix": "/new/eur"}]
        )
        result = _update_chunk_info_with_prepared(chunk_df, prepared_df)
        assert result.iloc[0]["sumstats_file"] == "/new/eur.sumstats.gz"
        assert result.iloc[1]["sumstats_file"] == "/old/asn.gz"

    def test_original_not_modified(self):
        """Original chunk_df should not be mutated (copy semantics)."""
        chunk_df = self._make_chunk_info(
            [
                {
                    "locus_id": "chr1_1_500",
                    "ancestry": "EUR",
                    "sumstats_file": "/old/file.gz",
                }
            ]
        )
        original_value = chunk_df.iloc[0]["sumstats_file"]
        prepared_df = self._make_prepared(
            [{"locus_id": "chr1_1_500", "popu": "EUR", "prefix": "/new/prefix"}]
        )
        _update_chunk_info_with_prepared(chunk_df, prepared_df)
        assert chunk_df.iloc[0]["sumstats_file"] == original_value

    def test_multi_ancestry_independent(self):
        """Each ancestry is matched independently by (locus_id, popu) key."""
        chunk_df = self._make_chunk_info(
            [
                {
                    "locus_id": "chr1_1_500",
                    "ancestry": "EUR",
                    "sumstats_file": "/old/eur.gz",
                },
                {
                    "locus_id": "chr1_1_500",
                    "ancestry": "ASN",
                    "sumstats_file": "/old/asn.gz",
                },
            ]
        )
        prepared_df = self._make_prepared(
            [
                {"locus_id": "chr1_1_500", "popu": "EUR", "prefix": "/new/eur"},
                {"locus_id": "chr1_1_500", "popu": "ASN", "prefix": "/new/asn"},
            ]
        )
        result = _update_chunk_info_with_prepared(chunk_df, prepared_df)
        assert result.iloc[0]["sumstats_file"] == "/new/eur.sumstats.gz"
        assert result.iloc[1]["sumstats_file"] == "/new/asn.sumstats.gz"


# ---------------------------------------------------------------------------
# TestCreateUpdatedChunkInfo
# ---------------------------------------------------------------------------


class TestCreateUpdatedChunkInfo:
    """Tests for cli.create_updated_chunk_info."""

    def test_basic_update(self, tmp_path):
        """Path column should be updated with chunk directory.

        The function matches via identifier = '{popu}_{cohort}', so the
        ancestry in chunk_info_df must be 'EUR_UKB' to match popu=EUR, cohort=UKB.
        """
        config_df = pd.DataFrame(
            [
                {
                    "popu": "EUR",
                    "cohort": "UKB",
                    "sample_size": 10000,
                    "path": "/old/eur.gz",
                    "ld_ref": "/ld/eur",
                },
            ]
        )
        chunk_info_df = pd.DataFrame(
            [
                {
                    "ancestry": "EUR_UKB",
                    "sumstats_file": "/chunks/EUR_UKB.chr1_1_500.sumstats.gz",
                },
            ]
        )
        out_path = str(tmp_path / "updated.txt")
        result_path = create_updated_chunk_info(config_df, chunk_info_df, out_path)
        assert result_path == out_path

        updated = pd.read_csv(result_path, sep="\t")
        assert updated.iloc[0]["path"] == "/chunks"

    def test_output_file_created(self, tmp_path):
        """Output file should be created on disk."""
        config_df = pd.DataFrame(
            [
                {
                    "popu": "EUR",
                    "cohort": "EUR",
                    "sample_size": 10000,
                    "path": "/old",
                    "ld_ref": "/ld",
                }
            ]
        )
        chunk_info_df = pd.DataFrame(
            [{"ancestry": "EUR", "sumstats_file": "/chunks/file.gz"}]
        )
        out_path = str(tmp_path / "updated.txt")
        create_updated_chunk_info(config_df, chunk_info_df, out_path)
        assert os.path.exists(out_path)

    def test_identifier_matching(self, tmp_path):
        """Identifier = {popu}_{cohort} must match ancestry in chunk_info."""
        config_df = pd.DataFrame(
            [
                {
                    "popu": "EUR",
                    "cohort": "UKB",
                    "sample_size": 10000,
                    "path": "/old",
                    "ld_ref": "/ld",
                },
            ]
        )
        chunk_info_df = pd.DataFrame(
            [
                {
                    "ancestry": "EUR_UKB",
                    "sumstats_file": "/chunks/EUR_UKB.chr1.sumstats.gz",
                },
            ]
        )
        out_path = str(tmp_path / "updated.txt")
        create_updated_chunk_info(config_df, chunk_info_df, out_path)
        updated = pd.read_csv(out_path, sep="\t")
        assert updated.iloc[0]["path"] == "/chunks"

    def test_returns_output_path(self, tmp_path):
        """Function should return the output_path string."""
        config_df = pd.DataFrame(
            [
                {
                    "popu": "EUR",
                    "cohort": "EUR",
                    "sample_size": 10000,
                    "path": "/old",
                    "ld_ref": "/ld",
                }
            ]
        )
        chunk_info_df = pd.DataFrame(
            [{"ancestry": "EUR", "sumstats_file": "/chunks/file.gz"}]
        )
        out_path = str(tmp_path / "result.txt")
        result = create_updated_chunk_info(config_df, chunk_info_df, out_path)
        assert result == out_path


# ---------------------------------------------------------------------------
# TestParsePopulationConfigFile
# ---------------------------------------------------------------------------


class TestParsePopulationConfigFile:
    """Tests for cli.parse_population_config_file."""

    def _setup_config(self, tmp_path, rows, create_files=True):
        """Create a config file and optionally the referenced files."""
        config_df = pd.DataFrame(rows)
        config_path = tmp_path / "config.txt"
        config_df.to_csv(config_path, sep="\t", index=False)

        if create_files:
            for _, row in config_df.iterrows():
                # Create sumstats file
                sumstats_path = row["path"]
                os.makedirs(os.path.dirname(sumstats_path), exist_ok=True)
                pd.DataFrame({"A": [1]}).to_csv(sumstats_path, sep="\t", index=False)

                # Create PLINK files
                ld_base = row["ld_ref"]
                os.makedirs(os.path.dirname(ld_base), exist_ok=True)
                for ext in ["bed", "bim", "fam"]:
                    open(f"{ld_base}.{ext}", "w").close()

        return str(config_path)

    def test_valid_config(self, tmp_path):
        """Valid config returns a three-tuple."""
        config_path = self._setup_config(
            tmp_path,
            [
                {
                    "popu": "EUR",
                    "cohort": "UKB",
                    "sample_size": 10000,
                    "path": str(tmp_path / "data" / "eur.txt"),
                    "ld_ref": str(tmp_path / "ld" / "eur"),
                },
            ],
        )
        sumstats_dict, ld_ref_dict, config_df = parse_population_config_file(
            config_path
        )
        assert "EUR_UKB" in sumstats_dict
        assert "EUR_UKB" in ld_ref_dict
        assert isinstance(config_df, pd.DataFrame)

    def test_identifier_format(self, tmp_path):
        """Identifier should be '{popu}_{cohort}'."""
        config_path = self._setup_config(
            tmp_path,
            [
                {
                    "popu": "ASN",
                    "cohort": "BBJ",
                    "sample_size": 5000,
                    "path": str(tmp_path / "data" / "asn.txt"),
                    "ld_ref": str(tmp_path / "ld" / "asn"),
                },
            ],
        )
        sumstats_dict, _, _ = parse_population_config_file(config_path)
        assert "ASN_BBJ" in sumstats_dict

    def test_file_not_found_raises(self, tmp_path):
        """Non-existent config file → FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            parse_population_config_file(str(tmp_path / "no_such_file.txt"))

    def test_missing_required_columns_raises(self, tmp_path):
        """Config file missing required columns → ValueError."""
        df = pd.DataFrame([{"popu": "EUR", "cohort": "UKB"}])
        config_path = tmp_path / "bad_config.txt"
        df.to_csv(config_path, sep="\t", index=False)
        with pytest.raises(ValueError, match="Missing required columns"):
            parse_population_config_file(str(config_path))

    def test_sumstats_file_not_found_raises(self, tmp_path):
        """Sumstats file doesn't exist → error."""
        config_path = self._setup_config(
            tmp_path,
            [
                {
                    "popu": "EUR",
                    "cohort": "UKB",
                    "sample_size": 10000,
                    "path": str(tmp_path / "nonexistent" / "eur.txt"),
                    "ld_ref": str(tmp_path / "ld" / "eur"),
                },
            ],
            create_files=False,
        )
        with pytest.raises(ValueError):
            parse_population_config_file(config_path)

    def test_ld_ref_not_found_raises(self, tmp_path):
        """LD ref files don't exist → error."""
        # Create sumstats but not LD files
        sumstats_path = str(tmp_path / "data" / "eur.txt")
        os.makedirs(os.path.dirname(sumstats_path), exist_ok=True)
        pd.DataFrame({"A": [1]}).to_csv(sumstats_path, sep="\t", index=False)

        config_path = self._setup_config(
            tmp_path,
            [
                {
                    "popu": "EUR",
                    "cohort": "UKB",
                    "sample_size": 10000,
                    "path": sumstats_path,
                    "ld_ref": str(tmp_path / "no_ld" / "missing"),
                },
            ],
            create_files=False,
        )
        with pytest.raises(ValueError):
            parse_population_config_file(config_path)
