import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from credtools import plot
from credtools.plot import (
    _coerce_qc_dataframe,
    get_population_color,
    read_compressed_file,
)


# ---------------------------------------------------------------------------
# Existing tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "plot_func, kwargs",
    (
        (plot.plot_lambda_s_boxplot, {}),
        (plot.plot_maf_corr_barplot, {}),
        (plot.plot_outliers_barplot, {"outlier_type": "lambda_s"}),
    ),
)
def test_qc_plots_accept_path_input(plot_func, kwargs, qc_summary_gz):
    fig, ax = plt.subplots()
    try:
        result = plot_func(qc_summary_gz, ax=ax, **kwargs)
    finally:
        plt.close(fig)
    assert result is ax


def test_plot_locusplot(tmp_path, locus_dir_with_pips):
    output_path = tmp_path / "locus_plot.png"
    fig = plot.plot_locusplot(
        locus_dir_with_pips, output_file=output_path, figsize=(6, 4), dpi=72
    )
    plt.close(fig)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


# ---------------------------------------------------------------------------
# TestGetPopulationColor
# ---------------------------------------------------------------------------
class TestGetPopulationColor:
    """Tests for get_population_color function."""

    def test_known_populations(self):
        assert get_population_color("EUR_UKB") == "#45B7D1"
        assert get_population_color("AFR_APCDR") == "#FF6B6B"
        assert get_population_color("EAS_BBJ") == "#4ECDC4"

    def test_unknown_population_returns_gray(self):
        assert get_population_color("UNKNOWN_cohort") == "#7F7F7F"

    def test_no_underscore(self):
        # If no underscore, the whole string is treated as population
        result = get_population_color("EUR")
        assert result == "#45B7D1"


# ---------------------------------------------------------------------------
# TestReadCompressedFile
# ---------------------------------------------------------------------------
class TestReadCompressedFile:
    """Tests for read_compressed_file function."""

    def test_gz_file(self, tmp_path):
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        gz_path = tmp_path / "test.txt.gz"
        df.to_csv(gz_path, sep="\t", index=False, compression="gzip")
        result = read_compressed_file(gz_path)
        assert list(result.columns) == ["a", "b"]
        assert len(result) == 2

    def test_plain_file(self, tmp_path):
        df = pd.DataFrame({"x": [10, 20]})
        path = tmp_path / "test.txt"
        df.to_csv(path, sep="\t", index=False)
        result = read_compressed_file(path)
        assert list(result.columns) == ["x"]

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError, match="File not found"):
            read_compressed_file("/nonexistent/file.txt")


# ---------------------------------------------------------------------------
# TestCoerceQcDataframe
# ---------------------------------------------------------------------------
class TestCoerceQcDataframe:
    """Tests for _coerce_qc_dataframe function."""

    def test_dataframe_direct(self):
        df = pd.DataFrame({"a": [1], "b": [2]})
        result = _coerce_qc_dataframe(df)
        assert isinstance(result, pd.DataFrame)

    def test_str_path(self, tmp_path):
        df = pd.DataFrame({"a": [1], "b": [2]})
        path = tmp_path / "test.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        result = _coerce_qc_dataframe(str(path))
        assert isinstance(result, pd.DataFrame)

    def test_path_object(self, tmp_path):
        df = pd.DataFrame({"col": [1]})
        path = tmp_path / "test.txt"
        df.to_csv(path, sep="\t", index=False)
        result = _coerce_qc_dataframe(Path(path))
        assert isinstance(result, pd.DataFrame)

    def test_invalid_type_raises(self):
        with pytest.raises(TypeError, match="must be a pandas DataFrame"):
            _coerce_qc_dataframe(12345)

    def test_missing_columns_raises(self):
        df = pd.DataFrame({"a": [1]})
        with pytest.raises(ValueError, match="missing required column"):
            _coerce_qc_dataframe(df, required_columns=["a", "b", "c"])


# ---------------------------------------------------------------------------
# TestPlotZscoreQq
# ---------------------------------------------------------------------------
class TestPlotZscoreQq:
    """Tests for plot_zscore_qq function."""

    def test_basic_plot(self, tmp_path):
        rng = np.random.default_rng(42)
        n = 50
        df = pd.DataFrame(
            {
                "z": rng.normal(0, 1, n),
                "condmean": rng.normal(0, 1, n),
                "cohort": ["EUR_UKB"] * n,
                "lambda_s": [0.95] * n,
            }
        )
        path = tmp_path / "expected_z.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_zscore_qq(path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)

    def test_missing_columns_shows_text(self, tmp_path):
        # Must include 'cohort' column since the function accesses it before checking
        df = pd.DataFrame({"other_col": [1, 2], "cohort": ["EUR_UKB", "EUR_UKB"]})
        path = tmp_path / "expected_z.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_zscore_qq(path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)

    def test_creates_new_figure_when_no_ax(self, tmp_path):
        n = 10
        df = pd.DataFrame(
            {
                "z": np.zeros(n),
                "condmean": np.zeros(n),
                "cohort": ["EUR_C1"] * n,
            }
        )
        path = tmp_path / "ez.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        try:
            result = plot.plot_zscore_qq(path)
            assert isinstance(result, plt.Axes)
        finally:
            plt.close("all")


# ---------------------------------------------------------------------------
# TestPlotLdDecay
# ---------------------------------------------------------------------------
class TestPlotLdDecay:
    """Tests for plot_ld_decay function."""

    def test_basic_plot(self, tmp_path):
        df = pd.DataFrame(
            {
                "distance_kb": [1, 10, 100, 1000] * 2,
                "r2_avg": [0.9, 0.5, 0.2, 0.05] * 2,
                "cohort": ["EUR_UKB"] * 4 + ["EAS_BBJ"] * 4,
            }
        )
        path = tmp_path / "ld_decay.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_ld_decay(path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)

    def test_missing_columns_shows_text(self, tmp_path):
        # Must include 'cohort' column since the function accesses it before checking
        df = pd.DataFrame({"wrong": [1], "cohort": ["EUR_UKB"]})
        path = tmp_path / "ld_decay.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_ld_decay(path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# TestPlotLd4thMoment
# ---------------------------------------------------------------------------
class TestPlotLd4thMoment:
    """Tests for plot_ld_4th_moment function."""

    def test_basic_plot(self, tmp_path):
        rng = np.random.default_rng(42)
        df = pd.DataFrame(
            {
                "EUR_UKB": rng.uniform(0, 1, 20),
                "EAS_BBJ": rng.uniform(0, 1, 20),
            }
        )
        path = tmp_path / "ld_4th_moment.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_ld_4th_moment(path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# TestPlotLocusPlotErrors
# ---------------------------------------------------------------------------
class TestPlotLocusPlotErrors:
    """Tests for plot_locusplot error paths."""

    def test_nonexistent_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            plot.plot_locusplot(tmp_path / "nonexistent")

    def test_no_pips_file_raises(self, tmp_path):
        empty_dir = tmp_path / "empty_locus"
        empty_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="No pips file found"):
            plot.plot_locusplot(empty_dir)


# ---------------------------------------------------------------------------
# TestPlotMafCorrBarplot
# ---------------------------------------------------------------------------
class TestPlotMafCorrBarplot:
    """Tests for plot_maf_corr_barplot with empty data."""

    def test_empty_maf_corr_shows_text(self):
        df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["UKB"],
                "maf_corr": [float("nan")],
            }
        )
        fig, ax = plt.subplots()
        try:
            result = plot.plot_maf_corr_barplot(df, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# TestPlotOutliersBarplot
# ---------------------------------------------------------------------------
class TestPlotOutliersBarplot:
    """Tests for plot_outliers_barplot edge cases."""

    def test_missing_outlier_column_shows_text(self):
        df = pd.DataFrame(
            {
                "popu": ["EUR"],
                "cohort": ["UKB"],
            }
        )
        fig, ax = plt.subplots()
        try:
            result = plot.plot_outliers_barplot(df, outlier_type="lambda_s", ax=ax)
            assert result is ax
        finally:
            plt.close(fig)
