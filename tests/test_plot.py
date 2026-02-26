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


# ---------------------------------------------------------------------------
# TestLoadLocusPips
# ---------------------------------------------------------------------------
class TestLoadLocusPips:
    """Tests for _load_locus_pips function."""

    def test_gz_file_returns_dataframe(self, locus_dir_with_pips):
        from credtools.plot import _load_locus_pips

        result = _load_locus_pips(locus_dir_with_pips)
        assert isinstance(result, pd.DataFrame)
        assert "BP" in result.columns
        assert "CRED" in result.columns
        assert "PIP" in result.columns
        assert len(result) > 0

    def test_uncompressed_txt_file(self, tmp_path):
        from credtools.plot import _load_locus_pips

        locus_dir = tmp_path / "locus_txt"
        locus_dir.mkdir()
        df = pd.DataFrame(
            {
                "BP": [100, 200, 300],
                "CRED": [0, 1, 0],
                "PIP": [0.1, 0.8, 0.05],
            }
        )
        df.to_csv(locus_dir / "pips.txt", sep="\t", index=False)
        result = _load_locus_pips(locus_dir)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3

    def test_directory_not_found(self, tmp_path):
        from credtools.plot import _load_locus_pips

        with pytest.raises(FileNotFoundError, match="Locus directory not found"):
            _load_locus_pips(tmp_path / "nonexistent_dir")

    def test_not_a_directory(self, tmp_path):
        from credtools.plot import _load_locus_pips

        file_path = tmp_path / "not_a_dir.txt"
        file_path.write_text("dummy")
        with pytest.raises(NotADirectoryError, match="not a directory"):
            _load_locus_pips(file_path)

    def test_missing_required_columns(self, tmp_path):
        from credtools.plot import _load_locus_pips

        locus_dir = tmp_path / "bad_locus"
        locus_dir.mkdir()
        df = pd.DataFrame({"BP": [100], "EXTRA": [0.5]})
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        with pytest.raises(ValueError, match="missing required column"):
            _load_locus_pips(locus_dir)

    def test_nan_bp_rows_are_dropped(self, tmp_path):
        from credtools.plot import _load_locus_pips

        locus_dir = tmp_path / "locus_nan"
        locus_dir.mkdir()
        df = pd.DataFrame(
            {
                "BP": [100, np.nan, 300],
                "CRED": [1, 0, np.nan],
                "PIP": [0.5, 0.3, 0.1],
            }
        )
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        result = _load_locus_pips(locus_dir)
        assert len(result) == 2
        assert result["CRED"].dtype == int


# ---------------------------------------------------------------------------
# TestPlotLocusPvalues
# ---------------------------------------------------------------------------
class TestPlotLocusPvalues:
    """Tests for plot_locus_pvalues function."""

    def test_basic_zscore_data(self, tmp_path):
        rng = np.random.default_rng(99)
        n = 20
        df = pd.DataFrame(
            {
                "BP": np.arange(1000, 1000 + n * 100, 100),
                "z": rng.normal(0, 3, n),
            }
        )
        path = tmp_path / "expected_z.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_locus_pvalues(path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)

    def test_missing_bp_or_z_shows_error_text(self, tmp_path):
        df = pd.DataFrame({"other_col": [1, 2, 3]})
        path = tmp_path / "expected_z.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        fig, ax = plt.subplots()
        try:
            result = plot.plot_locus_pvalues(path, ax=ax)
            assert result is ax
            # Verify error text was placed on the axes
            texts = [t.get_text() for t in ax.texts]
            assert any("Required columns" in t for t in texts)
        finally:
            plt.close(fig)

    def test_with_credible_sets_file(self, tmp_path):
        rng = np.random.default_rng(42)
        n = 10
        bps = np.arange(5000, 5000 + n * 100, 100)
        z_df = pd.DataFrame({"BP": bps, "z": rng.normal(0, 3, n)})
        z_path = tmp_path / "expected_z.txt.gz"
        z_df.to_csv(z_path, sep="\t", index=False, compression="gzip")

        cs_df = pd.DataFrame(
            {
                "BP": [bps[2], bps[5]],
                "PIP": [0.9, 0.7],
            }
        )
        cs_path = tmp_path / "credible_sets.txt.gz"
        cs_df.to_csv(cs_path, sep="\t", index=False, compression="gzip")

        fig, ax = plt.subplots()
        try:
            result = plot.plot_locus_pvalues(z_path, credible_sets_file=cs_path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)

    def test_creates_new_figure_when_no_ax(self, tmp_path):
        df = pd.DataFrame(
            {
                "BP": [1000, 2000, 3000],
                "z": [2.5, -1.3, 4.0],
            }
        )
        path = tmp_path / "ez.txt.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        try:
            result = plot.plot_locus_pvalues(path)
            assert isinstance(result, plt.Axes)
        finally:
            plt.close("all")


# ---------------------------------------------------------------------------
# TestPlotSummaryQc
# ---------------------------------------------------------------------------
class TestPlotSummaryQc:
    """Tests for plot_summary_qc function."""

    def test_valid_qc_file_returns_figure(self, qc_summary_gz):
        fig = plot.plot_summary_qc(qc_summary_gz)
        try:
            assert isinstance(fig, plt.Figure)
            axes = fig.get_axes()
            assert len(axes) == 4
        finally:
            plt.close(fig)

    def test_with_output_file(self, qc_summary_gz, tmp_path):
        output_path = tmp_path / "summary_qc.png"
        fig = plot.plot_summary_qc(qc_summary_gz, output_file=output_path)
        try:
            assert output_path.exists()
            assert output_path.stat().st_size > 0
        finally:
            plt.close(fig)

    def test_without_output_file(self, qc_summary_gz):
        fig = plot.plot_summary_qc(qc_summary_gz)
        try:
            assert isinstance(fig, plt.Figure)
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# TestPlotLocusPlotEdgeCases
# ---------------------------------------------------------------------------
class TestPlotLocusPlotEdgeCases:
    """Tests for plot_locusplot edge cases."""

    def test_no_credible_variants(self, tmp_path):
        """All CRED=0 should still produce a plot without credible set markers."""
        locus_dir = tmp_path / "chr1_100_200"
        locus_dir.mkdir()
        rng = np.random.default_rng(42)
        n = 6
        df = pd.DataFrame(
            {
                "BP": np.arange(100, 100 + n * 10, 10),
                "CRED": [0] * n,
                "PIP": rng.uniform(0, 0.3, n).round(3),
                "UKB_P": rng.uniform(1e-8, 0.05, n),
                "UKB_R2": rng.uniform(0, 1, n).round(3),
            }
        )
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        fig = plot.plot_locusplot(locus_dir, figsize=(6, 4), dpi=72)
        try:
            assert isinstance(fig, plt.Figure)
        finally:
            plt.close(fig)

    def test_multiple_cohorts(self, tmp_path):
        """Multiple _P columns should produce one subplot per cohort."""
        locus_dir = tmp_path / "chr2_300_400"
        locus_dir.mkdir()
        rng = np.random.default_rng(7)
        n = 5
        df = pd.DataFrame(
            {
                "BP": np.arange(300, 300 + n * 10, 10),
                "CRED": [0, 1, 0, 0, 1],
                "PIP": rng.uniform(0, 1, n).round(3),
                "UKB_P": rng.uniform(1e-10, 0.05, n),
                "UKB_R2": rng.uniform(0, 1, n).round(3),
                "BBJ_P": rng.uniform(1e-10, 0.05, n),
                "BBJ_R2": rng.uniform(0, 1, n).round(3),
            }
        )
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        fig = plot.plot_locusplot(locus_dir, figsize=(6, 8), dpi=72)
        try:
            assert isinstance(fig, plt.Figure)
            axes = fig.get_axes()
            assert len(axes) == 2
        finally:
            plt.close(fig)

    def test_no_output_file_returns_figure(self, locus_dir_with_pips):
        """plot_locusplot with output_file=None should return figure without saving."""
        fig = plot.plot_locusplot(locus_dir_with_pips, figsize=(6, 4), dpi=72)
        try:
            assert isinstance(fig, plt.Figure)
        finally:
            plt.close(fig)

    def test_no_p_columns_raises(self, tmp_path):
        """Raise ValueError when pips file lacks _P columns."""
        locus_dir = tmp_path / "locus_nop"
        locus_dir.mkdir()
        df = pd.DataFrame(
            {
                "BP": [100, 200],
                "CRED": [0, 1],
                "PIP": [0.2, 0.8],
            }
        )
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        with pytest.raises(ValueError, match="No cohort-specific p-value columns"):
            plot.plot_locusplot(locus_dir)

    def test_missing_r2_column_raises(self, tmp_path):
        """Having _P column but missing corresponding _R2 should raise ValueError."""
        locus_dir = tmp_path / "locus_nor2"
        locus_dir.mkdir()
        df = pd.DataFrame(
            {
                "BP": [100, 200],
                "CRED": [0, 1],
                "PIP": [0.2, 0.8],
                "UKB_P": [1e-5, 1e-8],
            }
        )
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        with pytest.raises(ValueError, match="Missing R2 column"):
            plot.plot_locusplot(locus_dir)

    def test_empty_cohort_data(self, tmp_path):
        """Cohort with all NaN p-values should display 'No data available' text."""
        locus_dir = tmp_path / "locus_empty_cohort"
        locus_dir.mkdir()
        df = pd.DataFrame(
            {
                "BP": [100, 200, 300],
                "CRED": [0, 1, 0],
                "PIP": [0.1, 0.5, 0.2],
                "UKB_P": [float("nan"), float("nan"), float("nan")],
                "UKB_R2": [0.5, 0.6, 0.3],
            }
        )
        df.to_csv(locus_dir / "pips.txt.gz", sep="\t", index=False, compression="gzip")
        fig = plot.plot_locusplot(locus_dir, figsize=(6, 4), dpi=72)
        try:
            assert isinstance(fig, plt.Figure)
            ax = fig.get_axes()[0]
            texts = [t.get_text() for t in ax.texts]
            assert any("No data available" in t for t in texts)
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# TestPlotLocusPvaluesEdgeCases
# ---------------------------------------------------------------------------
class TestPlotLocusPvaluesEdgeCases:
    """Additional edge-case tests for plot_locus_pvalues."""

    def test_bad_credible_sets_file_logs_warning(self, tmp_path):
        """A malformed credible-sets file should not crash, just log a warning."""
        rng = np.random.default_rng(11)
        n = 5
        z_df = pd.DataFrame(
            {
                "BP": np.arange(1000, 1000 + n * 100, 100),
                "z": rng.normal(0, 2, n),
            }
        )
        z_path = tmp_path / "expected_z.txt.gz"
        z_df.to_csv(z_path, sep="\t", index=False, compression="gzip")

        # Write an invalid credible-sets file (missing PIP column)
        cs_df = pd.DataFrame({"BP": [1000], "WRONG_COL": [0.5]})
        cs_path = tmp_path / "cs_bad.txt.gz"
        cs_df.to_csv(cs_path, sep="\t", index=False, compression="gzip")

        fig, ax = plt.subplots()
        try:
            result = plot.plot_locus_pvalues(z_path, credible_sets_file=cs_path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)

    def test_credible_sets_no_matching_bps(self, tmp_path):
        """Credible set BPs that do not match z-data should not crash."""
        rng = np.random.default_rng(55)
        n = 5
        z_df = pd.DataFrame(
            {
                "BP": np.arange(1000, 1000 + n * 100, 100),
                "z": rng.normal(0, 2, n),
            }
        )
        z_path = tmp_path / "expected_z.txt.gz"
        z_df.to_csv(z_path, sep="\t", index=False, compression="gzip")

        cs_df = pd.DataFrame({"BP": [99999], "PIP": [0.95]})
        cs_path = tmp_path / "cs_nomatch.txt.gz"
        cs_df.to_csv(cs_path, sep="\t", index=False, compression="gzip")

        fig, ax = plt.subplots()
        try:
            result = plot.plot_locus_pvalues(z_path, credible_sets_file=cs_path, ax=ax)
            assert result is ax
        finally:
            plt.close(fig)
