import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path

import pytest

from credtools import plot

QC_SUMMARY_PATH = Path("exampledata/testout/qc/qc.txt.gz")


@pytest.mark.parametrize(
    "plot_func, kwargs",
    (
        (plot.plot_lambda_s_boxplot, {}),
        (plot.plot_maf_corr_barplot, {}),
        (plot.plot_outliers_barplot, {"outlier_type": "lambda_s"}),
    ),
)
def test_qc_plots_accept_path_input(plot_func, kwargs):
    fig, ax = plt.subplots()
    try:
        result = plot_func(QC_SUMMARY_PATH, ax=ax, **kwargs)
    finally:
        plt.close(fig)
    assert result is ax


LOCUS_DIR = Path("exampledata/testout/susie/chr1_49782265_50282265")


def test_plot_locusplot(tmp_path):
    output_path = tmp_path / "locus_plot.png"
    fig = plot.plot_locusplot(
        LOCUS_DIR, output_file=output_path, figsize=(6, 4), dpi=72
    )
    plt.close(fig)
    assert output_path.exists()
    assert output_path.stat().st_size > 0
