import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import pytest

from credtools import plot


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
