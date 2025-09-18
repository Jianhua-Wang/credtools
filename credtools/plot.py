"""Plotting functions for credtools."""

import gzip
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Circle
from scipy import stats

try:
    from upsetplot import UpSet

    UPSETPLOT_AVAILABLE = True
except ImportError:
    UPSETPLOT_AVAILABLE = False

# Set up matplotlib style
plt.style.use("default")
sns.set_palette("husl")

# Population color mapping for consistency
POPULATION_COLORS = {
    "AFR": "#FF6B6B",  # Red
    "EAS": "#4ECDC4",  # Teal
    "EUR": "#45B7D1",  # Blue
    "SAS": "#FFA07A",  # Orange
    "AMR": "#98D8C8",  # Mint
    "HIS": "#DDA0DD",  # Plum
}


def read_compressed_file(file_path: Union[str, Path]) -> pd.DataFrame:
    """
    Read a compressed or uncompressed file.

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to the file to read.

    Returns
    -------
    pd.DataFrame
        Loaded dataframe.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    if file_path.suffix == ".gz":
        return pd.read_csv(file_path, sep="\t", compression="gzip")
    else:
        return pd.read_csv(file_path, sep="\t")


def get_population_color(cohort: str) -> str:
    """
    Get color for population based on cohort name.

    Parameters
    ----------
    cohort : str
        Cohort name in format "POPULATION_COHORT".

    Returns
    -------
    str
        Hex color code.
    """
    population = cohort.split("_")[0]
    return POPULATION_COLORS.get(population, "#7F7F7F")  # Default gray


def plot_lambda_s_boxplot(
    qc_data: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> plt.Axes:
    """
    Create boxplot of lambda-s values by cohort.

    Parameters
    ----------
    qc_data : pd.DataFrame
        QC summary data with 'lambda_s' and cohort columns.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create cohort identifier column
    qc_data = qc_data.copy()
    qc_data["cohort_id"] = qc_data["popu"] + "_" + qc_data["cohort"]

    # Create colors for each population
    colors = [get_population_color(cohort) for cohort in qc_data["cohort_id"]]

    # Create boxplot
    box_plot = ax.boxplot(
        [qc_data[qc_data["cohort_id"] == cohort]["lambda_s"].values for cohort in qc_data["cohort_id"].unique()],
        labels=qc_data["cohort_id"].unique(),
        patch_artist=True,
    )

    # Color the boxes
    for patch, color in zip(
        box_plot["boxes"],
        [get_population_color(c) for c in qc_data["cohort_id"].unique()],
    ):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel("lambda_s")
    ax.set_xlabel("cohort")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(True, alpha=0.3)

    # Add legend for populations
    populations = sorted(set(c.split("_")[0] for c in qc_data["cohort_id"].unique()))
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=POPULATION_COLORS.get(pop, "#7F7F7F"),
            alpha=0.7,
            label=pop,
        )
        for pop in populations
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    return ax


def plot_maf_corr_barplot(
    qc_data: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> plt.Axes:
    """
    Create barplot of MAF correlations by cohort.

    Parameters
    ----------
    qc_data : pd.DataFrame
        QC summary data with 'maf_corr' and cohort columns.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create cohort identifier column
    qc_data = qc_data.copy()
    qc_data["cohort_id"] = qc_data["popu"] + "_" + qc_data["cohort"]

    # Remove NaN values
    plot_data = qc_data.dropna(subset=["maf_corr"])

    if plot_data.empty:
        ax.text(
            0.5,
            0.5,
            "No MAF correlation data available",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Create colors for each cohort
    colors = [get_population_color(cohort) for cohort in plot_data["cohort_id"]]

    # Create barplot
    bars = ax.bar(range(len(plot_data)), plot_data["maf_corr"], color=colors, alpha=0.7)

    ax.set_ylabel("MAF correlation")
    ax.set_xlabel("cohort")
    ax.set_xticks(range(len(plot_data)))
    ax.set_xticklabels(plot_data["cohort_id"], rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_ylim(0, 1)

    # Add legend for populations
    populations = sorted(set(c.split("_")[0] for c in plot_data["cohort_id"].unique()))
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=POPULATION_COLORS.get(pop, "#7F7F7F"),
            alpha=0.7,
            label=pop,
        )
        for pop in populations
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    return ax


def plot_outliers_barplot(
    qc_data: pd.DataFrame,
    outlier_type: str = "lambda_s",
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> plt.Axes:
    """
    Create barplot of outlier counts by cohort.

    Parameters
    ----------
    qc_data : pd.DataFrame
        QC summary data with outlier count columns.
    outlier_type : str
        Type of outlier to plot ('lambda_s' or 'dentist_s').
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create cohort identifier column
    qc_data = qc_data.copy()
    qc_data["cohort_id"] = qc_data["popu"] + "_" + qc_data["cohort"]

    outlier_col = f"n_{outlier_type}_outlier"
    if outlier_col not in qc_data.columns:
        ax.text(
            0.5,
            0.5,
            f"No {outlier_type} outlier data available",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Create colors for each cohort
    colors = [get_population_color(cohort) for cohort in qc_data["cohort_id"]]

    # Create barplot
    bars = ax.bar(range(len(qc_data)), qc_data[outlier_col], color=colors, alpha=0.7)

    ax.set_ylabel(f"n_{outlier_type}_outlier")
    ax.set_xlabel("cohort")
    ax.set_xticks(range(len(qc_data)))
    ax.set_xticklabels(qc_data["cohort_id"], rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    # Add legend for populations
    populations = sorted(set(c.split("_")[0] for c in qc_data["cohort_id"].unique()))
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=POPULATION_COLORS.get(pop, "#7F7F7F"),
            alpha=0.7,
            label=pop,
        )
        for pop in populations
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    return ax


def plot_summary_qc(
    qc_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (16, 12),
    dpi: int = 300,
) -> plt.Figure:
    """
    Create 2x2 summary QC plot from aggregated QC data.

    Parameters
    ----------
    qc_file : Union[str, Path]
        Path to QC summary file (qc.txt.gz).
    output_file : Optional[Union[str, Path]]
        Output file path. If None, displays plot.
    figsize : Tuple[float, float]
        Figure size.
    dpi : int
        DPI for output file.

    Returns
    -------
    plt.Figure
        Matplotlib figure object.
    """
    # Read QC data
    qc_data = read_compressed_file(qc_file)

    # Create 2x2 subplot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: Lambda-s boxplot
    plot_lambda_s_boxplot(qc_data, ax=ax1)
    ax1.set_title("Lambda-s Distribution by Cohort")

    # Plot 2: MAF correlation barplot
    plot_maf_corr_barplot(qc_data, ax=ax2)
    ax2.set_title("MAF Correlation by Cohort")

    # Plot 3: Lambda-s outliers
    plot_outliers_barplot(qc_data, outlier_type="lambda_s", ax=ax3)
    ax3.set_title("Lambda-s Outliers by Cohort")

    # Plot 4: Dentist-s outliers
    plot_outliers_barplot(qc_data, outlier_type="dentist_s", ax=ax4)
    ax4.set_title("Dentist-s Outliers by Cohort")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")

    return fig


def plot_locus_pvalues(
    expected_z_file: Union[str, Path],
    credible_sets_file: Optional[Union[str, Path]] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (12, 6),
) -> plt.Axes:
    """
    Create locus plot with p-values and credible set annotations.

    Parameters
    ----------
    expected_z_file : Union[str, Path]
        Path to expected_z.txt.gz file.
    credible_sets_file : Optional[Union[str, Path]]
        Path to credible sets file for annotations.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Read expected z data
    z_data = read_compressed_file(expected_z_file)

    if "BP" not in z_data.columns or "z" not in z_data.columns:
        ax.text(
            0.5,
            0.5,
            "Required columns (BP, z) not found in data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Calculate p-values from z-scores
    p_values = 2 * (1 - stats.norm.cdf(np.abs(z_data["z"])))
    neg_log_p = -np.log10(p_values)

    # Plot points
    ax.scatter(z_data["BP"], neg_log_p, s=20, alpha=0.7, color="#1f77b4")

    # Add significance lines
    ax.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.7, label="5e-8")
    ax.axhline(y=-np.log10(1e-5), color="orange", linestyle="--", alpha=0.7, label="1e-5")

    # Add credible set annotations if available
    if credible_sets_file and Path(credible_sets_file).exists():
        try:
            cs_data = read_compressed_file(credible_sets_file)
            if "BP" in cs_data.columns and "PIP" in cs_data.columns:
                # Highlight credible set variants
                cs_variants = cs_data[cs_data["PIP"] > 0.01]  # Filter for meaningful PIP
                for _, variant in cs_variants.iterrows():
                    bp = variant["BP"]
                    # Find corresponding p-value
                    z_match = z_data[z_data["BP"] == bp]
                    if not z_match.empty:
                        p_val = 2 * (1 - stats.norm.cdf(np.abs(z_match["z"].iloc[0])))
                        y_pos = -np.log10(p_val)
                        circle = Circle(
                            (bp, y_pos),
                            radius=(ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.005,
                            color="red",
                            fill=False,
                            linewidth=2,
                        )
                        ax.add_patch(circle)
        except Exception as e:
            logging.warning(f"Could not add credible set annotations: {e}")

    ax.set_xlabel("Position (BP)")
    ax.set_ylabel("-log10(P-value)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    return ax


def plot_zscore_qq(
    expected_z_file: Union[str, Path],
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (6, 6),
) -> plt.Axes:
    """
    Create QQ plot of observed vs expected z-scores.

    Parameters
    ----------
    expected_z_file : Union[str, Path]
        Path to expected_z.txt.gz file.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Read expected z data
    z_data = read_compressed_file(expected_z_file)

    if "z" not in z_data.columns or "z_expected" not in z_data.columns:
        ax.text(
            0.5,
            0.5,
            "Required columns (z, z_expected) not found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Create QQ plot
    ax.scatter(z_data["z_expected"], z_data["z"], s=20, alpha=0.6, color="#1f77b4")

    # Add diagonal line
    min_val = min(z_data["z_expected"].min(), z_data["z"].min())
    max_val = max(z_data["z_expected"].max(), z_data["z"].max())
    ax.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8, linewidth=2)

    # Add lambda value if available
    if "lambda_s" in z_data.columns:
        lambda_val = z_data["lambda_s"].iloc[0]
        ax.text(
            0.05,
            0.95,
            f"lambda={lambda_val:.3f}",
            transform=ax.transAxes,
            fontsize=12,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    ax.set_xlabel("Expected z-score")
    ax.set_ylabel("Observed z-score")
    ax.grid(True, alpha=0.3)

    # Make axes equal
    ax.set_aspect("equal", adjustable="box")

    return ax


def plot_ld_decay(
    ld_decay_file: Union[str, Path],
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> plt.Axes:
    """
    Create LD decay line plot.

    Parameters
    ----------
    ld_decay_file : Union[str, Path]
        Path to ld_decay.txt.gz file.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Read LD decay data
    decay_data = read_compressed_file(ld_decay_file)

    required_cols = ["distance_kb", "r2_avg", "cohort"]
    if not all(col in decay_data.columns for col in required_cols):
        ax.text(
            0.5,
            0.5,
            f"Required columns {required_cols} not found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Plot each cohort
    cohorts = decay_data["cohort"].unique()
    for cohort in cohorts:
        cohort_data = decay_data[decay_data["cohort"] == cohort]
        color = get_population_color(cohort)
        ax.plot(
            cohort_data["distance_kb"],
            cohort_data["r2_avg"],
            color=color,
            label=cohort,
            linewidth=2,
            alpha=0.8,
        )

    ax.set_xlabel("Distance (kb)")
    ax.set_ylabel("r²")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_ylim(0, None)

    return ax


def plot_ld_4th_moment(
    ld_4th_file: Union[str, Path],
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> plt.Axes:
    """
    Create LD 4th moment boxplot.

    Parameters
    ----------
    ld_4th_file : Union[str, Path]
        Path to ld_4th_moment.txt.gz file.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Read LD 4th moment data
    ld_4th_data = read_compressed_file(ld_4th_file)

    # The file should have cohort columns
    cohort_cols = [col for col in ld_4th_data.columns if "_" in col]
    if not cohort_cols:
        ax.text(
            0.5,
            0.5,
            "No cohort columns found in LD 4th moment data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Prepare data for boxplot
    data_for_plot = []
    labels = []
    colors = []

    for col in cohort_cols:
        data_for_plot.append(ld_4th_data[col].dropna().values)
        labels.append(col)
        colors.append(get_population_color(col))

    # Create boxplot
    box_plot = ax.boxplot(data_for_plot, labels=labels, patch_artist=True)

    # Color the boxes
    for patch, color in zip(box_plot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel("LD 4th moment")
    ax.set_xlabel("Cohort")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(True, alpha=0.3)

    # Add legend for populations
    populations = sorted(set(c.split("_")[0] for c in cohort_cols))
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=POPULATION_COLORS.get(pop, "#7F7F7F"),
            alpha=0.7,
            label=pop,
        )
        for pop in populations
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    return ax


def plot_snp_missingness_upset(
    snp_missingness_file: Union[str, Path],
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (12, 8),
) -> plt.Axes:
    """
    Create SNP missingness upset plot.

    Parameters
    ----------
    snp_missingness_file : Union[str, Path]
        Path to snp_missingness.txt.gz file.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on. If None, creates new figure.
    figsize : Tuple[float, float]
        Figure size if creating new figure.

    Returns
    -------
    plt.Axes
        Matplotlib axes object.
    """
    if not UPSETPLOT_AVAILABLE:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5,
            0.5,
            "upsetplot package not available\nInstall with: pip install upsetplot",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Read SNP missingness data
    miss_data = read_compressed_file(snp_missingness_file)

    # Expected format: SNPID column + cohort columns with 0/1 values
    if "SNPID" not in miss_data.columns:
        ax.text(
            0.5,
            0.5,
            "SNPID column not found in missingness data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Get cohort columns (exclude SNPID)
    cohort_cols = [col for col in miss_data.columns if col != "SNPID"]
    if not cohort_cols:
        ax.text(
            0.5,
            0.5,
            "No cohort columns found in missingness data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    # Convert to boolean for upset plot (1 = present, 0 = missing)
    upset_data = miss_data[cohort_cols].astype(bool)
    upset_data.index = miss_data["SNPID"]

    try:
        # Create upset plot
        upset = UpSet(upset_data, subset_size="count", show_counts=True)
        upset.plot(fig=ax.figure if hasattr(ax, "figure") else plt.gcf())
    except Exception as e:
        ax.text(
            0.5,
            0.5,
            f"Error creating upset plot: {str(e)}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )

    return ax


def plot_locus_qc(
    locus_dir: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (16, 12),
    dpi: int = 300,
    include_upset: bool = True,
) -> plt.Figure:
    """
    Create 2x2 (or 2x3 with upset) locus-specific QC plot.

    Parameters
    ----------
    locus_dir : Union[str, Path]
        Directory containing locus QC files.
    output_file : Optional[Union[str, Path]]
        Output file path. If None, displays plot.
    figsize : Tuple[float, float]
        Figure size.
    dpi : int
        DPI for output file.
    include_upset : bool
        Whether to include upset plot (requires 2x3 layout).

    Returns
    -------
    plt.Figure
        Matplotlib figure object.
    """
    locus_dir = Path(locus_dir)

    # File paths
    expected_z_file = locus_dir / "expected_z.txt.gz"
    ld_decay_file = locus_dir / "ld_decay.txt.gz"
    ld_4th_file = locus_dir / "ld_4th_moment.txt.gz"
    snp_miss_file = locus_dir / "snp_missingness.txt.gz"

    # Create subplot layout
    if include_upset and UPSETPLOT_AVAILABLE and snp_miss_file.exists():
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(figsize[0], figsize[1] * 1.5))
        ax6.axis("off")  # Hide unused subplot
    else:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        ax5 = None

    # Plot 1: Locus p-values
    if expected_z_file.exists():
        plot_locus_pvalues(expected_z_file, ax=ax1)
        ax1.set_title("Locus P-values")
    else:
        ax1.text(
            0.5,
            0.5,
            "expected_z.txt.gz not found",
            ha="center",
            va="center",
            transform=ax1.transAxes,
        )

    # Plot 2: Z-score QQ plot
    if expected_z_file.exists():
        plot_zscore_qq(expected_z_file, ax=ax2)
        ax2.set_title("Observed vs Expected Z-scores")
    else:
        ax2.text(
            0.5,
            0.5,
            "expected_z.txt.gz not found",
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )

    # Plot 3: LD decay
    if ld_decay_file.exists():
        plot_ld_decay(ld_decay_file, ax=ax3)
        ax3.set_title("LD Decay")
    else:
        ax3.text(
            0.5,
            0.5,
            "ld_decay.txt.gz not found",
            ha="center",
            va="center",
            transform=ax3.transAxes,
        )

    # Plot 4: LD 4th moment
    if ld_4th_file.exists():
        plot_ld_4th_moment(ld_4th_file, ax=ax4)
        ax4.set_title("LD 4th Moment")
    else:
        ax4.text(
            0.5,
            0.5,
            "ld_4th_moment.txt.gz not found",
            ha="center",
            va="center",
            transform=ax4.transAxes,
        )

    # Plot 5: SNP missingness upset plot (if space and data available)
    if ax5 and snp_miss_file.exists():
        plot_snp_missingness_upset(snp_miss_file, ax=ax5)
        ax5.set_title("SNP Missingness Overlap")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")

    return fig
