"""Wrapper for CARMA fine-mapping."""

import io
import json
import logging
from contextlib import redirect_stdout
from typing import Optional

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import Locus, intersect_sumstat_ld
from credtools.utils import check_r_package, io_in_tempdir

logger = logging.getLogger("CARMA")


@io_in_tempdir("./tmp/CARMA")
def run_carma(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    effect_size_prior: str = "Spike-slab",
    input_alpha: float = 0.0,
    y_var: float = 1.0,
    bf_threshold: float = 10.0,
    outlier_bf_threshold: float = 1 / 3.2,
    outlier_switch: bool = True,
    max_model_dim: int = 200000,
    all_inner_iter: int = 10,
    all_iter: int = 3,
    tau: float = 0.04,
    epsilon_threshold: float = 1e-5,
    printing_log: bool = False,
    em_dist: str = "logistic",
    temp_dir: Optional[str] = None,
) -> CredibleSet:
    """
    Run CARMA fine-mapping using R through rpy2.

    CARMA (Causal And Robust Multi-ancestry) performs Bayesian fine-mapping
    using a spike-and-slab prior with automatic model selection. The method
    can handle multiple causal variants and provides robust inference in the
    presence of outliers and model misspecification.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix.
        Must have matched summary statistics and LD data.
    max_causal : int, optional
        Maximum number of causal variants assumed per locus, by default 1.
        Higher values increase computational complexity but allow for
        more realistic multi-variant models.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
        This determines the probability mass included in each credible set.
    effect_size_prior : str, optional
        Prior distribution for effect sizes, by default "Spike-slab".
        Options: "Cauchy" or "Spike-slab". Spike-slab is generally recommended
        for fine-mapping as it provides better sparsity properties.
    input_alpha : float, optional
        Elastic net mixing parameter (0 ≤ input_alpha ≤ 1), by default 0.0.
        Controls the balance between ridge (α=0) and lasso (α=1) regularization.
    y_var : float, optional
        Variance of the summary statistics, by default 1.0.
        Typically set to 1.0 for standardized effect sizes.
    bf_threshold : float, optional
        Bayes factor threshold for credible models, by default 10.0.
        Models with BF > threshold are included in the credible set.
    outlier_bf_threshold : float, optional
        Bayes factor threshold for outlier detection, by default 1/3.2.
        Variants with BF < threshold may be flagged as outliers.
    outlier_switch : bool, optional
        Whether to perform outlier detection, by default True.
        Enables robust inference against outlying variants.
    max_model_dim : int, optional
        Maximum number of top candidate models to consider, by default 200000.
        Larger values provide more comprehensive search but increase runtime.
    all_inner_iter : int, optional
        Maximum iterations for Shotgun algorithm within EM, by default 10.
        Controls the inner optimization loop convergence.
    all_iter : int, optional
        Maximum iterations for EM algorithm, by default 3.
        Controls the outer EM algorithm convergence.
    tau : float, optional
        Prior precision parameter of effect size, by default 0.04.
        Smaller values correspond to more diffuse priors on effect sizes.
    epsilon_threshold : float, optional
        Convergence threshold for Bayes factors, by default 1e-5.
        Algorithm stops when relative changes fall below this threshold.
    printing_log : bool, optional
        Whether to print CARMA running log, by default False.
        Enables verbose output for debugging purposes.
    em_dist : str, optional
        Distribution for modeling prior probability, by default "logistic".
        Currently only "logistic" is supported.
    temp_dir : Optional[str], optional
        Temporary directory for intermediate files, by default None.
        Automatically provided by the @io_in_tempdir decorator.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Posterior inclusion probabilities for all variants
        - Credible set assignments for different configurations
        - Lead SNPs with highest PIPs in each credible set
        - Model selection and convergence information

    Warnings
    --------
    If the summary statistics and LD matrix are not matched, they will be
    automatically intersected and reordered with a warning message.

    Notes
    -----
    CARMA implements a Bayesian fine-mapping framework with several key features:

    1. **Spike-and-slab priors**: Provides automatic variable selection with
       interpretable posterior inclusion probabilities

    2. **Multi-ancestry support**: Can handle multiple ancestries with
       appropriate LD matrices (though this wrapper focuses on single ancestry)

    3. **Outlier detection**: Robust inference against outlying variants
       that may violate model assumptions

    4. **EM algorithm**: Efficient optimization using expectation-maximization
       with shotgun stochastic search

    The algorithm workflow:
    1. Convert summary statistics to Z-scores
    2. Set up R environment and import CARMA package
    3. Configure priors and algorithm parameters
    4. Run EM algorithm with shotgun search
    5. Extract posterior inclusion probabilities and credible sets
    6. Identify lead variants within each credible set

    Model assumptions:
    - Summary statistics follow multivariate normal distribution
    - LD matrix accurately reflects population structure
    - Effect sizes follow spike-and-slab or Cauchy priors
    - At most max_causal variants are truly causal

    Reference:
    Jin, J. et al. A powerful approach to estimating a composite null model in
    large-scale genomics data. arXiv preprint (2021).

    Examples
    --------
    >>> # Basic CARMA analysis with default parameters
    >>> credible_set = run_carma(locus)
    >>> print(f"Found {credible_set.n_cs} credible sets")
    >>> print(f"Lead SNPs: {credible_set.lead_snps}")
    Found 1 credible sets
    Lead SNPs: ['rs123456']

    >>> # CARMA with multiple causal variants and Cauchy prior
    >>> credible_set = run_carma(
    ...     locus,
    ...     max_causal=3,
    ...     effect_size_prior="Cauchy",
    ...     coverage=0.99,
    ...     outlier_switch=True
    ... )
    >>> print(f"Credible set sizes: {credible_set.cs_sizes}")
    >>> print(f"Top PIP: {credible_set.pips.max():.4f}")
    Credible set sizes: [4, 7, 12]
    Top PIP: 0.8542

    >>> # Access detailed results
    >>> pips_df = credible_set.pips.sort_values(ascending=False)
    >>> print("Top 5 variants by PIP:")
    >>> print(pips_df.head())
    Top 5 variants by PIP:
    rs123456    0.8542
    rs789012    0.6431
    rs345678    0.4329
    rs456789    0.2108
    rs567890    0.1054
    """
    # if not check_r_package("CARMA"):
    #     raise RuntimeError("CARMA is not installed or R version is not supported.")
    if not locus.is_matched:
        logger.warning(
            "The sumstat and LD are not matched, will match them in same order."
        )
        locus = intersect_sumstat_ld(locus)
    logger.info(f"Running CARMA on {locus}")

    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "effect_size_prior": effect_size_prior,
        "input_alpha": input_alpha,
        "y_var": y_var,
        "bf_threshold": bf_threshold,
        "outlier_bf_threshold": outlier_bf_threshold,
        "outlier_switch": outlier_switch,
        "max_model_dim": max_model_dim,
        "all_inner_iter": all_inner_iter,
        "all_iter": all_iter,
        "tau": tau,
        "epsilon_threshold": epsilon_threshold,
        "em_dist": em_dist,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    sumstats = locus.sumstats.copy()
    ld = locus.ld.r.copy()
    sumstats[ColName.Z] = sumstats[ColName.BETA] / sumstats[ColName.SE]

    # Import required R packages
    carma = importr("CARMA")

    # Create R lists for input
    z_list = ro.ListVector({"1": ro.FloatVector(sumstats[ColName.Z].values)})
    ld_matrix = ro.r.matrix(ro.FloatVector(ld.flatten()), nrow=ld.shape[0], ncol=ld.shape[1])  # type: ignore
    ld_list = ro.ListVector({"1": ld_matrix})
    lambda_list = ro.ListVector({"1": ro.FloatVector([1.0])})

    # Run CARMA with all parameters
    f = io.StringIO()
    with redirect_stdout(f):
        carma_results = carma.CARMA(
            z_list,
            ld_list,
            lambda_list=lambda_list,
            effect_size_prior=effect_size_prior,
            input_alpha=input_alpha,
            y_var=y_var,
            rho_index=coverage,
            BF_index=bf_threshold,
            outlier_BF_index=outlier_bf_threshold,
            outlier_switch=outlier_switch,
            num_causal=max_causal,
            Max_Model_Dim=max_model_dim,
            all_inner_iter=all_inner_iter,
            all_iter=all_iter,
            tau=tau,
            epsilon_threshold=epsilon_threshold,
            printing_log=printing_log,
            EM_dist=em_dist,
            output_labels=temp_dir,
        )
    logger.debug(f.getvalue())

    # Extract PIPs
    pips = np.array(carma_results[0].rx2("PIPs"))

    # Extract credible sets
    cs = np.zeros(len(sumstats))
    credible_sets = carma_results[0].rx2("Credible set")[1]
    if credible_sets:
        for i, cs_indices in enumerate(credible_sets, 1):
            cs[np.array(cs_indices, dtype=int) - 1] = i  # R uses 1-based indexing

    # Add results to summary statistics
    result_df = sumstats.copy()
    result_df["PIP"] = pips
    result_df["CS"] = cs.astype(int)

    pips = pd.Series(
        data=result_df["PIP"].to_numpy(), index=result_df["SNPID"].to_numpy()
    )

    cs_snps = []
    lead_snps = []
    for cs_i, sub_df in result_df.groupby("CS"):
        if cs_i == 0:
            continue
        cs_snps.append(sub_df["SNPID"].values.tolist())
        lead_snps.append(pips[pips.index.isin(sub_df["SNPID"].values)].idxmax())
    return CredibleSet(
        tool=Method.CARMA,
        n_cs=len(cs_snps),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=[len(i) for i in cs_snps],
        pips=pips,
        parameters=parameters,
    )
