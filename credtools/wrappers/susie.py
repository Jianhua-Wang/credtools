"""Wrapper of SuSiE fine-mapping method."""

import json
import logging
from typing import List

import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import Locus, intersect_sumstat_ld
from credtools.wrappers.susie_rss import susie_rss

logger = logging.getLogger("SuSiE")


def run_susie(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    max_iter: int = 100,
    estimate_residual_variance: bool = False,
    min_abs_corr: float = 0.1,
    convergence_tol: float = 1e-3,
) -> CredibleSet:
    """
    Run SuSiE (Sum of Single Effects) fine-mapping with shotgun stochastic search.

    SuSiE is a Bayesian variable selection method that fits a regression model
    with multiple single-effect components. Each component captures one causal
    signal, allowing the method to identify multiple independent associations
    while accounting for linkage disequilibrium.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix data.
        Must have matched summary statistics and LD matrix.
    max_causal : int, optional
        Maximum number of causal variants (L parameter), by default 1.
        This determines the number of single-effect components in the model.
        Higher values allow detection of more independent signals.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
        This determines the cumulative posterior probability mass
        included in each credible set.
    max_iter : int, optional
        Maximum number of iterations for the IBSS algorithm, by default 100.
        More iterations may improve convergence but increase runtime.
    estimate_residual_variance : bool, optional
        Whether to estimate residual variance from data, by default False.
        If False, residual variance is set to 1 (appropriate for z-scores).
    min_abs_corr : float, optional
        Minimum absolute correlation threshold for credible set purity, by default 0.1.
        Credible sets with pairwise correlations below this threshold may be filtered.
    convergence_tol : float, optional
        Convergence tolerance for the ELBO (Evidence Lower BOund), by default 1e-3.
        Algorithm stops when ELBO change falls below this threshold.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Posterior inclusion probabilities for all variants
        - Credible sets for each detected signal
        - Lead SNPs (highest PIP in each credible set)
        - Coverage probabilities and purity measures
        - Algorithm parameters and convergence information

    Warnings
    --------
    If the summary statistics and LD matrix are not matched, they will be
    automatically intersected and reordered with a warning message.

    Notes
    -----
    SuSiE implements the "Sum of Single Effects" model:

    y = Σ(l=1 to L) X * b_l + ε

    where:
    - y is the phenotype vector
    - X is the genotype matrix
    - b_l is the l-th single-effect vector (sparse, at most one non-zero element)
    - ε is the residual error
    - L is the maximum number of causal variants

    Key features of the algorithm:

    1. **Iterative Bayesian Stepwise Selection (IBSS)**: Iteratively updates
       each single-effect component while holding others fixed

    2. **Automatic multiplicity control**: Each component models exactly one
       causal variant, naturally controlling for multiple testing

    3. **Credible sets with coverage guarantees**: Provides sets of variants
       with specified coverage probability for each signal

    4. **LD-aware inference**: Properly accounts for correlation structure
       through the LD matrix

    The algorithm workflow:
    1. Initialize L single-effect regressions with uniform priors
    2. Iteratively update each effect using variational Bayes
    3. Monitor convergence using Evidence Lower Bound (ELBO)
    4. Construct credible sets based on posterior inclusion probabilities
    5. Filter credible sets based on purity criteria

    Advantages:
    - Principled uncertainty quantification through credible sets
    - Natural handling of multiple causal variants
    - Robust to LD structure and population stratification
    - Computationally efficient for genome-wide analysis

    Reference:
    Wang, G. et al. A simple new approach to variable selection in regression,
    with application to genetic fine mapping. J. R. Stat. Soc. B 82, 1273-1300 (2020).

    Examples
    --------
    >>> # Basic SuSiE analysis with single causal variant
    >>> credible_set = run_susie(locus)
    >>> print(f"Found {credible_set.n_cs} credible sets")
    >>> print(f"Top PIP: {credible_set.pips.max():.4f}")
    Found 1 credible sets
    Top PIP: 0.8234

    >>> # SuSiE with multiple causal variants and strict convergence
    >>> credible_set = run_susie(
    ...     locus,
    ...     max_causal=5,
    ...     coverage=0.99,
    ...     max_iter=500,
    ...     convergence_tol=1e-6
    ... )
    >>> print(f"Detected {credible_set.n_cs} independent signals")
    >>> print(f"Credible set sizes: {credible_set.cs_sizes}")
    Detected 3 independent signals
    Credible set sizes: [12, 8, 5]

    >>> # Access posterior inclusion probabilities
    >>> top_variants = credible_set.pips.nlargest(10)
    >>> print("Top 10 variants by PIP:")
    >>> print(top_variants)
    Top 10 variants by PIP:
    rs123456    0.8234
    rs789012    0.7456
    rs345678    0.6789
    ...

    >>> # Examine credible sets
    >>> for i, snps in enumerate(credible_set.snps):
    ...     print(f"Credible set {i+1}: {len(snps)} variants")
    ...     print(f"Lead SNP: {credible_set.lead_snps[i]}")
    Credible set 1: 12 variants
    Lead SNP: rs123456
    Credible set 2: 8 variants
    Lead SNP: rs789012
    """
    if not locus.is_matched:
        logger.warning(
            "The sumstat and LD are not matched, will match them in same order."
        )
        locus = intersect_sumstat_ld(locus)
    logger.info(f"Running SuSiE on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "max_iter": max_iter,
        "estimate_residual_variance": estimate_residual_variance,
        "min_abs_corr": min_abs_corr,
        "convergence_tol": convergence_tol,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    s = susie_rss(
        bhat=locus.sumstats[ColName.BETA].to_numpy(),
        shat=locus.sumstats[ColName.SE].to_numpy(),
        n=locus.sample_size,
        R=locus.ld.r,
        L=max_causal,
        coverage=coverage,
        max_iter=max_iter,
        estimate_residual_variance=estimate_residual_variance,
        min_abs_corr=min_abs_corr,
        tol=convergence_tol,
    )
    pip = s["pip"]
    if s["sets"]["cs"]:
        cs_idx = list(s["sets"]["cs"].values())
        n_cs = len(cs_idx)
        cs_sizes = [len(idx) for idx in cs_idx]
        cred_snps = [locus.sumstats[ColName.SNPID].iloc[idx].tolist() for idx in cs_idx]
    else:
        n_cs = 0
        cs_sizes = []
        cred_snps = []
    pips = pd.Series(data=pip, index=locus.sumstats[ColName.SNPID].tolist())
    lead_snps = [str(pips[pips.index.isin(cred_snps[i])].idxmax()) for i in range(n_cs)]
    logger.info(f"Finished SuSiE on {locus}")
    logger.info(f"N of credible set: {n_cs}")
    logger.info(f"Credible set size: {cs_sizes}")
    return CredibleSet(
        tool=Method.SUSIE,
        n_cs=n_cs,
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cred_snps,
        cs_sizes=cs_sizes,
        pips=pips,
        parameters=parameters,
    )
