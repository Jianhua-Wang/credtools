"""
SuSiE (Sum of Single Effects) model implementation for fine-mapping.

This module provides comprehensive functionality for fitting the SuSiE model
to summary statistics and full data, including credible set construction,
posterior inference, and model diagnostics.
"""

import logging
import math
import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy import optimize, sparse, stats

logger = logging.getLogger("SuSiE")


def susie_get_cs(
    res: Dict[str, np.ndarray],
    X: Optional[np.ndarray] = None,
    Xcorr: Optional[np.ndarray] = None,
    coverage: float = 0.95,
    min_abs_corr: float = 0.5,
    dedup: bool = True,
    squared: bool = False,
    check_symmetric: bool = True,
    n_purity: int = 100,
) -> Dict[str, Any]:
    """
    Extract credible sets from SuSiE results with quality assessment.

    This function processes the posterior inclusion probabilities from a SuSiE
    analysis to construct credible sets with specified coverage probability.
    It optionally filters credible sets based on correlation structure to
    ensure high-quality fine-mapping results.

    Parameters
    ----------
    res : Dict[str, np.ndarray]
        Dictionary containing SuSiE results. Must include:
        - 'alpha': L×P matrix of posterior inclusion probabilities
        - 'V': (optional) L-vector of prior variances for filtering
        - 'null_index': (optional) index of null component
    X : np.ndarray, optional
        n×p matrix of standardized genotype values.
        Either X or Xcorr must be provided for purity calculation.
    Xcorr : np.ndarray, optional
        p×p correlation matrix between variables.
        Either X or Xcorr must be provided for purity calculation.
    coverage : float, optional
        Target coverage probability for credible sets, by default 0.95.
        Each credible set will contain variants with cumulative PIP ≥ coverage.
    min_abs_corr : float, optional
        Minimum absolute correlation threshold for credible set purity, by default 0.5.
        Credible sets with maximum pairwise correlation below this threshold
        are filtered out as potentially unreliable.
    dedup : bool, optional
        Whether to remove duplicate credible sets, by default True.
        Duplicate sets can arise when multiple components identify the same signal.
    squared : bool, optional
        Whether to use squared correlations for purity assessment, by default False.
        If True, min_abs_corr is interpreted as minimum squared correlation.
    check_symmetric : bool, optional
        Whether to verify and enforce symmetry of Xcorr matrix, by default True.
        Non-symmetric correlation matrices are symmetrized by averaging.
    n_purity : int, optional
        Maximum number of variants per credible set for purity calculation, by default 100.
        Large credible sets are randomly subsampled for computational efficiency.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing credible set results:
        - 'cs': Dict mapping credible set names to variant indices
        - 'coverage': Array of actual coverage for each credible set
        - 'purity': Dict with correlation statistics (min, mean, median)
        - 'cs_index': Indices of SuSiE components corresponding to credible sets
        - 'requested_coverage': The input coverage parameter

        Returns None values if no credible sets pass quality filters.

    Raises
    ------
    ValueError
        If both X and Xcorr are provided (only one should be specified).

    Notes
    -----
    Credible set construction follows these steps:

    1. **Component filtering**: Remove components with very small prior variance
       (if 'V' is available) and null components

    2. **Coverage-based inclusion**: For each component, include variants in
       order of decreasing PIP until cumulative probability ≥ coverage

    3. **Deduplication**: Remove identical credible sets that may arise from
       multiple components capturing the same signal

    4. **Purity assessment**: Calculate correlation statistics within each
       credible set using provided genotype data or correlation matrix

    5. **Quality filtering**: Retain only credible sets meeting the minimum
       correlation threshold

    Purity metrics provide insight into credible set quality:
    - **min_abs_corr**: Minimum pairwise correlation (key quality metric)
    - **mean_abs_corr**: Average pairwise correlation
    - **median_abs_corr**: Median pairwise correlation

    High-quality credible sets should have:
    - min_abs_corr ≥ 0.5 (variants in strong LD)
    - Compact size (< 100 variants typically)
    - High coverage (> 0.95 for reliable inference)

    Examples
    --------
    >>> # Basic credible set extraction
    >>> cs_results = susie_get_cs(susie_fit)
    >>> print(f"Found {len(cs_results['cs'])} credible sets")
    >>> for name, variants in cs_results['cs'].items():
    ...     print(f"{name}: {len(variants)} variants")
    Found 2 credible sets
    L0: 12 variants
    L1: 8 variants

    >>> # High-stringency credible sets with strict purity requirements
    >>> cs_results = susie_get_cs(
    ...     susie_fit,
    ...     Xcorr=ld_matrix,
    ...     coverage=0.99,
    ...     min_abs_corr=0.8
    ... )
    >>> if cs_results['cs'] is not None:
    ...     purity = cs_results['purity']
    ...     print(f"Min correlations: {purity['min_abs_corr']}")
    ...     print(f"Mean correlations: {purity['mean_abs_corr']}")
    Min correlations: [0.85, 0.92]
    Mean correlations: [0.91, 0.95]

    >>> # Access specific credible set details
    >>> for cs_name, variants in cs_results['cs'].items():
    ...     idx = int(cs_name[1:])  # Extract component index
    ...     coverage = cs_results['coverage'][idx]
    ...     purity = cs_results['purity']['min_abs_corr'][idx]
    ...     print(f"Set {cs_name}: {len(variants)} variants, "
    ...           f"coverage={coverage:.3f}, purity={purity:.3f}")
    Set L0: 12 variants, coverage=0.982, purity=0.851
    Set L1: 8 variants, coverage=0.991, purity=0.923
    """
    if X is not None and Xcorr is not None:
        raise ValueError("Only one of X or Xcorr should be specified")

    if check_symmetric and Xcorr is not None:
        if not np.allclose(Xcorr, Xcorr.T):
            logger.warning("Xcorr is not symmetric; forcing Xcorr to be symmetric")
            Xcorr = (Xcorr + Xcorr.T) / 2

    null_index = 0
    include_idx = np.ones(res["alpha"].shape[0], dtype=bool)
    if "null_index" in res:
        null_index = res["null_index"]
    if "V" in res and isinstance(res["V"], np.ndarray):
        include_idx = res["V"] > 1e-9

    # L x P binary matrix
    status = in_CS(res["alpha"], coverage)

    # L-list of CS positions
    cs = [np.where(status[i] != 0)[0] for i in range(status.shape[0])]
    claimed_coverage = np.array([res["alpha"][i, cs[i]].sum() for i in range(len(cs))])
    include_idx = include_idx & (np.array([len(c) for c in cs]) > 0)

    if dedup:
        cs_tuples = [tuple(arr) for arr in cs]

        seen = set()
        not_duplicated = []
        for item in cs_tuples:
            if item not in seen:
                seen.add(item)
                not_duplicated.append(True)
            else:
                not_duplicated.append(False)

        not_duplicated = np.array(not_duplicated)

        if len(include_idx) != len(not_duplicated):
            include_idx = np.ones_like(not_duplicated, dtype=int)

        include_idx = include_idx * not_duplicated

    include_idx = include_idx.astype(bool)
    if not np.any(include_idx):
        return {"cs": None, "coverage": None, "requested_coverage": coverage}

    cs = [cs[i] for i in np.where(include_idx)[0]]
    claimed_coverage = claimed_coverage[include_idx]

    if Xcorr is None and X is None:
        cs_dict = {f"L{i}": cs[i] for i in range(len(cs))}
        return {
            "cs": cs_dict,
            "coverage": claimed_coverage,
            "requested_coverage": coverage,
        }
    else:
        purity = []
        for i, c in enumerate(cs):
            if null_index > 0 and null_index in c:
                purity.append([-9, -9, -9])
            else:
                purity.append(get_purity(c, X, Xcorr, squared, n_purity))  # type: ignore

        purity = np.array(purity)
        purity_df = {
            "min_sq_corr" if squared else "min_abs_corr": purity[:, 0],
            "mean_sq_corr" if squared else "mean_abs_corr": purity[:, 1],
            "median_sq_corr" if squared else "median_abs_corr": purity[:, 2],
        }

        threshold = min_abs_corr**2 if squared else min_abs_corr
        is_pure = np.where(purity[:, 0] >= threshold)[0]

        if len(is_pure) > 0:
            cs = [cs[i] for i in is_pure]
            purity = {k: v[is_pure] for k, v in purity_df.items()}
            row_names = [f"L{i}" for i in np.where(include_idx)[0][is_pure]]
            cs_dict = dict(zip(row_names, cs))

            # Re-order based on purity
            ordering = np.argsort(purity["min_sq_corr" if squared else "min_abs_corr"])[
                ::-1
            ]
            return {
                "cs": {row_names[i]: cs[i] for i in ordering},
                "purity": {k: v[ordering] for k, v in purity.items()},
                "cs_index": np.where(include_idx)[0][is_pure][ordering],
                "coverage": claimed_coverage[is_pure][ordering],
                "requested_coverage": coverage,
            }
        else:
            return {"cs": None, "coverage": None, "requested_coverage": coverage}


def in_CS_x(x: np.ndarray, coverage: float = 0.9) -> np.ndarray:
    """
    Determine credible set membership for a single component.

    Creates a binary indicator vector specifying which variants should be
    included in the credible set for a single SuSiE component, based on
    the specified coverage probability.

    Parameters
    ----------
    x : np.ndarray
        Vector of posterior inclusion probabilities (PIPs) for all variants.
        Should sum to approximately 1 for a well-calibrated component.
    coverage : float, optional
        Target coverage probability, by default 0.9.
        Variants are included until cumulative PIP reaches this threshold.

    Returns
    -------
    np.ndarray
        Binary vector of same length as x, where 1 indicates inclusion
        in the credible set and 0 indicates exclusion.

    Notes
    -----
    The algorithm:
    1. Ranks variants by decreasing posterior inclusion probability
    2. Includes top variants until cumulative probability ≥ coverage
    3. Returns binary vector indicating membership

    For well-calibrated components, the credible set should contain
    the true causal variant with probability ≥ coverage.

    Examples
    --------
    >>> # Simple example with 5 variants
    >>> pips = np.array([0.1, 0.6, 0.2, 0.05, 0.05])
    >>> cs_membership = in_CS_x(pips, coverage=0.8)
    >>> print(f"Credible set: variants {np.where(cs_membership)[0]}")
    >>> print(f"Coverage: {pips[cs_membership].sum():.2f}")
    Credible set: variants [1 2]
    Coverage: 0.80
    """
    n = n_in_CS_x(x, coverage)
    o = np.argsort(x)[::-1]
    result = np.zeros_like(x)
    result[o[:n]] = 1
    return result


def in_CS(
    res: Union[Dict[str, np.ndarray], np.ndarray], coverage: float = 0.9
) -> np.ndarray:
    """
    Determine credible set membership across all SuSiE components.

    Constructs credible sets for each component in a SuSiE model by applying
    the coverage criterion to the posterior inclusion probabilities.

    Parameters
    ----------
    res : Union[Dict[str, np.ndarray], np.ndarray]
        Either a SuSiE results dictionary containing 'alpha' key,
        or directly the L×P alpha matrix of posterior inclusion probabilities.
    coverage : float, optional
        Target coverage probability for each credible set, by default 0.9.

    Returns
    -------
    np.ndarray
        L×P binary matrix where entry (i,j) = 1 if variant j is included
        in the credible set for component i, and 0 otherwise.

    Notes
    -----
    This function applies in_CS_x to each row of the alpha matrix,
    creating credible sets for all L components simultaneously.

    The resulting matrix can be used to:
    - Identify which variants belong to each credible set
    - Calculate actual coverage achieved by each set
    - Check for overlaps between credible sets from different components

    Examples
    --------
    >>> # Extract credible sets from SuSiE results
    >>> cs_matrix = in_CS(susie_results, coverage=0.95)
    >>> print(f"Shape: {cs_matrix.shape}")  # L x P
    >>> print(f"Component 0 credible set size: {cs_matrix[0].sum()}")
    Shape: (5, 1000)
    Component 0 credible set size: 12

    >>> # Find variants in any credible set
    >>> any_cs = cs_matrix.any(axis=0)
    >>> print(f"Total variants in any credible set: {any_cs.sum()}")
    Total variants in any credible set: 45

    >>> # Check for overlapping credible sets
    >>> for i in range(cs_matrix.shape[0]):
    ...     for j in range(i+1, cs_matrix.shape[0]):
    ...         overlap = (cs_matrix[i] & cs_matrix[j]).sum()
    ...         if overlap > 0:
    ...             print(f"Components {i} and {j} overlap: {overlap} variants")
    """
    if isinstance(res, dict):
        res = res["alpha"]
    return np.apply_along_axis(lambda x: in_CS_x(x, coverage), 1, res)


def n_in_CS(
    res: Union[Dict[str, np.ndarray], np.ndarray], coverage: float = 0.9
) -> np.ndarray:
    """
    Compute the number of variables in each credible set.

    Parameters
    ----------
    res : Union[Dict[str, np.ndarray], np.ndarray]
        SuSiE results or alpha matrix.
    coverage : float, default 0.9
        The desired coverage.

    Returns
    -------
    np.ndarray
        Number of variables in each CS.
    """
    if isinstance(res, dict):
        res = res["alpha"]
    return np.apply_along_axis(lambda x: n_in_CS_x(x, coverage), 1, res)


def n_in_CS_x(x: np.ndarray, coverage: float = 0.9) -> int:
    """
    Compute the number of variables needed to achieve the desired coverage.

    Parameters
    ----------
    x : np.ndarray
        A probability vector.
    coverage : float, default 0.9
        The desired coverage.

    Returns
    -------
    int
        Number of variables needed.
    """
    cs = np.cumsum(np.sort(x)[::-1])
    return np.searchsorted(cs, coverage) + 1  # type: ignore


def get_purity(
    pos: np.ndarray,
    X: np.ndarray,
    Xcorr: Optional[np.ndarray],
    squared: bool = False,
    n: int = 100,
) -> np.ndarray:
    """
    Compute purity statistics for a credible set.

    Parameters
    ----------
    pos : np.ndarray
        Indices of variables in the credible set.
    X : np.ndarray, optional
        Data matrix.
    Xcorr : np.ndarray, optional
        Correlation matrix.
    squared : bool, default False
        If True, compute squared correlations.
    n : int, default 100
        Maximum number of variables to use.

    Returns
    -------
    np.ndarray
        Array of [min, mean, median] correlation statistics.
    """
    if len(pos) == 1:
        return np.array([1, 1, 1])

    if len(pos) > n:
        pos = np.random.choice(pos, n, replace=False)

    if Xcorr is None:
        X_sub = X[:, pos]
        value = np.abs(np.triu(np.corrcoef(X_sub.T), k=1))
    else:
        value = np.abs(np.triu(Xcorr[np.ix_(pos, pos)], k=1))

    value = value[value != 0]

    if squared:
        value = value**2

    return np.array([np.min(value), np.mean(value), np.median(value)])


def susie_get_objective(
    res: Dict[str, Union[List[float], np.ndarray]],
    last_only: bool = True,
    warning_tol: float = 1e-6,
) -> Union[float, np.ndarray]:
    """
    Get the evidence lower bound (ELBO) achieved by the fitted SuSiE model.

    Parameters
    ----------
    res : dict
        A dictionary containing SuSiE results, must include 'elbo' key.
    last_only : bool, default True
        If True, return only the ELBO from the last iteration.
    warning_tol : float, default 1e-6
        Tolerance level for warning about decreasing ELBO.

    Returns
    -------
    Union[float, np.ndarray]
        The ELBO value(s) from the SuSiE model.

    Raises
    ------
    Warning
        If the objective is decreasing beyond the specified tolerance.

    Notes
    -----
    This function returns the evidence lower bound (ELBO) achieved by
    the fitted SuSiE model. If last_only is False, it returns the ELBO
    from all iterations.
    """
    elbo = np.array(res["elbo"])
    if not np.all(np.diff(elbo) >= -warning_tol):
        logger.warning("Objective is decreasing")

    if last_only:
        return elbo[-1]
    else:
        return elbo


def susie_get_pip(
    res: Union[Dict[str, np.ndarray], np.ndarray],
    prune_by_cs: bool = False,
    prior_tol: float = 1e-9,
) -> np.ndarray:
    """
    Get posterior inclusion probabilities (PIPs) for all variables.

    Parameters
    ----------
    res : Union[Dict[str, np.ndarray], np.ndarray]
        A SuSiE fit result or an alpha matrix.
    prune_by_cs : bool, default False
        Whether to ignore single effects not in a reported CS when calculating PIP.
    prior_tol : float, default 1e-9
        Filter out effects having estimated prior variance smaller than this threshold.

    Returns
    -------
    np.ndarray
        A vector containing the posterior inclusion probabilities (PIPs) for all variables.

    Notes
    -----
    This function calculates the posterior inclusion probabilities (PIPs)
    for all variables based on the SuSiE model results.
    """
    if isinstance(res, dict):
        alpha = res["alpha"]

        # Drop null weight columns
        if "null_index" in res and res["null_index"] > 0:
            alpha = np.delete(alpha, res["null_index"] - 1, axis=1)

        # Drop the single-effects with estimated prior of zero
        if "V" in res and isinstance(res["V"], np.ndarray):
            include_idx = np.where(res["V"] > prior_tol)[0]
        else:
            include_idx = np.arange(alpha.shape[0])

        # Only consider variables in reported CS
        if prune_by_cs and "sets" in res and "cs_index" in res["sets"]:
            include_idx = np.intersect1d(include_idx, res["sets"]["cs_index"])
        elif prune_by_cs:
            include_idx = np.array([], dtype=int)

        # Extract relevant rows from alpha matrix
        if len(include_idx) > 0:
            alpha = alpha[include_idx]
        else:
            alpha = np.zeros((1, alpha.shape[1]))
    else:
        alpha = res

    return 1 - np.prod(1 - alpha, axis=0)


def estimate_residual_variance_ss(
    XtX: np.ndarray,
    Xty: np.ndarray,
    s: Dict[str, Union[np.ndarray, float]],
    yty: float,
    n: int,
) -> float:
    """
    Estimate residual variance for summary statistics.

    Parameters
    ----------
    XtX : np.ndarray
        A p by p matrix (X'X).
    Xty : np.ndarray
        A p vector (X'y).
    s : dict
        A SuSiE fit object.
    yty : float
        y'y, where y is centered to have mean 0.
    n : int
        Sample size.

    Returns
    -------
    float
        Estimated residual variance.
    """
    return (1 / n) * get_ER2_ss(XtX, Xty, s, yty)


def get_objective_ss(
    XtX: np.ndarray,
    Xty: np.ndarray,
    s: Dict[str, Union[np.ndarray, float]],
    yty: float,
    n: int,
) -> float:
    """
    Get objective function from data and SuSiE fit object.

    Parameters
    ----------
    XtX : np.ndarray
        A p by p matrix (X'X).
    Xty : np.ndarray
        A p vector (X'y).
    s : dict
        A SuSiE fit object.
    yty : float
        y'y, where y is centered to have mean 0.
    n : int
        Sample size.

    Returns
    -------
    float
        Objective function value.
    """
    return Eloglik_ss(XtX, Xty, s, yty, n) - np.sum(s["KL"])


def Eloglik_ss(
    XtX: np.ndarray,
    Xty: np.ndarray,
    s: Dict[str, Union[np.ndarray, float]],
    yty: float,
    n: int,
) -> float:
    """
    Compute expected log-likelihood for summary statistics.

    Parameters
    ----------
    XtX : np.ndarray
        A p by p matrix (X'X).
    Xty : np.ndarray
        A p vector (X'y).
    s : dict
        A SuSiE fit object.
    yty : float
        y'y, where y is centered to have mean 0.
    n : int
        Sample size.

    Returns
    -------
    float
        Expected log-likelihood.
    """
    return -n / 2 * np.log(2 * np.pi * s["sigma2"]) - 1 / (2 * s["sigma2"]) * get_ER2_ss(XtX, Xty, s, yty)  # type: ignore


def get_ER2_ss(
    XtX: np.ndarray, Xty: np.ndarray, s: Dict[str, Union[np.ndarray, float]], yty: float
) -> float:
    """
    Get expected squared residuals for summary statistics.

    Parameters
    ----------
    XtX : np.ndarray
        A p by p matrix (X'X).
    Xty : np.ndarray
        A p vector (X'y).
    s : dict
        A SuSiE fit object.
    yty : float
        y'y, where y is centered to have mean 0.

    Returns
    -------
    float
        Expected squared residuals.
    """
    B = s["alpha"] * s["mu"]
    XB2 = np.sum((B @ XtX) * B)
    betabar = np.sum(B, axis=0)
    d = np.diag(XtX)
    postb2 = s["alpha"] * s["mu2"]  # Posterior second moment

    result = (
        yty
        - 2 * np.sum(betabar * Xty)
        + np.sum(betabar * (XtX @ betabar))
        - XB2
        + np.sum(d * postb2)
    )  # .T)
    return result


def SER_posterior_e_loglik_ss(
    dXtX: np.ndarray, Xty: np.ndarray, s2: float, Eb: np.ndarray, Eb2: np.ndarray
) -> float:
    """
    Posterior expected log-likelihood for a single effect regression.

    Parameters
    ----------
    dXtX : np.ndarray
        A p vector of diagonal elements of XtX.
    Xty : np.ndarray
        A p vector.
    s2 : float
        The residual variance.
    Eb : np.ndarray
        The posterior mean of b (p vector) (alpha * mu).
    Eb2 : np.ndarray
        The posterior second moment of b (p vector) (alpha * mu2).

    Returns
    -------
    float
        Posterior expected log-likelihood.
    """
    return -0.5 / s2 * (-2 * np.sum(Eb * Xty) + np.sum(dXtX * Eb2))  # type: ignore


def est_V_uniroot(
    betahat: np.ndarray, shat2: np.ndarray, prior_weights: np.ndarray
) -> float:
    """
    Estimate prior variance using uniroot method.

    Parameters
    ----------
    betahat : np.ndarray
        Estimated beta coefficients.
    shat2 : np.ndarray
        Estimated variances.
    prior_weights : np.ndarray
        Prior weights.

    Returns
    -------
    float
        Estimated prior variance.
    """
    V_u = optimize.root_scalar(
        negloglik_grad_logscale,
        args=(betahat, shat2, prior_weights),
        bracket=(-10, 10),
        method="brentq",
    )
    return np.exp(V_u.root)


def neg_loglik_logscale(
    lV: float, betahat: np.ndarray, shat2: np.ndarray, prior_weights: np.ndarray
) -> float:
    """
    Negative log-likelihood on log scale.

    Parameters
    ----------
    lV : float
        Log of prior variance.
    betahat : np.ndarray
        Estimated beta coefficients.
    shat2 : np.ndarray
        Estimated variances.
    prior_weights : np.ndarray
        Prior weights.

    Returns
    -------
    float
        Negative log-likelihood.
    """
    return -loglik(np.exp(lV), betahat, shat2, prior_weights)


def loglik_grad(
    V: float, betahat: np.ndarray, shat2: np.ndarray, prior_weights: np.ndarray
) -> float:
    """
    Gradient of log-likelihood.

    Parameters
    ----------
    V : float
        Prior variance.
    betahat : np.ndarray
        Estimated beta coefficients.
    shat2 : np.ndarray
        Estimated variances.
    prior_weights : np.ndarray
        Prior weights.

    Returns
    -------
    float
        Gradient of log-likelihood.
    """
    lbf = stats.norm.logpdf(betahat, 0, np.sqrt(V + shat2)) - stats.norm.logpdf(
        betahat, 0, np.sqrt(shat2)
    )
    lpo = lbf + np.log(prior_weights + np.sqrt(np.finfo(float).eps))

    # Deal with special case of infinite shat2
    lbf[np.isinf(shat2)] = 0
    lpo[np.isinf(shat2)] = 0

    maxlpo = np.max(lpo)
    w_weighted = np.exp(lpo - maxlpo)
    weighted_sum_w = np.sum(w_weighted)
    alpha = w_weighted / weighted_sum_w
    return np.sum(alpha * lbf_grad(V, shat2, betahat**2 / shat2))


def negloglik_grad_logscale(
    lV: float, betahat: np.ndarray, shat2: np.ndarray, prior_weights: np.ndarray
) -> float:
    """
    Negative gradient of log-likelihood on log scale.

    Parameters
    ----------
    lV : float
        Log of prior variance.
    betahat : np.ndarray
        Estimated beta coefficients.
    shat2 : np.ndarray
        Estimated variances.
    prior_weights : np.ndarray
        Prior weights.

    Returns
    -------
    float
        Negative gradient of log-likelihood.
    """
    return -np.exp(lV) * loglik_grad(np.exp(lV), betahat, shat2, prior_weights)


def lbf_grad(V: float, shat2: np.ndarray, T2: np.ndarray) -> np.ndarray:
    """
    Gradient of log Bayes factor.

    Parameters
    ----------
    V : float
        Prior variance.
    shat2 : np.ndarray
        Estimated variances.
    T2 : np.ndarray
        Squared t-statistics.

    Returns
    -------
    np.ndarray
        Gradient of log Bayes factor.
    """
    lbf = 0.5 * (1 / (V + shat2)) * ((shat2 / (V + shat2)) * T2 - 1)
    lbf[np.isnan(lbf)] = 0
    return lbf


def loglik(
    V: float, betahat: np.ndarray, shat2: np.ndarray, prior_weights: np.ndarray
) -> float:
    """
    Log-likelihood function for SER model.

    Parameters
    ----------
    V : float
        Prior variance.
    betahat : np.ndarray
        Estimated beta coefficients.
    shat2 : np.ndarray
        Estimated variances.
    prior_weights : np.ndarray
        Prior weights.

    Returns
    -------
    float
        Log-likelihood.
    """
    lbf = stats.norm.logpdf(betahat, 0, np.sqrt(V + shat2)) - stats.norm.logpdf(
        betahat, 0, np.sqrt(shat2)
    )
    lpo = lbf + np.log(prior_weights + np.sqrt(np.finfo(float).eps))

    # Deal with special case of infinite shat2
    lbf[np.isinf(shat2)] = 0
    lpo[np.isinf(shat2)] = 0

    maxlpo = np.max(lpo)
    w_weighted = np.exp(lpo - maxlpo)
    weighted_sum_w = np.sum(w_weighted)
    return np.log(weighted_sum_w) + maxlpo


def optimize_prior_variance(
    optimize_V: str,
    betahat: np.ndarray,
    shat2: np.ndarray,
    prior_weights: np.ndarray,
    alpha: Optional[np.ndarray] = None,
    post_mean2: Optional[np.ndarray] = None,
    V_init: Optional[float] = None,
    check_null_threshold: float = 0,
) -> float:
    """
    Optimize the prior variance using different methods.

    Parameters
    ----------
    optimize_V : str
        Method to optimize V. Options are "simple", "optim", "uniroot", or "EM".
    betahat : np.ndarray
        Estimated beta coefficients.
    shat2 : np.ndarray
        Estimated variances.
    prior_weights : np.ndarray
        Prior weights.
    alpha : np.ndarray, optional
        Alpha values for EM method.
    post_mean2 : np.ndarray, optional
        Posterior mean squared for EM method.
    V_init : float, optional
        Initial value for V.
    check_null_threshold : float, default 0
        Threshold for setting V to zero.

    Returns
    -------
    float
        Optimized prior variance.

    Raises
    ------
    ValueError
        If an invalid optimization method is specified.
    """
    V = V_init

    if optimize_V != "simple":
        if optimize_V == "optim":
            result = optimize.minimize_scalar(
                neg_loglik_logscale,
                args=(betahat, shat2, prior_weights),
                method="bounded",
                bounds=(-30, 15),
            )
            lV = result.x  # type: ignore

            # If the estimated one is worse than current one, don't change it
            if neg_loglik_logscale(
                lV, betahat, shat2, prior_weights
            ) > neg_loglik_logscale(
                np.log(V) if V != 0 else -np.inf, betahat, shat2, prior_weights  # type: ignore
            ):
                lV = np.log(V) if V != 0 else -np.inf  # type: ignore
            V = np.exp(lV)

        elif optimize_V == "uniroot":
            V = est_V_uniroot(betahat, shat2, prior_weights)

        elif optimize_V == "EM":
            if alpha is None or post_mean2 is None:
                raise ValueError("Alpha and post_mean2 must be provided for EM method")
            V = np.sum(alpha * post_mean2)

        else:
            raise ValueError("Invalid option for optimize_V method")

    # Set V exactly 0 if that beats the numerical value by check_null_threshold in loglik
    if loglik(0, betahat, shat2, prior_weights) + check_null_threshold >= loglik(V, betahat, shat2, prior_weights):  # type: ignore
        V = 0

    return V  # type: ignore


def single_effect_regression_ss(
    Xty: np.ndarray,
    dXtX: np.ndarray,
    V: float = 1,
    residual_variance: float = 1,
    prior_weights: Optional[np.ndarray] = None,
    optimize_V: str = "none",
    check_null_threshold: float = 0,
) -> Dict[str, Any]:
    """
    Perform single effect regression using summary statistics.

    Parameters
    ----------
    Xty : np.ndarray
        A p-vector.
    dXtX : np.ndarray
        A p-vector containing the diagonal elements of X'X.
    V : float, default 1
        Prior variance.
    residual_variance : float, default 1
        Residual variance.
    prior_weights : np.ndarray, optional
        Prior weights for each variable.
    optimize_V : str, default "none"
        Method to optimize V. Options are "none", "optim", "uniroot", "EM", or "simple".
    check_null_threshold : float, default 0
        Threshold for setting V to zero.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing regression results:
        - alpha: posterior inclusion probabilities
        - mu: posterior mean
        - mu2: posterior second moment
        - lbf: log Bayes factors
        - V: optimized prior variance
        - lbf_model: log Bayes factor for the model

    Raises
    ------
    ValueError
        If an invalid optimization method is specified.
    """
    valid_optimize_V = ["none", "optim", "uniroot", "EM", "simple"]
    if optimize_V not in valid_optimize_V:
        raise ValueError(
            f"Invalid optimize_V method. Must be one of {valid_optimize_V}"
        )

    betahat = Xty / dXtX
    shat2 = residual_variance / dXtX

    if prior_weights is None:
        prior_weights = np.ones_like(dXtX) / len(dXtX)

    if optimize_V != "EM" and optimize_V != "none":
        V = optimize_prior_variance(
            optimize_V,
            betahat,
            shat2,
            prior_weights,
            V_init=V,
            check_null_threshold=check_null_threshold,
        )

    # log(po) = log(BF * prior) for each SNP
    lbf = stats.norm.logpdf(betahat, 0, np.sqrt(V + shat2)) - stats.norm.logpdf(
        betahat, 0, np.sqrt(shat2)
    )
    lpo = lbf + np.log(prior_weights + np.sqrt(np.finfo(float).eps))

    # Deal with special case of infinite shat2
    lbf[np.isinf(shat2)] = 0
    lpo[np.isinf(shat2)] = 0
    maxlpo = np.max(lpo)

    # w is proportional to posterior odds = BF * prior,
    # but subtract max for numerical stability
    w_weighted = np.exp(lpo - maxlpo)
    weighted_sum_w = np.sum(w_weighted)
    alpha = w_weighted / weighted_sum_w
    if V == 0:
        # post_var is all zeros if V=0, shape is the same as alpha
        post_var = np.zeros_like(alpha)
    else:
        post_var = 1 / (1 / V + dXtX / residual_variance)  # Posterior variance
    post_mean = (1 / residual_variance) * post_var * Xty
    post_mean2 = post_var + post_mean**2  # Second moment
    lbf_model = maxlpo + np.log(
        weighted_sum_w
    )  # Analogue of loglik in the non-summary case

    if optimize_V == "EM":
        V = optimize_prior_variance(
            optimize_V,
            betahat,
            shat2,
            prior_weights,
            alpha=alpha,
            post_mean2=post_mean2,
            check_null_threshold=check_null_threshold,
        )

    return {
        "alpha": alpha,
        "mu": post_mean,
        "mu2": post_mean2,
        "lbf": lbf,
        "V": V,
        "lbf_model": lbf_model,
    }


def update_each_effect_ss(
    XtX: np.ndarray,
    Xty: np.ndarray,
    s_init: Dict[str, Any],
    estimate_prior_variance: bool = False,
    estimate_prior_method: str = "optim",
    check_null_threshold: float = 0,
) -> Dict[str, Any]:
    """
    Update each effect once.

    Parameters
    ----------
    XtX : np.ndarray
        A p by p matrix, X'X.
    Xty : np.ndarray
        A p vector.
    s_init : Dict[str, Any]
        A dictionary with elements sigma2, V, alpha, mu, XtXr.
    estimate_prior_variance : bool, default False
        Boolean indicating whether to estimate prior variance.
    estimate_prior_method : str, default "optim"
        The method used for estimating prior variance, 'optim' or 'EM'.
    check_null_threshold : float, default 0
        A threshold on the log scale to compare likelihood between current estimate and zero the null.

    Returns
    -------
    Dict[str, Any]
        Updated dictionary of model parameters.
    """
    if not estimate_prior_variance:
        estimate_prior_method = "none"

    s = s_init.copy()
    L = s["alpha"].shape[0]

    if L > 0:
        for effect_index in range(L):
            # Remove effect_index-th effect from fitted values
            s["XtXr"] = s["XtXr"] - XtX @ (
                s["alpha"][effect_index] * s["mu"][effect_index]
            )
            logger.debug(f"single_effect_regression_ss: {effect_index}")
            # Compute residuals
            XtR = Xty - s["XtXr"]
            res = single_effect_regression_ss(
                XtR,
                np.diag(XtX),
                s["V"][effect_index],
                s["sigma2"],
                s["pi"],
                estimate_prior_method,
                check_null_threshold,
            )

            # Update the variational estimate of the posterior mean
            s["mu"][effect_index] = res["mu"]
            s["alpha"][effect_index] = res["alpha"]
            s["mu2"][effect_index] = res["mu2"]
            s["V"][effect_index] = res["V"]
            s["lbf"][effect_index] = res["lbf_model"]
            s["lbf_variable"][effect_index] = res["lbf"]
            s["KL"][effect_index] = -res["lbf_model"] + SER_posterior_e_loglik_ss(
                np.diag(XtX),
                XtR,
                s["sigma2"],
                res["alpha"] * res["mu"],
                res["alpha"] * res["mu2"],
            )

            s["XtXr"] = s["XtXr"] + XtX @ (
                s["alpha"][effect_index] * s["mu"][effect_index]
            )

    s["XtXr"] = s["XtXr"].astype(float)
    return s


def susie_slim(res: Dict[str, Union[np.ndarray, int, float]]) -> Dict[str, Any]:
    """
    Slim the result of fitted SuSiE model.

    Parameters
    ----------
    res : Dict[str, Union[np.ndarray, int, float]]
        The result of a fitted SuSiE model.

    Returns
    -------
    Dict[str, Union[np.ndarray, int, float]]
        A slimmed version of the SuSiE model result.
    """
    return {
        "alpha": res["alpha"],
        "niter": res["niter"],
        "V": res["V"],
        "sigma2": res["sigma2"],
    }


def susie_prune_single_effects(
    s: Dict[str, Any],
    L: int = 0,
    V: Optional[Union[float, np.ndarray]] = None,
) -> Dict[str, Union[np.ndarray, int, float, Dict]]:
    """
    Prune single effects to given number L in SuSiE model object.

    Parameters
    ----------
    s : Dict[str, Union[np.ndarray, int, float, Dict]]
        A SuSiE model object.
    L : int, default 0
        The number of effects to keep. If 0, it will be determined based on non-zero elements in s['V'].
    V : Optional[Union[float, np.ndarray]], default None
        Prior variances for the effects.

    Returns
    -------
    Dict[str, Union[np.ndarray, int, float, Dict]]
        The pruned SuSiE model object.
    """
    num_effects = s["alpha"].shape[0]

    if L == 0:
        # Filtering will be based on non-zero elements in s['V']
        if s.get("V") is not None:
            L = np.sum(s["V"] > 0)
        else:
            L = num_effects

    if L == num_effects:
        s["sets"] = None
        return s

    if s.get("sets", {}).get("cs_index") is not None:
        effects_rank = np.concatenate(
            [
                s["sets"]["cs_index"],
                np.setdiff1d(np.arange(num_effects), s["sets"]["cs_index"]),
            ]
        )
    else:
        effects_rank = np.arange(num_effects)

    if L > num_effects:
        logger.warning(
            f"Specified number of effects L = {L} is greater than the number of effects {num_effects} "
            f"in input SuSiE model. The SuSiE model will be expanded to have {L} effects."
        )

        s["alpha"] = np.vstack(
            [
                s["alpha"][effects_rank],
                np.full(
                    (L - num_effects, s["alpha"].shape[1]), 1 / s["alpha"].shape[1]
                ),
            ]
        )

        for n in ["mu", "mu2", "lbf_variable"]:
            if n in s:
                s[n] = np.vstack(
                    [s[n][effects_rank], np.zeros((L - num_effects, s[n].shape[1]))]
                )

        for n in ["KL", "lbf"]:
            if n in s:
                s[n] = np.concatenate(
                    [s[n][effects_rank], np.full(L - num_effects, np.nan)]
                )

        if V is not None:
            if isinstance(V, np.ndarray) and V.size > 1:
                V[:num_effects] = s["V"][effects_rank]
            else:
                V = np.full(L, V)
        s["V"] = V
    else:
        # Prune to L effects
        s["alpha"] = s["alpha"][effects_rank[:L]]
        for n in ["mu", "mu2", "lbf_variable"]:
            if n in s:
                s[n] = s[n][effects_rank[:L]]
        for n in ["KL", "lbf", "V"]:
            if n in s:
                s[n] = s[n][effects_rank[:L]]

    s["sets"] = None
    return s


def init_setup(
    n: int,
    p: int,
    L: int,
    scaled_prior_variance: float,
    residual_variance: float,
    prior_weights: np.ndarray,
    null_weight: float,
    varY: float,
    standardize: bool,
) -> Dict[str, Any]:
    """
    Set default SuSiE initialization.

    Parameters
    ----------
    n : int
        Number of samples.
    p : int
        Number of variables.
    L : int
        Number of effects.
    scaled_prior_variance : float
        Scaled prior variance.
    residual_variance : float
        Residual variance.
    prior_weights : np.ndarray
        Prior weights for each variable.
    null_weight : float
        Weight for the null component.
    varY : float
        Variance of Y.
    standardize : bool
        Whether the data is standardized.

    Returns
    -------
    Dict[str, Union[np.ndarray, float, int]]
        Initialized SuSiE model object.

    Raises
    ------
    ValueError
        If input parameters are invalid.
    """
    if not isinstance(scaled_prior_variance, (int, float)) or scaled_prior_variance < 0:
        raise ValueError("Scaled prior variance should be positive number")
    if scaled_prior_variance > 1 and standardize:
        raise ValueError(
            "Scaled prior variance should be no greater than 1 when standardize = True"
        )

    if residual_variance is None:
        residual_variance = varY

    if prior_weights is None:
        prior_weights = np.full(p, 1 / p)
    else:
        if np.all(prior_weights == 0):
            raise ValueError(
                "Prior weight should be greater than 0 for at least one variable."
            )
        prior_weights = prior_weights / np.sum(prior_weights)

    if len(prior_weights) != p:
        raise ValueError("Prior weights must have length p")

    if p < L:
        L = p

    s = {
        "alpha": np.full((L, p), 1 / p),
        "mu": np.zeros((L, p)),
        "mu2": np.zeros((L, p)),
        "Xr": np.zeros(n),
        "KL": np.full(L, np.nan),
        "lbf": np.full(L, np.nan),
        "lbf_variable": np.full((L, p), np.nan),
        "sigma2": residual_variance,
        "V": [scaled_prior_variance * varY],
        "pi": prior_weights,
    }

    if null_weight is None:
        s["null_index"] = 0
    else:
        s["null_index"] = p

    return s


def init_finalize(
    s: Dict[str, Any],
    X: Optional[np.ndarray] = None,
    Xr: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Update a SuSiE fit object in order to initialize SuSiE model.

    Parameters
    ----------
    s : Dict[str, Union[np.ndarray, float, int]]
        SuSiE model object.
    X : Optional[np.ndarray]
        Design matrix.
    Xr : Optional[np.ndarray]
        Residuals.

    Returns
    -------
    Dict[str, Union[np.ndarray, float, int]]
        Updated SuSiE model object.

    Raises
    ------
    ValueError
        If input parameters are invalid.
    """
    if len(s["V"]) == 1:
        s["V"] = np.full(s["alpha"].shape[0], s["V"])

    # Check sigma2
    if not isinstance(s["sigma2"], (int, float)):
        raise ValueError("Input residual variance sigma2 must be numeric")
    s["sigma2"] = float(s["sigma2"])
    if s["sigma2"] <= 0:
        raise ValueError(
            "Residual variance sigma2 must be positive (is your var(Y) zero?)"
        )

    # Check prior variance
    if not isinstance(s["V"], np.ndarray):
        raise ValueError("Input prior variance must be numeric")
    if not np.all(s["V"] >= 0):
        raise ValueError("Prior variance must be non-negative")
    if s["mu"].shape != s["mu2"].shape:
        raise ValueError("Dimension of mu and mu2 in input object do not match")
    if s["mu"].shape != s["alpha"].shape:
        raise ValueError("Dimension of mu and alpha in input object do not match")
    if s["alpha"].shape[0] != len(s["V"]):
        raise ValueError(
            "Input prior variance V must have length of nrow of alpha in input object"
        )

    # Update Xr
    if Xr is not None:
        s["Xr"] = Xr
    elif X is not None:
        s["Xr"] = compute_Xb(X, np.sum(s["mu"] * s["alpha"], axis=0))

    # Reset KL and lbf
    s["KL"] = np.full(s["alpha"].shape[0], np.nan)
    s["lbf"] = np.full(s["alpha"].shape[0], np.nan)

    return s


def compute_tf_Xb(order: int, b: np.ndarray) -> np.ndarray:
    """
    Compute unscaled X @ b using the special structure of trend filtering.

    This function applies the trend filtering operation of a given order
    to a vector b.

    Parameters
    ----------
    order : int
        The order of trend filtering.
    b : np.ndarray
        An n-vector to be filtered.

    Returns
    -------
    np.ndarray
        An n-vector resulting from the trend filtering operation.

    Notes
    -----
    This function applies the trend filtering operation iteratively,
    where each iteration involves a reverse cumulative sum operation.
    """
    for _ in range(order + 1):
        b = -np.cumsum(b[::-1])[::-1]
    return b


def compute_Xb(X: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute standardized X @ b using sparse multiplication trick.

    This function computes the product of a standardized matrix X and a vector b,
    where X is provided in its unstandardized form along with standardization parameters.

    Parameters
    ----------
    X : np.ndarray
        An n by p unstandardized matrix. It should have the following attributes:
        - "scaled:center": The centering vector used for standardization.
        - "scaled:scale": The scaling vector used for standardization.
        - "matrix.type": Optional. If present, X is treated as a trend filtering matrix.
        - "order": Required if X is a trend filtering matrix. The order of the trend filtering.

    b : np.ndarray
        A p-vector to be multiplied with X.

    Returns
    -------
    np.ndarray
        An n-vector resulting from the multiplication of standardized X and b.

    Notes
    -----
    This function assumes that X has been unstandardized and that the standardization
    parameters are provided as attributes. It performs the multiplication efficiently
    by avoiding explicit standardization of X.
    """
    cm = getattr(X, "scaled:center", None)
    csd = getattr(X, "scaled:scale", None)

    if cm is None or csd is None:
        raise ValueError("X must have 'scaled:center' and 'scaled:scale' attributes")

    # Scale Xb
    if hasattr(X, "matrix.type"):
        # When X is a trend filtering matrix
        order = getattr(X, "order", None)
        if order is None:
            raise ValueError("Trend filtering matrix X must have an 'order' attribute")
        scaled_Xb = compute_tf_Xb(order, b / csd)
    else:
        # When X is an ordinary sparse/dense matrix
        if sparse.issparse(X):
            scaled_Xb = X.dot(b / csd)
        else:
            scaled_Xb = X @ (b / csd)

    # Center Xb
    Xb = scaled_Xb - np.sum(cm * b / csd)

    return Xb


def muffled_cov2cor(x):
    """
    Calculate the correlation matrix from a covariance matrix.

    This function calculates the correlation matrix from a covariance matrix,
    suppressing warnings related to division by zero.

    Parameters
    ----------
    x : np.ndarray
        A p by p covariance matrix.

    Returns
    -------
    np.ndarray
        A p by p correlation matrix.

    Notes
    -----
    This function suppresses warnings related to division by zero when calculating
    the correlation matrix from a covariance matrix.
    """

    def custom_formatwarning(msg, *args, **kwargs):
        if "invalid value encountered in true_divide" in str(msg):
            return ""
        return warnings.formatwarning(msg, *args, **kwargs)

    original_formatwarning = warnings.formatwarning
    warnings.formatwarning = custom_formatwarning

    try:
        # 计算相关矩阵
        with warnings.catch_warnings():
            warnings.simplefilter("always")
            std = np.sqrt(np.diag(x))
            corr = x / np.outer(std, std)
            corr[~np.isfinite(corr)] = 0
        return corr
    finally:
        warnings.formatwarning = original_formatwarning


def susie_suff_stat(
    XtX: np.ndarray,
    Xty: np.ndarray,
    yty: float,
    n: int,
    X_colmeans: Optional[np.ndarray] = None,
    y_mean: Optional[float] = None,
    L: int = 10,
    scaled_prior_variance: float = 0.2,
    residual_variance: Optional[float] = None,
    estimate_residual_variance: bool = True,
    estimate_prior_variance: bool = True,
    estimate_prior_method: str = "optim",
    check_null_threshold: float = 0,
    prior_tol: float = 1e-9,
    r_tol: float = 1e-8,
    prior_weights: Optional[np.ndarray] = None,
    null_weight: float = 0,
    standardize: bool = True,
    max_iter: int = 100,
    s_init: Optional[Dict[str, Any]] = None,
    coverage: float = 0.95,
    min_abs_corr: float = 0.5,
    tol: float = 1e-3,
    track_fit: bool = False,
    refine: bool = False,
    check_prior: bool = False,
    n_purity: int = 100,
) -> Dict[str, Any]:
    """
    Perform Sum of Single Effects (SuSiE) regression using sufficient statistics.

    Parameters
    ----------
    XtX : np.ndarray
        A p by p matrix X'X where columns of X are centered to have mean zero.
    Xty : np.ndarray
        A p-vector X'y where y and columns of X are centered to have mean zero.
    yty : float
        A scalar y'y where y is centered to have mean zero.
    n : int
        The sample size.
    X_colmeans : np.ndarray, optional
        A p-vector of column means of X.
    y_mean : float, optional
        A scalar containing the mean of y.
    L : int, default=10
        Number of effects to fit.
    scaled_prior_variance : float, default=0.2
        Prior variance of effects.
    residual_variance : float, optional
        Residual variance. If None, it will be estimated.
    estimate_residual_variance : bool, default=True
        Whether to estimate residual variance.
    estimate_prior_variance : bool, default=True
        Whether to estimate prior variance.
    estimate_prior_method : str, default="optim"
        Method to estimate prior variance. Options: "optim", "EM", "simple".
    check_null_threshold : float, default=0
        Threshold for null model check.
    prior_tol : float, default=1e-9
        Tolerance for prior calculations.
    r_tol : float, default=1e-8
        Tolerance for eigenvalue check of positive semidefinite matrix.
    prior_weights : np.ndarray, optional
        Prior weights for each variable.
    null_weight : float, default=0
        Weight for the null model.
    standardize : bool, default=True
        Whether to standardize the variables.
    max_iter : int, default=100
        Maximum number of iterations.
    s_init : dict, optional
        Initial SuSiE fit to start from.
    coverage : float, default=0.95
        Coverage for credible sets.
    min_abs_corr : float, default=0.5
        Minimum absolute correlation for credible sets.
    tol : float, default=1e-3
        Convergence tolerance.
    track_fit : bool, default=False
        Whether to track the fit at each iteration.
    refine : bool, default=False
        Whether to refine the fit.
    check_prior : bool, default=False
        Whether to check if the estimated prior variance is reasonable.
    n_purity : int, default=100
        Number of samples for purity calculations.

    Returns
    -------
    dict
        A dictionary containing the SuSiE fit results.

    Raises
    ------
    ValueError
        If inputs are invalid or inconsistent.

    Notes
    -----
    This function is a complex implementation of the SuSiE algorithm. It includes
    many options and features that may not be fully documented here. Please refer
    to the original R implementation and associated documentation for more details.
    """
    # Input validation
    if XtX.shape[1] != len(Xty):
        raise ValueError(
            f"The dimension of XtX ({XtX.shape[0]} by {XtX.shape[1]}) "
            f"does not agree with expected ({len(Xty)} by {len(Xty)})"
        )

    # Ensure XtX is symmetric
    if not np.allclose(XtX, XtX.T):
        logger.warning("XtX is not symmetric; forcing XtX to be symmetric")
        XtX = (XtX + XtX.T) / 2

    if np.isinf(Xty).any():
        raise ValueError("Input Xty contains infinite values")

    if np.isnan(XtX).any():
        raise ValueError("Input XtX matrix contains NAs")

    # Replace NAs in Xty with zeros
    if np.isnan(Xty).any():
        logger.warning("NA values in Xty are replaced with 0")
        Xty = np.nan_to_num(Xty)

    # Process null_weight
    if null_weight == 0:
        null_weight = None  # type: ignore
    if null_weight is not None:
        if not 0 <= null_weight < 1:
            raise ValueError("Null weight must be between 0 and 1")
        if prior_weights is None:
            prior_weights = np.full(XtX.shape[1], (1 - null_weight) / XtX.shape[1])
        else:
            prior_weights = np.concatenate(
                [prior_weights * (1 - null_weight), [null_weight]]
            )
        XtX = np.pad(XtX, ((0, 1), (0, 1)))
        Xty = np.append(Xty, 0)

    p = XtX.shape[1]

    # Standardization
    if standardize:
        dXtX = np.diag(XtX)
        csd = np.sqrt(dXtX / (n - 1))
        csd[csd == 0] = 1
        XtX = (XtX / csd).T / csd
        Xty = Xty / csd
    else:
        csd = np.ones(p)

    # Check X_colmeans
    if X_colmeans is None:
        X_colmeans = np.full(p, np.nan)
    elif np.isscalar(X_colmeans):
        X_colmeans = np.full(p, X_colmeans)
    elif len(X_colmeans) != p:
        raise ValueError(
            "The length of X_colmeans does not agree with number of variables"
        )

    # Initialize SuSiE fit
    s = init_setup(
        0,
        p,
        L,
        scaled_prior_variance,
        residual_variance,  # type: ignore
        prior_weights,
        null_weight,
        yty / (n - 1),
        standardize,
    )
    s["Xr"] = None
    s["XtXr"] = np.zeros(p)

    # Process s_init if provided
    if s_init is not None:
        if not isinstance(s_init, dict):
            raise ValueError("s_init should be a susie object")
        if np.max(s_init["alpha"]) > 1 or np.min(s_init["alpha"]) < 0:
            raise ValueError(
                "s_init['alpha'] has invalid values outside range [0,1]; please check your input"
            )

        # First, remove effects with s_init['V'] = 0
        s_init = susie_prune_single_effects(s_init)
        num_effects = s_init["alpha"].shape[0]
        if L is None:
            L = num_effects
        elif min(p, L) < num_effects:
            logger.warning(
                f"Specified number of effects L = {min(p, L)} is smaller than the number of effects "
                f"{num_effects} in input SuSiE model. The initialized SuSiE model will have {num_effects} effects."
            )
            L = num_effects
        # Expand s_init if L > num_effects
        s_init = susie_prune_single_effects(s_init, min(p, L), s["V"])
        s.update(s_init)
        s = init_finalize(s, XtX)
        s["XtXr"] = s["Xr"]
        s["Xr"] = None  # type: ignore
    else:
        s = init_finalize(s)

    # Main iteration loop
    elbo = np.full(max_iter + 1, np.nan)
    elbo[0] = -np.inf
    tracking: List[Dict[str, Any]] = []

    bhat = Xty / np.diag(XtX)
    shat = np.sqrt(s["sigma2"] / np.diag(XtX))
    z = bhat / shat
    zm = np.max(np.abs(z[~np.isnan(z)]))

    for i in range(max_iter):
        if track_fit:
            tracking.append(susie_slim(s))
        logger.info(f"Iteration {i}/{max_iter}")
        s = update_each_effect_ss(
            XtX,
            Xty,
            s,
            estimate_prior_variance,
            estimate_prior_method,
            check_null_threshold,
        )

        if check_prior:
            if np.any(s["V"] > 100 * (zm**2)):
                raise ValueError(
                    "The estimated prior variance is unreasonably large. "
                    "This is usually caused by mismatch between the summary "
                    "statistics and the LD matrix. Please check the input."
                )

        logger.debug(f"objective: {get_objective_ss(XtX, Xty, s, yty, n)}")

        elbo[i + 1] = get_objective_ss(XtX, Xty, s, yty, n)
        if np.isinf(elbo[i + 1]):
            raise ValueError("The objective becomes infinite. Please check the input.")

        if (elbo[i + 1] - elbo[i]) < tol:
            logger.info("Converged, stopping iteration")
            s["converged"] = True
            break

        if estimate_residual_variance:
            est_sigma2 = estimate_residual_variance_ss(XtX, Xty, s, yty, n)
            if est_sigma2 < 0:
                raise ValueError(
                    "Estimating residual variance failed: the estimated value is negative"
                )
            s["sigma2"] = est_sigma2
            logger.debug(f"objective: {get_objective_ss(XtX, Xty, s, yty, n)}")

    elbo = elbo[1 : i + 2]  # type: ignore
    s["elbo"] = elbo
    s["niter"] = i + 1

    if "converged" not in s:
        logger.warning(
            f"IBSS algorithm did not converge in {max_iter} iterations! "
            "Please check consistency between summary statistics and LD matrix."
        )
        s["converged"] = False

    s["X_column_scale_factors"] = csd

    # Compute intercept
    if y_mean is None:
        s["intercept"] = np.nan
    else:
        s["intercept"] = y_mean - np.sum(
            X_colmeans
            * (np.sum(s["alpha"] * s["mu"], axis=0) / s["X_column_scale_factors"])
        )

    if track_fit:
        s["trace"] = tracking
    # SuSiE CS and PIP
    if coverage is not None and min_abs_corr is not None:
        if not np.all(np.isin(np.diag(XtX), [0, 1])):
            s["sets"] = susie_get_cs(
                s,
                coverage=coverage,
                Xcorr=muffled_cov2cor(XtX),
                min_abs_corr=min_abs_corr,
                check_symmetric=False,
                n_purity=n_purity,
            )
        else:
            s["sets"] = susie_get_cs(
                s,
                coverage=coverage,
                Xcorr=XtX,
                min_abs_corr=min_abs_corr,
                check_symmetric=False,
                n_purity=n_purity,
            )
        s["pip"] = susie_get_pip(s, prune_by_cs=False, prior_tol=prior_tol)

    if refine:
        if s_init is not None:
            logger.warning("The given s_init is not used in refinement")

        if null_weight is not None and null_weight != 0:
            # If null_weight is specified, we remove the extra 0 column
            XtX = XtX[:-1, :-1]
            Xty = Xty[:-1]
            pw_s = s["pi"][:-1] / (1 - null_weight)
        else:
            pw_s = s["pi"]

        conti = True
        while conti and len(s["sets"]["cs"]) > 0:
            m = []
            for cs in range(len(s["sets"]["cs"])):
                pw_cs = pw_s.copy()
                pw_cs[s["sets"]["cs"][cs]] = 0
                if np.all(pw_cs == 0):
                    break

                s2 = susie_suff_stat(
                    XtX=XtX,
                    Xty=Xty,
                    yty=yty,
                    n=n,
                    L=L,
                    X_colmeans=X_colmeans,
                    y_mean=y_mean,
                    prior_weights=pw_cs,
                    s_init=None,
                    scaled_prior_variance=scaled_prior_variance,
                    residual_variance=residual_variance,
                    estimate_residual_variance=estimate_residual_variance,
                    estimate_prior_variance=estimate_prior_variance,
                    estimate_prior_method=estimate_prior_method,
                    check_null_threshold=check_null_threshold,
                    prior_tol=prior_tol,
                    r_tol=r_tol,
                    max_iter=max_iter,
                    null_weight=null_weight,
                    standardize=standardize,
                    coverage=coverage,
                    min_abs_corr=min_abs_corr,
                    tol=tol,
                    track_fit=False,
                    refine=False,
                )

                sinit2 = {key: s2[key] for key in ["alpha", "mu", "mu2"]}

                s3 = susie_suff_stat(
                    XtX=XtX,
                    Xty=Xty,
                    yty=yty,
                    n=n,
                    L=L,
                    X_colmeans=X_colmeans,
                    y_mean=y_mean,
                    prior_weights=pw_s,
                    s_init=sinit2,
                    scaled_prior_variance=scaled_prior_variance,
                    residual_variance=residual_variance,
                    estimate_residual_variance=estimate_residual_variance,
                    estimate_prior_variance=estimate_prior_variance,
                    estimate_prior_method=estimate_prior_method,
                    check_null_threshold=check_null_threshold,
                    prior_tol=prior_tol,
                    r_tol=r_tol,
                    max_iter=max_iter,
                    null_weight=null_weight,
                    standardize=standardize,
                    coverage=coverage,
                    min_abs_corr=min_abs_corr,
                    tol=tol,
                    track_fit=False,
                    refine=False,
                )

                m.append(s3)

            if len(m) == 0:
                conti = False
            else:
                elbo = [susie_get_objective(x) for x in m]
                if max(elbo) - susie_get_objective(s) <= 0:  # type: ignore
                    conti = False
                else:
                    s = m[np.argmax(elbo)]  # type: ignore

    return s


def susie_rss(
    z: Optional[np.ndarray] = None,
    R: Optional[np.ndarray] = None,
    n: Optional[int] = None,
    bhat: Optional[np.ndarray] = None,
    shat: Optional[np.ndarray] = None,
    var_y: Optional[float] = None,
    z_ld_weight: float = 0,
    estimate_residual_variance: bool = False,
    prior_variance: float = 50,
    check_prior: bool = True,
    **kwargs,
) -> Dict[str, Any]:
    """
    Sum of Single Effects (SuSiE) Regression using Summary Statistics.

    This function performs variable selection under a sparse Bayesian multiple linear regression
    of Y on X using the z-scores from standard univariate regression of Y on each column of X,
    an estimate, R, of the correlation matrix for the columns of X, and optionally (but strongly
    recommended) the sample size n.

    Parameters
    ----------
    z : np.ndarray, optional
        p-vector of z-scores.
    R : np.ndarray
        p x p correlation matrix.
    n : int, optional
        The sample size.
    bhat : np.ndarray, optional
        Alternative summary data giving the estimated effects (a vector of length p).
    shat : np.ndarray, optional
        Alternative summary data giving the standard errors of the estimated effects (a vector of length p).
    var_y : float, optional
        The sample variance of y, defined as y'y/(n-1).
    z_ld_weight : float, default=0
        Weight for adjusting the LD matrix (not recommended to use non-zero value).
    estimate_residual_variance : bool, default=False
        Whether to estimate residual variance.
    prior_variance : float, default=50
        The prior variance(s) for the non-zero noncentrality parameters.
    check_prior : bool, default=True
        Whether to check if the estimated prior variance becomes unreasonably large.
    **kwargs : dict
        Additional arguments to be passed to susie_suff_stat.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the SuSiE fit results, including:
        - alpha: An L by p matrix of posterior inclusion probabilities.
        - mu: An L by p matrix of posterior means, conditional on inclusion.
        - mu2: An L by p matrix of posterior second moments, conditional on inclusion.
        - lbf: log-Bayes Factor for each single effect.
        - lbf_variable: log-Bayes Factor for each variable and single effect.
        - V: Prior variance of the non-zero elements of b.
        - elbo: The value of the variational lower bound achieved at each iteration.
        - sets: Credible sets estimated from model fit.
        - pip: A vector of length p giving the posterior inclusion probabilities.
        - niter: Number of IBSS iterations performed.
        - converged: Whether the IBSS converged to a solution.

    Raises
    ------
    ValueError
        If inputs are invalid or inconsistent.

    Notes
    -----
    This function is a Python translation of the R function susie_rss.
    """
    if estimate_residual_variance:
        logger.warning(
            "For estimate_residual_variance = TRUE, please check "
            "that R is the 'in-sample' LD matrix; that is, the "
            "correlation matrix obtained using the exact same data "
            "matrix X that was used for the other summary "
            "statistics. Also note, when covariates are included in "
            "the univariate regressions that produced the summary "
            "statistics, also consider removing these effects from "
            "X before computing R."
        )

    # Check input R
    if z is None:
        p = len(bhat)  # type: ignore
    else:
        p = len(z)
    if R.shape[0] != p:  # type: ignore
        raise ValueError(
            f"The dimension of R ({R.shape[0]} x {R.shape[1]}) does not "  # type: ignore
            f"agree with expected ({p} x {p})"  # type: ignore
        )

    # Check input n
    if n is not None and n <= 1:
        raise ValueError("n must be greater than 1")

    # Check inputs z, bhat and shat
    if (z is None) == (bhat is None or shat is None):
        raise ValueError("Please provide either z or (bhat, shat), but not both")
    if z is None:
        if np.isscalar(shat):
            shat = np.full_like(bhat, shat)
        if len(bhat) != len(shat):  # type: ignore
            raise ValueError("The lengths of bhat and shat do not agree")
        if np.isnan(bhat).any() or np.isnan(shat).any():  # type: ignore
            raise ValueError("bhat, shat cannot have missing values")
        if (shat <= 0).any():  # type: ignore
            raise ValueError("shat cannot have zero or negative elements")
        z = bhat / shat  # type: ignore
    if len(z) < 1:  # type: ignore
        raise ValueError("Input vector z should have at least one element")
    z = np.nan_to_num(z)  # type: ignore

    # When n is provided, compute the PVE-adjusted z-scores
    if n is not None:
        adj = (n - 1) / (z**2 + n - 2)  # type: ignore
        z = np.sqrt(adj) * z

    # Modify R by z_ld_weight
    if z_ld_weight > 0:
        # warnings.warn("As of version 0.11.0, use of non-zero z_ld_weight is no longer recommended")
        R = muffled_cov2cor((1 - z_ld_weight) * R + z_ld_weight * np.outer(z, z))  # type: ignore
        R = (R + R.T) / 2  # type: ignore

    # Call susie_suff_stat
    if n is None:
        logger.warning(
            "Providing the sample size (n), or even a rough estimate of n, "
            "is highly recommended. Without n, the implicit assumption is "
            "n is large (Inf) and the effect sizes are small (close to zero)."
        )
        s = susie_suff_stat(
            XtX=R,  # type: ignore
            Xty=z,  # type: ignore
            n=2,
            yty=1,
            scaled_prior_variance=prior_variance,
            estimate_residual_variance=estimate_residual_variance,
            standardize=False,
            check_prior=check_prior,
            **kwargs,
        )
    else:
        if shat is not None and var_y is not None:
            XtXdiag = var_y * adj / (shat**2)
            XtX = (R * np.sqrt(XtXdiag)).T * np.sqrt(XtXdiag)
            XtX = (XtX + XtX.T) / 2
            Xty = z * np.sqrt(adj) * var_y / shat
        else:
            XtX = (n - 1) * R  # type: ignore
            Xty = np.sqrt(n - 1) * z
            var_y = 1
        s = susie_suff_stat(
            XtX=XtX,
            Xty=Xty,
            n=n,
            yty=(n - 1) * var_y,
            estimate_residual_variance=estimate_residual_variance,
            check_prior=check_prior,
            **kwargs,
        )

    return s


def summary_susie(object, **kwargs):
    """
    Summarize Susie results.

    Parameters
    ----------
    object : dict
        A susie fit object.
    **kwargs :
        Additional arguments passed to the generic summary or print.summary method.

    Returns
    -------
    dict
        A dictionary containing a DataFrame of variables and a DataFrame of credible sets.

    Raises
    ------
    ValueError
        If credible set information is not available in the susie object.
    """
    if "sets" not in object or object["sets"] is None:
        raise ValueError(
            "Cannot summarize SuSiE object because credible set information is not available"
        )

    variables = pd.DataFrame(
        {
            "variable": range(0, len(object["pip"])),
            "variable_prob": object["pip"],
            "cs": [-1] * len(object["pip"]),
        }
    )

    if object["null_index"] > 0:
        variables = variables.drop(object["null_index"] - 1)

    if "cs" in object["sets"] and object["sets"]["cs"] is not None:
        cs = pd.DataFrame(
            columns=["cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "variable"]
        )
        idx = 0
        for i, cs_set in object["sets"]["cs"].items():
            ith = int(i[1:])
            variables.loc[variables["variable"].isin(cs_set), "cs"] = ith

            new_row = {
                "cs": int(i[1:]),
                "cs_log10bf": object["lbf"][ith] / math.log(10),
                "cs_avg_r2": object["sets"]["purity"]["mean_abs_corr"][idx] ** 2,
                "cs_min_r2": object["sets"]["purity"]["min_abs_corr"][idx] ** 2,
                "variable": ",".join(map(str, cs_set)),
            }
            cs.loc[idx] = new_row  # type: ignore
            idx += 1

        variables = variables.sort_values("variable_prob", ascending=False)
    else:
        cs = None

    out = {"vars": variables, "cs": cs}
    return out


if __name__ == "__main__":
    np.random.seed(1)
    # read beta from betahat.txt
    beta = np.loadtxt("/project/voight_nextflow/wjh/WJH_packages/susiepy/betahat.txt")
    sebetahat = np.loadtxt(
        "/project/voight_nextflow/wjh/WJH_packages/susiepy/sebetahat.txt"
    )
    R = np.loadtxt("/project/voight_nextflow/wjh/WJH_packages/susiepy/R.txt")
    var_y = 7.84240878887082
    n = 574
    s = susie_rss(
        bhat=beta,
        shat=sebetahat,
        n=n,
        R=R,
        # var_y=var_y,
        L=10,
        # estimate_residual_variance=True,
    )
    res = summary_susie(s)
    res["vars"].to_csv("/project/voight_nextflow/wjh/WJH_packages/susiepy/pyvars.csv")
    res["cs"].to_csv("/project/voight_nextflow/wjh/WJH_packages/susiepy/pycs.csv")
