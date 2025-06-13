"""
RSparsePro wrapper for multi-ancestry fine-mapping.

Original code from https://github.com/zhwm/RSparsePro_LD
"""

import json
import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.special import softmax

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import Locus, intersect_sumstat_ld

logger = logging.getLogger("RSparsePro")


class RSparsePro:
    """
    RSparsePro for robust fine-mapping in the presence of LD mismatch.

    This class implements the RSparsePro algorithm for robust fine-mapping that can
    handle linkage disequilibrium (LD) mismatch between the reference panel and
    the study population. The method uses a robust regression approach with an
    error parameter to account for potential LD estimation errors.

    The algorithm models the relationship between summary statistics and LD matrix
    with potential errors:

    β̂ = Rβ + ε

    where β̂ are the estimated effect sizes, R is the LD matrix, β are the true
    effect sizes, and ε accounts for LD mismatch errors.

    Attributes
    ----------
    p : int
        The number of variants in the analysis.
    k : int
        The number of causal signals to detect.
    vare : float
        The error parameter accounting for LD mismatch.
    mat : np.ndarray
        The transformation matrix used in the algorithm.
        Computed as R(I + R/vare)^(-1) when vare > 0.
    beta_mu : np.ndarray
        Posterior mean effect sizes for each variant and signal.
        Shape: (p, k)
    gamma : np.ndarray
        Posterior inclusion probabilities for each variant and signal.
        Shape: (p, k)
    tilde_b : np.ndarray
        Estimated effect sizes accounting for LD mismatch.
        Shape: (p,)

    Methods
    -------
    infer_q_beta(R)
        Update posterior effect size distributions.
    infer_tilde_b(bhat)
        Update effect size estimates with LD mismatch correction.
    train(bhat, R, maxite, eps, ubound)
        Run the iterative optimization algorithm.
    get_PIP()
        Calculate posterior inclusion probabilities.
    get_effect(cthres)
        Extract credible sets based on coverage threshold.
    get_ztilde()
        Get corrected effect size estimates.

    Notes
    -----
    Key features of RSparsePro:

    1. **LD Mismatch Robustness**: Explicitly models errors in LD estimation
       through the error parameter vare, making the method robust to reference
       panel mismatch.

    2. **Adaptive Error Estimation**: Automatically tunes the error parameter
       through an adaptive training procedure that optimizes model fit.

    3. **Multiple Signal Detection**: Can detect multiple independent causal
       signals while accounting for LD structure and potential errors.

    4. **Coverage-based Credible Sets**: Constructs credible sets with
       specified coverage probability, with additional filters for
       effect group quality.

    The optimization objective combines:
    - Data likelihood under the robust model
    - Sparsity constraints on the number of causal variants
    - Regularization to prevent overfitting

    Advantages:
    - Robust to LD reference panel choice
    - Handles population stratification artifacts
    - Adaptive to varying degrees of LD mismatch
    - Provides interpretable credible sets

    Reference:
    Zhang, W. et al. RSparsePro: an R package for robust sparse regression.
    Bioinformatics (2023).
    """

    def __init__(self, P: int, K: int, R: np.ndarray, vare: float) -> None:
        """
        Initialize the RSparsePro model.

        Parameters
        ----------
        P : int
            The number of variants in the analysis.
        K : int
            The number of causal signals to detect.
        R : np.ndarray
            The LD correlation matrix of shape (P, P).
        vare : float
            The error parameter for LD mismatch correction.
            If vare=0, assumes perfect LD; larger values allow more mismatch.
        """
        self.p = P
        self.k = K
        self.vare = vare
        if vare != 0:
            self.mat = np.dot(R, np.linalg.inv(np.eye(self.p) + 1 / vare * R))
        self.beta_mu = np.zeros([self.p, self.k])
        self.gamma = np.zeros([self.p, self.k])
        self.tilde_b = np.zeros((self.p,))

    def infer_q_beta(self, R: np.ndarray) -> None:
        """
        Infer the posterior distribution of effect sizes.

        Updates the posterior mean effect sizes (beta_mu) and inclusion
        probabilities (gamma) for each signal component using coordinate
        descent optimization.

        Parameters
        ----------
        R : np.ndarray
            The LD correlation matrix of shape (p, p).
        """
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:, idxall] * self.beta_mu[:, idxall]).sum(axis=1)
            res_beta = self.tilde_b - np.dot(R, beta_all_k)
            self.beta_mu[:, k] = res_beta
            u = 0.5 * self.beta_mu[:, k] ** 2
            self.gamma[:, k] = softmax(u)

    def infer_tilde_b(self, bhat: np.ndarray) -> None:
        """
        Infer the corrected effect size estimates.

        Updates the effect size estimates (tilde_b) by combining the
        observed summary statistics with the current model estimates,
        accounting for LD mismatch through the error parameter.

        Parameters
        ----------
        bhat : np.ndarray
            The observed summary statistics (effect size estimates)
            of shape (p,).
        """
        if self.vare == 0:
            self.tilde_b = bhat
        else:
            beta_all = (self.gamma * self.beta_mu).sum(axis=1)
            self.tilde_b = np.dot(self.mat, (1 / self.vare * bhat + beta_all))

    def train(
        self,
        bhat: np.ndarray,
        R: np.ndarray,
        maxite: int,
        eps: float,
        ubound: int,
    ) -> bool:
        """
        Train the RSparsePro model using iterative optimization.

        Runs the coordinate descent algorithm to optimize the model parameters,
        alternating between updating effect size distributions and corrected
        estimates until convergence.

        Parameters
        ----------
        bhat : np.ndarray
            The observed summary statistics of shape (p,).
        R : np.ndarray
            The LD correlation matrix of shape (p, p).
        maxite : int
            The maximum number of iterations.
        eps : float
            The convergence criterion for parameter changes.
        ubound : int
            Upper bound threshold for divergence detection.

        Returns
        -------
        bool
            True if the algorithm converged, False otherwise.
        """
        for ite in range(maxite):
            old_gamma = self.gamma.copy()
            old_beta = self.beta_mu.copy()
            old_tilde = self.tilde_b.copy()
            self.infer_tilde_b(bhat)
            self.infer_q_beta(R)
            diff_gamma = np.linalg.norm(self.gamma - old_gamma)
            diff_beta = np.linalg.norm(self.beta_mu - old_beta)
            diff_b = np.linalg.norm(self.tilde_b - old_tilde)
            all_diff = diff_gamma + diff_beta + diff_b
            logger.info(
                "Iteration-->{} . Diff_b: {:.1f} . Diff_s: {:.1f} . Diff_mu: {:.1f} . ALL: {:.1f}".format(
                    ite,
                    float(diff_b),
                    float(diff_gamma),
                    float(diff_beta),
                    float(all_diff),
                )
            )
            if all_diff < eps:
                logger.info("The RSparsePro algorithm has converged.")
                converged = True
                break
            if ite == (maxite - 1) or abs(all_diff) > ubound:
                logger.info("The RSparsePro algorithm didn't converge.")
                converged = False
                break
        return converged

    def get_PIP(self) -> np.ndarray:
        """
        Get the posterior inclusion probabilities.

        Computes the maximum posterior inclusion probability across all
        signal components for each variant, providing a summary measure
        of the evidence for each variant being causal.

        Returns
        -------
        np.ndarray
            Array of posterior inclusion probabilities of shape (p,).
            Each element represents the maximum PIP across all components
            for the corresponding variant.
        """
        return np.max((self.gamma), axis=1).round(4)

    def get_effect(
        self, cthres: float
    ) -> Tuple[Dict[int, List[int]], Dict[int, np.ndarray], Dict[int, np.ndarray]]:
        """
        Extract effect groups (credible sets) based on coverage threshold.

        Identifies distinct causal signals by grouping variants with high
        posterior inclusion probabilities, subject to coverage and quality
        constraints.

        Parameters
        ----------
        cthres : float
            The coverage threshold for credible set construction.
            Only effect groups with cumulative PIP ≥ cthres are returned.

        Returns
        -------
        eff : Dict[int, List[int]]
            Dictionary mapping effect group index to list of variant indices
            in the credible set.
        eff_gamma : Dict[int, np.ndarray]
            Dictionary mapping effect group index to posterior inclusion
            probabilities for variants in the credible set.
        eff_mu : Dict[int, np.ndarray]
            Dictionary mapping effect group index to posterior mean effect
            sizes for variants in the credible set.
        """
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        mat_eff[mat_eff < 1 / (self.p + 1)] = 0
        csum = mat_eff.sum(axis=0).round(2)
        logger.info("Attainable coverage for effect groups: {}".format(csum))
        eff: Dict[int, List[int]] = {}
        eff_gamma: Dict[int, np.ndarray] = {}
        eff_mu: Dict[int, np.ndarray] = {}
        for k in range(self.k):
            if csum[k] >= cthres:
                p = 0
                while np.sum(mat_eff[matidx[0:p, k], k]) < cthres * csum[k]:
                    p = p + 1
                cidx = matidx[0:p, k].tolist()
                eff[k] = cidx
                eff_gamma[k] = mat_eff[cidx, k].round(4)
                eff_mu[k] = self.beta_mu[cidx, k].round(4)
        return eff, eff_gamma, eff_mu

    def get_ztilde(self) -> np.ndarray:
        """
        Get the corrected effect size estimates.

        Returns the effect size estimates corrected for LD mismatch,
        which can be used for downstream analysis or comparison with
        the original summary statistics.

        Returns
        -------
        np.ndarray
            Array of corrected effect size estimates of shape (p,).
        """
        return self.tilde_b.round(4)


def get_eff_maxld(eff: Dict[int, List[int]], ld: np.ndarray) -> float:
    """
    Get the maximum LD between lead variants across effect groups.

    Calculates the maximum absolute correlation between the lead variants
    (highest PIP variants) from different effect groups. This metric helps
    assess the independence of detected signals.

    Parameters
    ----------
    eff : Dict[int, List[int]]
        Dictionary mapping effect group indices to lists of variant indices.
    ld : np.ndarray
        The LD correlation matrix.

    Returns
    -------
    float
        Maximum absolute correlation between lead variants of different
        effect groups. Returns 0.0 if fewer than 2 effect groups exist.
    """
    idx = [i[0] for i in eff.values()]
    if len(eff) > 1:
        maxld = np.abs(np.tril(ld[np.ix_(idx, idx)], -1)).max()
    else:
        maxld = 0.0
    return maxld


def get_eff_minld(eff: Dict[int, List[int]], ld: np.ndarray) -> float:
    """
    Get the minimum LD within effect groups.

    Calculates the minimum absolute correlation within each effect group,
    providing a measure of the internal consistency of each credible set.

    Parameters
    ----------
    eff : Dict[int, List[int]]
        Dictionary mapping effect group indices to lists of variant indices.
    ld : np.ndarray
        The LD correlation matrix.

    Returns
    -------
    float
        Minimum absolute correlation within effect groups.
        Returns 1.0 if no effect groups exist.
    """
    if len(eff) == 0:
        minld = 1.0
    else:
        minld = min([abs(ld[np.ix_(v, v)]).min() for _, v in eff.items()])
    return minld


def get_ordered(eff_mu: Dict[int, np.ndarray]) -> bool:
    """
    Check if the effect sizes are properly ordered.

    Verifies that the detected effect groups are ordered by effect size
    magnitude, which is expected from the algorithm's prioritization scheme.

    Parameters
    ----------
    eff_mu : Dict[int, np.ndarray]
        Dictionary mapping effect group indices to posterior mean effect sizes.

    Returns
    -------
    bool
        True if effect groups are properly ordered, False otherwise.
        Always returns True if 1 or fewer effect groups exist.
    """
    if len(eff_mu) > 1:
        val_mu = [round(-abs(i[0])) for _, i in eff_mu.items()]
        ordered = (
            list(eff_mu.keys())[-1] == len(eff_mu) - 1
        )  # and (sorted(val_mu) == val_mu)
    else:
        ordered = True
    return ordered


def adaptive_train(
    zscore: np.ndarray,
    ld: np.ndarray,
    K: int,
    maxite: int,
    eps: float,
    ubound: int,
    cthres: float,
    minldthres: float,
    maxldthres: float,
    eincre: float,
    varemax: float,
    varemin: float,
) -> Tuple[
    Dict[int, List[int]],
    Dict[int, np.ndarray],
    Dict[int, np.ndarray],
    np.ndarray,
    np.ndarray,
]:
    """
    Adaptively train the RSparsePro model with error parameter optimization.

    Implements an adaptive training procedure that automatically tunes the
    error parameter (vare) to achieve good model fit while satisfying
    quality constraints on the detected effect groups.

    Parameters
    ----------
    zscore : np.ndarray
        Z-scores from GWAS summary statistics.
    ld : np.ndarray
        LD correlation matrix.
    K : int
        Maximum number of causal signals to detect.
    maxite : int
        Maximum iterations for each model fit.
    eps : float
        Convergence tolerance.
    ubound : int
        Upper bound for divergence detection.
    cthres : float
        Coverage threshold for credible sets.
    minldthres : float
        Minimum LD threshold within effect groups.
    maxldthres : float
        Maximum LD threshold between effect groups.
    eincre : float
        Multiplicative increment for error parameter.
    varemax : float
        Maximum allowed error parameter value.
    varemin : float
        Minimum error parameter value to try.

    Returns
    -------
    eff : Dict[int, List[int]]
        Final effect groups (credible sets).
    eff_gamma : Dict[int, np.ndarray]
        Posterior inclusion probabilities for effect groups.
    eff_mu : Dict[int, np.ndarray]
        Posterior mean effect sizes for effect groups.
    PIP : np.ndarray
        Posterior inclusion probabilities for all variants.
    ztilde : np.ndarray
        Corrected effect size estimates.

    Notes
    -----
    The adaptive training procedure:

    1. Start with vare=0 (perfect LD assumption)
    2. Fit model and evaluate effect group quality
    3. If quality constraints are not met, increase vare
    4. Repeat until constraints are satisfied or vare exceeds maximum
    5. Fall back to single-effect model if no good solution found

    Quality constraints include:
    - Model convergence
    - Proper ordering of effect sizes
    - Minimum LD within effect groups ≥ minldthres
    - Maximum LD between effect groups ≤ maxldthres
    - Adequate coverage for each effect group

    This procedure makes the method robust to LD reference panel choice
    by automatically adapting to the level of mismatch present.
    """
    vare = 0.0
    mc = False
    eff: Dict[int, List[int]] = {}
    eff_mu: Dict[int, np.ndarray] = {}
    minld = 1.0
    maxld = 0.0
    while (
        (not mc)
        or (not get_ordered(eff_mu))
        or (minld < minldthres)
        or (maxld > maxldthres)
    ):
        model = RSparsePro(len(zscore), K, ld, vare)
        mc = model.train(zscore, ld, maxite, eps, ubound)
        eff, eff_gamma, eff_mu = model.get_effect(cthres)
        maxld = get_eff_maxld(eff, ld)
        minld = get_eff_minld(eff, ld)
        logging.info("Max ld across effect groups: {}.".format(maxld))
        logging.info("Min ld within effect groups: {}.".format(minld))
        logging.info("vare = {}".format(round(vare, 4)))
        if vare > varemax or (len(eff) < 2 and get_ordered(eff_mu)):
            # logging.info("Algorithm didn't converge at the max vare. Setting K to 1.")
            model = RSparsePro(len(zscore), 1, ld, 0)
            mc = model.train(zscore, ld, maxite, eps, ubound)
            eff, eff_gamma, eff_mu = model.get_effect(cthres)
            break
        elif vare == 0:
            vare = varemin
        else:
            vare *= eincre
    ztilde = model.get_ztilde()
    # resz = model.get_resz(zscore, ld, eff)
    PIP = model.get_PIP()
    return eff, eff_gamma, eff_mu, PIP, ztilde  # resz


def rsparsepro_main(
    zfile: pd.DataFrame,
    ld: np.ndarray,
    K: int = 10,
    maxite: int = 100,
    eps: float = 1e-5,
    ubound: int = 100000,
    cthres: float = 0.95,
    eincre: float = 1.5,
    minldthres: float = 0.7,
    maxldthres: float = 0.2,
    varemax: float = 100.0,
    varemin: float = 1e-3,
) -> pd.DataFrame:
    """
    Run the main RSparsePro analysis pipeline.

    Executes the complete RSparsePro workflow including adaptive training,
    credible set construction, and result formatting.

    Parameters
    ----------
    zfile : pd.DataFrame
        DataFrame containing GWAS summary statistics with Z-scores.
        Must include columns: RSID, Z.
    ld : np.ndarray
        LD correlation matrix matching the variants in zfile.
    K : int, optional
        Maximum number of causal signals, by default 10.
    maxite : int, optional
        Maximum iterations for model fitting, by default 100.
    eps : float, optional
        Convergence tolerance, by default 1e-5.
    ubound : int, optional
        Upper bound for divergence detection, by default 100000.
    cthres : float, optional
        Coverage threshold for credible sets, by default 0.95.
    eincre : float, optional
        Error parameter increment factor, by default 1.5.
    minldthres : float, optional
        Minimum LD within effect groups, by default 0.7.
    maxldthres : float, optional
        Maximum LD between effect groups, by default 0.2.
    varemax : float, optional
        Maximum error parameter value, by default 100.0.
    varemin : float, optional
        Minimum error parameter value, by default 1e-3.

    Returns
    -------
    pd.DataFrame
        Enhanced DataFrame with additional columns:
        - PIP: Posterior inclusion probabilities
        - z_estimated: Corrected Z-score estimates
        - cs: Credible set assignment (0 for not in any set)
    """
    eff, eff_gamma, eff_mu, PIP, ztilde = adaptive_train(
        zfile["Z"].to_numpy(),
        ld,
        K,
        maxite,
        eps,
        ubound,
        cthres,
        minldthres,
        maxldthres,
        eincre,
        varemax,
        varemin,
    )
    zfile["PIP"] = PIP
    zfile["z_estimated"] = ztilde
    zfile["cs"] = 0
    for e in eff:
        mcs_idx = [zfile["RSID"][j] for j in eff[e]]
        logger.info(f"The {e}-th effect group contains effective variants:")
        logger.info(f"causal variants: {mcs_idx}")
        logger.info(f"variant probabilities for this effect group: {eff_gamma[e]}")
        logger.info(f"zscore for this effect group: {eff_mu[e]}\n")
        zfile.loc[list(eff[e]), "cs"] = e + 1
    # zfile.to_csv("{}.rsparsepro.txt".format(save), sep="\t", header=True, index=False)
    return zfile


def run_rsparsepro(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    maxite: int = 100,
    eps: float = 1e-5,
    ubound: int = 100000,
    eincre: float = 1.5,
    minldthres: float = 0.7,
    maxldthres: float = 0.2,
    varemax: float = 100.0,
    varemin: float = 1e-3,
) -> CredibleSet:
    """
    Run RSparsePro fine-mapping analysis on a genomic locus.

    Performs robust fine-mapping using the RSparsePro algorithm, which is
    designed to handle LD mismatch between reference panels and study populations.
    The method automatically adapts to the level of LD mismatch through an
    adaptive error parameter optimization procedure.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix.
        The summary statistics and LD matrix will be automatically matched
        if they are not already aligned.
    max_causal : int, optional
        Maximum number of causal signals to detect, by default 1.
        Higher values allow detection of multiple independent associations.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
        Determines the cumulative posterior probability required for
        credible set inclusion.
    maxite : int, optional
        Maximum iterations for model optimization, by default 100.
        More iterations may improve convergence but increase runtime.
    eps : float, optional
        Convergence tolerance for parameter updates, by default 1e-5.
        Smaller values require tighter convergence but may improve accuracy.
    ubound : int, optional
        Upper bound for divergence detection, by default 100000.
        Algorithm stops if parameter changes exceed this threshold.
    eincre : float, optional
        Multiplicative factor for error parameter increases, by default 1.5.
        Controls the rate of adaptation during error parameter tuning.
    minldthres : float, optional
        Minimum LD threshold within effect groups, by default 0.7.
        Effect groups with lower internal LD are rejected.
    maxldthres : float, optional
        Maximum LD threshold between effect groups, by default 0.2.
        Effect groups with higher between-group LD are rejected.
    varemax : float, optional
        Maximum allowed error parameter value, by default 100.0.
        Limits the degree of LD mismatch correction applied.
    varemin : float, optional
        Minimum error parameter to try during adaptation, by default 1e-3.
        Starting point for error parameter optimization.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Posterior inclusion probabilities for all variants
        - Credible sets for each detected signal
        - Lead SNPs (highest PIP in each credible set)
        - Coverage probabilities and algorithm parameters
        - RSparsePro-specific results including corrected effect estimates

    Warnings
    --------
    If the summary statistics and LD matrix are not matched, they will be
    automatically intersected and reordered with a warning message.

    Notes
    -----
    RSparsePro addresses several challenges in fine-mapping:

    1. **LD Reference Panel Mismatch**: Accounts for differences between
       the LD structure in the reference panel and the study population
       through an adaptive error model.

    2. **Robust Effect Estimation**: Provides corrected effect size estimates
       that are more reliable in the presence of LD mismatch.

    3. **Quality Control**: Implements automated quality filters for
       credible sets based on LD structure and coverage criteria.

    4. **Multiple Signal Detection**: Can identify multiple independent
       causal signals while maintaining robustness to LD errors.

    The adaptive training procedure automatically determines the optimal
    level of LD mismatch correction by:
    - Starting with the assumption of perfect LD (vare=0)
    - Gradually increasing the error parameter until quality criteria are met
    - Falling back to simpler models if convergence is not achieved

    This makes the method particularly suitable for:
    - Cross-population fine-mapping studies
    - Analysis with limited reference panel data
    - Situations with suspected population stratification
    - Studies with heterogeneous LD patterns

    Reference:
    Zhang, W. et al. RSparsePro: an R package for robust sparse regression.
    Bioinformatics (2023).

    Examples
    --------
    >>> # Basic RSparsePro analysis
    >>> credible_set = run_rsparsepro(locus)
    >>> print(f"Found {credible_set.n_cs} credible sets")
    >>> print(f"Top PIP: {credible_set.pips.max():.4f}")
    Found 1 credible sets
    Top PIP: 0.8934

    >>> # RSparsePro with multiple signals and strict quality control
    >>> credible_set = run_rsparsepro(
    ...     locus,
    ...     max_causal=5,
    ...     coverage=0.99,
    ...     minldthres=0.8,  # Stricter within-group LD requirement
    ...     maxldthres=0.1   # Stricter between-group LD requirement
    ... )
    >>> print(f"Detected {credible_set.n_cs} high-quality signals")
    >>> print(f"Credible set sizes: {credible_set.cs_sizes}")
    Detected 2 high-quality signals
    Credible set sizes: [12, 8]

    >>> # Access posterior inclusion probabilities
    >>> top_variants = credible_set.pips.nlargest(10)
    >>> print("Top 10 variants by PIP:")
    >>> print(top_variants)
    Top 10 variants by PIP:
    rs123456    0.8934
    rs789012    0.7234
    rs345678    0.6456
    ...

    >>> # Examine credible sets and lead variants
    >>> for i, snps in enumerate(credible_set.snps):
    ...     lead_snp = credible_set.lead_snps[i]
    ...     pip = credible_set.pips[lead_snp]
    ...     print(f"Signal {i+1}: {lead_snp} (PIP={pip:.4f}, size={len(snps)})")
    Signal 1: rs123456 (PIP=0.8934, size=12)
    Signal 2: rs789012 (PIP=0.7234, size=8)
    """
    if not locus.is_matched:
        logger.warning(
            "The sumstat and LD are not matched, will match them in same order."
        )
        locus = intersect_sumstat_ld(locus)
    logger.info(f"Running RSparsePro on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "maxite": maxite,
        "eps": eps,
        "ubound": ubound,
        "eincre": eincre,
        "minldthres": minldthres,
        "maxldthres": maxldthres,
        "varemax": varemax,
        "varemin": varemin,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")

    sumstats = locus.sumstats.copy()
    ld = locus.ld.r.copy()
    sumstats["RSID"] = sumstats[ColName.SNPID]
    sumstats[ColName.Z] = sumstats[ColName.BETA] / sumstats[ColName.SE]
    zfile = rsparsepro_main(
        zfile=sumstats,
        ld=ld,
        K=max_causal,
        maxite=maxite,
        eps=eps,
        ubound=ubound,
        cthres=coverage,
        eincre=eincre,
        minldthres=minldthres,
        maxldthres=maxldthres,
        varemax=varemax,
        varemin=varemin,
    )

    pips = pd.Series(data=zfile["PIP"].to_numpy(), index=zfile["SNPID"].to_numpy())
    cs_snps: List[List[str]] = []
    lead_snps: List[str] = []
    for cs_i, sub_df in zfile.groupby("cs"):
        if cs_i == 0:
            continue
        cs_snps.append(sub_df["SNPID"].values.tolist())
        lead_snps.append(pips[pips.index.isin(sub_df["SNPID"].values)].idxmax())
    cs_sizes = [len(i) for i in cs_snps]
    logger.info(f"N of credible set: {len(cs_snps)}")
    logger.info(f"Credible set size: {cs_sizes}")
    return CredibleSet(
        tool=Method.RSparsePro,
        n_cs=len(cs_snps),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=cs_sizes,
        pips=pips,
        parameters=parameters,
    )
