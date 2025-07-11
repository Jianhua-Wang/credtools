"""Quality control functions for CREDTOOLS data."""

import logging
import os
from multiprocessing import Pool
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeRemainingColumn,
)
from scipy import stats
from scipy.optimize import curve_fit, minimize_scalar
from sklearn.mixture import GaussianMixture

from credtools.constants import ColName
from credtools.locus import Locus, LocusSet, intersect_sumstat_ld, load_locus_set

logger = logging.getLogger("QC")


def get_eigen(ldmatrix: np.ndarray) -> Dict[str, np.ndarray]:
    """
    Compute eigenvalues and eigenvectors of LD matrix.

    Parameters
    ----------
    ldmatrix : np.ndarray
        A p by p symmetric, positive semidefinite correlation matrix.

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary containing eigenvalues and eigenvectors with keys:
        - 'eigvals': eigenvalues array
        - 'eigvecs': eigenvectors matrix

    Notes
    -----
    TODO: accelerate with joblib for large matrices.

    This function uses numpy.linalg.eigh which is optimized for symmetric matrices
    and returns eigenvalues in ascending order.
    """
    # ldmatrix = ldmatrix.astype(np.float32)
    eigvals, eigvecs = np.linalg.eigh(ldmatrix)
    return {"eigvals": eigvals, "eigvecs": eigvecs}


def estimate_s_rss(
    locus: Locus,
    r_tol: float = 1e-8,
    method: str = "null-mle",
    eigvens: Optional[Dict[str, np.ndarray]] = None,
) -> float:
    """
    Estimate s parameter in the susie_rss Model Using Regularized LD.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix.
    r_tol : float, optional
        Tolerance level for eigenvalue check of positive semidefinite matrix, by default 1e-8.
    method : str, optional
        Method to estimate s, by default "null-mle".
        Options: "null-mle", "null-partialmle", or "null-pseudomle".
    eigvens : Optional[Dict[str, np.ndarray]], optional
        Pre-computed eigenvalues and eigenvectors, by default None.

    Returns
    -------
    float
        Estimated s value between 0 and 1 (or potentially > 1 for "null-partialmle").

    Raises
    ------
    ValueError
        If n <= 1 or if the method is not implemented.

    Notes
    -----
    This function estimates the parameter s, which provides information about the
    consistency between z-scores and the LD matrix. A larger s indicates a strong
    inconsistency between z-scores and the LD matrix.

    The function implements three estimation methods:

    - "null-mle": Maximum likelihood estimation under the null model
    - "null-partialmle": Partial MLE using null space projection
    - "null-pseudomle": Pseudo-likelihood estimation

    The z-scores are transformed using the formula:
    z_transformed = sqrt(sigma2) * z
    where sigma2 = (n-1) / (z^2 + n-2)
    """
    # make sure the LD matrix and sumstats file are matched
    input_locus = locus.copy()
    input_locus = intersect_sumstat_ld(input_locus)
    z = (
        input_locus.sumstats[ColName.BETA] / input_locus.sumstats[ColName.SE]
    ).to_numpy()
    n = input_locus.sample_size
    # Check and process input arguments z, R
    z = np.where(np.isnan(z), 0, z)
    if eigvens is not None:
        eigvals = eigvens["eigvals"]
        eigvecs = eigvens["eigvecs"]
    else:
        eigens = get_eigen(input_locus.ld.r)
        eigvals = eigens["eigvals"]
        eigvecs = eigens["eigvecs"]

    # if np.any(eigvals < -r_tol):
    #     logger.warning("The LD matrix is not positive semidefinite. Negative eigenvalues are set to zero")
    eigvals[eigvals < r_tol] = 0

    if n <= 1:
        raise ValueError("n must be greater than 1")

    sigma2 = (n - 1) / (z**2 + n - 2)
    z = np.sqrt(sigma2) * z

    if method == "null-mle":

        def negloglikelihood(s, ztv, d):
            denom = (1 - s) * d + s
            term1 = 0.5 * np.sum(np.log(denom))
            term2 = 0.5 * np.sum((ztv / denom) * ztv)
            return term1 + term2

        ztv = eigvecs.T @ z
        result = minimize_scalar(
            negloglikelihood,
            bounds=(0, 1),
            method="bounded",
            args=(ztv, eigvals),
            options={"xatol": np.sqrt(np.finfo(float).eps)},
        )
        s = result.x  # type: ignore

    elif method == "null-partialmle":
        colspace = np.where(eigvals > 0)[0]
        if len(colspace) == len(z):
            s = 0
        else:
            znull = eigvecs[:, ~np.isin(np.arange(len(z)), colspace)].T @ z
            s = np.sum(znull**2) / len(znull)

    elif method == "null-pseudomle":

        def pseudolikelihood(
            s: float, z: np.ndarray, eigvals: np.ndarray, eigvecs: np.ndarray
        ) -> float:
            precision = eigvecs @ (eigvecs.T / ((1 - s) * eigvals + s))
            postmean = np.zeros_like(z)
            postvar = np.zeros_like(z)
            for i in range(len(z)):
                postmean[i] = -(1 / precision[i, i]) * precision[i, :].dot(z) + z[i]
                postvar[i] = 1 / precision[i, i]
            return -np.sum(stats.norm.logpdf(z, loc=postmean, scale=np.sqrt(postvar)))

        result = minimize_scalar(
            pseudolikelihood,
            bounds=(0, 1),
            method="bounded",
            args=(z, eigvals, eigvecs),
        )
        s = result.x  # type: ignore

    else:
        raise ValueError("The method is not implemented")

    return s  # type: ignore


def kriging_rss(
    locus: Locus,
    r_tol: float = 1e-8,
    s: Optional[float] = None,
    eigvens: Optional[Dict[str, np.ndarray]] = None,
) -> pd.DataFrame:
    """
    Compute distribution of z-scores of variant j given other z-scores, and detect possible allele switch issues.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix.
    r_tol : float, optional
        Tolerance level for eigenvalue check of positive semidefinite matrix, by default 1e-8.
    s : Optional[float], optional
        An estimated s parameter from estimate_s_rss function, by default None.
        If None, s will be estimated automatically.
    eigvens : Optional[Dict[str, np.ndarray]], optional
        Pre-computed eigenvalues and eigenvectors, by default None.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the results of the kriging RSS test with columns:
        - SNPID: SNP identifier
        - z: transformed z-score
        - condmean: conditional mean
        - condvar: conditional variance
        - z_std_diff: standardized difference
        - logLR: log likelihood ratio

    Raises
    ------
    ValueError
        If n <= 1.

    Notes
    -----
    Under the null hypothesis, the RSS model with regularized LD matrix assumes:
    z|R,s ~ N(0, (1-s)R + sI)

    This function uses a mixture of normals to model the conditional distribution
    of z_j given other z-scores. The method can help detect:

    - Allele switch issues
    - Outlier variants
    - LD inconsistencies

    The algorithm:
    1. Computes conditional means and variances for each variant
    2. Fits a Gaussian mixture model to capture heterogeneity
    3. Calculates likelihood ratios for allele switch detection
    """
    # Check and process input arguments z, R
    input_locus = locus.copy()
    input_locus = intersect_sumstat_ld(input_locus)
    z = (
        input_locus.sumstats[ColName.BETA] / input_locus.sumstats[ColName.SE]
    ).to_numpy()
    n = input_locus.sample_size
    z = np.where(np.isnan(z), 0, z)

    # Compute eigenvalues and eigenvectors
    if eigvens is not None:
        eigvals = eigvens["eigvals"]
        eigvecs = eigvens["eigvecs"]
    else:
        eigens = get_eigen(input_locus.ld.r)
        eigvals = eigens["eigvals"]
        eigvecs = eigens["eigvecs"]
    if s is None:
        s = estimate_s_rss(locus, eigvens={"eigvals": eigvals, "eigvecs": eigvecs})
    eigvals = eigvals[::-1]
    eigvecs = eigvecs[:, ::-1]

    eigvals[eigvals < r_tol] = 0

    if n <= 1:
        raise ValueError("n must be greater than 1")

    sigma2 = (n - 1) / (z**2 + n - 2)
    z = np.sqrt(sigma2) * z

    dinv = 1 / ((1 - s) * eigvals + s)
    dinv[np.isinf(dinv)] = 0
    precision = eigvecs @ (eigvecs * dinv).T
    condmean = np.zeros_like(z)
    condvar = np.zeros_like(z)
    for i in range(len(z)):
        condmean[i] = -(1 / precision[i, i]) * precision[i, :i].dot(z[:i]) - (
            1 / precision[i, i]
        ) * precision[i, i + 1 :].dot(z[i + 1 :])
        condvar[i] = 1 / precision[i, i]
    z_std_diff = (z - condmean) / np.sqrt(condvar)

    # Obtain grid
    a_min = 0.8
    a_max = 2 if np.max(z_std_diff**2) < 1 else 2 * np.sqrt(np.max(z_std_diff**2))
    npoint = int(np.ceil(np.log2(a_max / a_min) / np.log2(1.05)))
    # Ensure npoint doesn't exceed number of samples
    npoint = min(npoint, len(z) - 1)
    a_grid = 1.05 ** np.arange(-npoint, 1) * a_max

    # Compute likelihood
    sd_mtx = np.outer(np.sqrt(condvar), a_grid)
    matrix_llik = stats.norm.logpdf(
        z[:, np.newaxis] - condmean[:, np.newaxis], scale=sd_mtx
    )
    lfactors = np.max(matrix_llik, axis=1)
    matrix_llik = matrix_llik - lfactors[:, np.newaxis]

    # Estimate weight using Gaussian Mixture Model
    gmm = GaussianMixture(
        n_components=len(a_grid), covariance_type="diag", max_iter=1000
    )
    gmm.fit(matrix_llik)
    w = gmm.weights_

    # Compute denominators in likelihood ratios
    logl0mix = np.log(np.sum(np.exp(matrix_llik) * (w + 1e-15), axis=1)) + lfactors  # type: ignore

    # Compute numerators in likelihood ratios
    matrix_llik = stats.norm.logpdf(
        z[:, np.newaxis] + condmean[:, np.newaxis], scale=sd_mtx
    )
    lfactors = np.max(matrix_llik, axis=1)
    matrix_llik = matrix_llik - lfactors[:, np.newaxis]
    logl1mix = np.log(np.sum(np.exp(matrix_llik) * (w + 1e-15), axis=1)) + lfactors  # type: ignore

    # Compute (log) likelihood ratios
    logLRmix = logl1mix - logl0mix

    res = pd.DataFrame(
        {
            "SNPID": input_locus.sumstats[ColName.SNPID].to_numpy(),
            "z": z,
            "condmean": condmean,
            "condvar": condvar,
            "z_std_diff": z_std_diff,
            "logLR": logLRmix,
        },
        # index=input_locus.sumstats[ColName.SNPID].to_numpy(),
    )
    # TODO: remove variants with logLR > 2 and abs(z) > 2

    return res


def compute_dentist_s(locus: Locus) -> pd.DataFrame:
    """
    Compute Dentist-S statistic and p-value for outlier detection.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the results of the Dentist-S test with columns:
        - SNPID: SNP identifier
        - t_dentist_s: Dentist-S test statistic
        - p_dentist_s: log p-value

    Notes
    -----
    Reference: https://github.com/mkanai/slalom/blob/854976f8e19e6fad2db3123eb9249e07ba0e1c1b/slalom.py#L254

    The Dentist-S statistic tests for outliers by comparing each variant's z-score
    to what would be expected based on its LD with the lead variant:

    t_dentist_s = (z_j - r_jk * z_k)^2 / (1 - r_jk^2)

    where:
    - z_j: z-score for variant j
    - z_k: z-score for lead variant k
    - r_jk: LD correlation between variants j and k

    TODO: Use ABF to select lead variant, although in most cases the lead variant
    is the one with the smallest p-value.
    """
    input_locus = locus.copy()
    input_locus = intersect_sumstat_ld(input_locus)
    df = input_locus.sumstats.copy()
    df["Z"] = df[ColName.BETA] / df[ColName.SE]
    lead_idx = df[ColName.P].idxmin()
    # TODO: use abf to select lead variant, although in most cases the lead variant is the one with the smallest p-value
    lead_z = df.loc[lead_idx, "Z"]
    df["r"] = input_locus.ld.r[lead_idx]

    df["t_dentist_s"] = (df["Z"] - df["r"] * lead_z) ** 2 / (1 - df["r"] ** 2)  # type: ignore
    df["t_dentist_s"] = np.where(df["t_dentist_s"] < 0, np.inf, df["t_dentist_s"])
    df.at[lead_idx, "t_dentist_s"] = np.nan
    df["p_dentist_s"] = stats.chi2.logsf(df["t_dentist_s"], df=1)

    df = df[[ColName.SNPID, "t_dentist_s", "p_dentist_s"]].copy()
    # df.set_index(ColName.SNPID, inplace=True)
    # df.index.name = None
    return df


def compare_maf(locus: Locus) -> pd.DataFrame:
    """
    Compare allele frequencies between summary statistics and LD reference.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the comparison results with columns:
        - SNPID: SNP identifier
        - MAF_sumstats: MAF from summary statistics
        - MAF_ld: MAF from LD reference

        Returns empty DataFrame if AF2 column is not available in LD matrix.

    Warnings
    --------
    If AF2 column is not present in the LD matrix, a warning is logged.

    Notes
    -----
    This function compares minor allele frequencies (MAF) between:

    1. Summary statistics (derived from EAF)
    2. LD reference panel (from AF2 column)

    Large discrepancies may indicate:
    - Population stratification
    - Allele frequency differences between studies
    - Potential data quality issues

    MAF is calculated as min(AF, 1-AF) for both sources.
    """
    input_locus = locus.copy()
    if "AF2" not in input_locus.ld.map.columns:
        logger.warning("AF2 is not in the LD matrix.")
        return pd.DataFrame()
    input_locus = intersect_sumstat_ld(input_locus)
    df = input_locus.sumstats[[ColName.SNPID, ColName.MAF]].copy()
    df.rename(columns={ColName.MAF: "MAF_sumstats"}, inplace=True)
    df.set_index(ColName.SNPID, inplace=True)
    af_ld = pd.Series(
        index=input_locus.ld.map[ColName.SNPID].tolist(),
        data=input_locus.ld.map["AF2"].values,
    )
    maf_ld = np.minimum(af_ld, 1 - af_ld)
    df["MAF_ld"] = maf_ld
    df[ColName.SNPID] = df.index
    df.reset_index(drop=True, inplace=True)
    return df


def snp_missingness(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute the missingness rate of each cohort across variants.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple loci/cohorts.

    Returns
    -------
    pd.DataFrame
        DataFrame with variants as rows and cohorts as columns, where 1 indicates
        presence and 0 indicates absence of the variant in that cohort.

    Warnings
    --------
    If any cohort has a missing rate > 0.1, a warning is logged.

    Notes
    -----
    This function:

    1. Identifies all unique variants across cohorts
    2. Creates a binary matrix indicating variant presence/absence
    3. Calculates and logs missing rates for each cohort
    4. Issues warnings for cohorts with high missing rates (>10%)

    High missing rates may indicate:
    - Different genotyping platforms
    - Different quality control criteria
    - Population-specific variants
    """
    missingness_df = []
    for locus in locus_set.loci:
        loc = intersect_sumstat_ld(locus)
        loc = loc.sumstats[[ColName.SNPID]].copy()
        loc[f"{locus.popu}_{locus.cohort}"] = 1
        loc.set_index(ColName.SNPID, inplace=True)
        missingness_df.append(loc)
    missingness_df = pd.concat(missingness_df, axis=1)
    missingness_df.fillna(0, inplace=True)
    # log warning if missing rate > 0.1
    for col in missingness_df.columns:
        missing_rate = float(
            round(1 - missingness_df[col].sum() / missingness_df.shape[0], 3)
        )
        if missing_rate > 0.1:
            logger.warning(f"The missing rate of {col} is {missing_rate}")
        else:
            logger.info(f"The missing rate of {col} is {missing_rate}")

    return missingness_df


def ld_4th_moment(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute the 4th moment of the LD matrix as a measure of LD structure.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple loci/cohorts.

    Returns
    -------
    pd.DataFrame
        DataFrame with variants as rows and cohorts as columns, containing
        the 4th moment values (sum of r^4 - 1 for each variant).

    Notes
    -----
    The 4th moment is calculated as:
    Σ(r_ij^4) - 1 for each variant i

    This metric provides information about:
    - LD structure complexity
    - Potential issues with LD matrix quality
    - Differences in LD patterns between populations

    Higher values may indicate:
    - Strong local LD structure
    - Potential genotyping errors
    - Population stratification effects

    The function intersects variants across all cohorts to ensure fair comparison.
    """
    ld_4th_res = []
    # intersect between loci
    overlap_snps = set(locus_set.loci[0].sumstats[ColName.SNPID])
    for locus in locus_set.loci[1:]:
        overlap_snps = overlap_snps.intersection(set(locus.sumstats[ColName.SNPID]))
    for locus in locus_set.loci:
        locus = locus.copy()
        locus.sumstats = locus.sumstats[
            locus.sumstats[ColName.SNPID].isin(overlap_snps)
        ]
        locus = intersect_sumstat_ld(locus)
        r_4th = pd.Series(
            index=locus.ld.map[ColName.SNPID], data=np.power(locus.ld.r, 4).sum(axis=0)
        )
        r_4th = r_4th - 1
        r_4th.name = f"{locus.popu}_{locus.cohort}"
        ld_4th_res.append(r_4th)
    return pd.concat(ld_4th_res, axis=1)


def ld_decay(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute LD decay patterns across cohorts.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple loci/cohorts.

    Returns
    -------
    pd.DataFrame
        DataFrame containing LD decay information with columns:
        - distance_kb: distance in kilobases
        - r2_avg: average r² value at that distance
        - decay_rate: fitted exponential decay rate parameter
        - cohort: cohort identifier

    Notes
    -----
    This function analyzes LD decay by:

    1. Computing pairwise distances between all variants
    2. Binning distances into 1kb windows
    3. Calculating average r² within each distance bin
    4. Fitting an exponential decay model: r² = a * exp(-b * distance)

    LD decay patterns can reveal:
    - Population-specific recombination patterns
    - Effective population size differences
    - Demographic history effects

    Different populations typically show different decay rates due to:
    - Historical effective population sizes
    - Admixture patterns
    - Founder effects
    """

    def fit_exp(x: np.ndarray, a: float, b: float) -> np.ndarray:
        with np.errstate(over="ignore"):
            return a * np.exp(-b * x)

    binsize = 1000
    decay_res = []
    for locus in locus_set.loci:
        ldmap = locus.ld.map.copy()
        r = locus.ld.r.copy()
        distance_mat = np.array(
            [ldmap["BP"] - ldmap["BP"].values[i] for i in range(len(ldmap))]
        )
        distance_mat = distance_mat[np.tril_indices_from(distance_mat, k=-1)].flatten()
        distance_mat = np.abs(distance_mat)
        r = r[np.tril_indices_from(r, k=-1)].flatten()
        r = np.square(r)
        bins = np.arange(0, ldmap["BP"].max() - ldmap["BP"].min() + binsize, binsize)

        r_sum, _ = np.histogram(distance_mat, bins=bins, weights=r)
        count, _ = np.histogram(distance_mat, bins=bins)

        with np.errstate(divide="ignore", invalid="ignore"):
            r2_avg = np.where(count > 0, r_sum / count, 0)
        popt, _ = curve_fit(fit_exp, bins[1:] / binsize, r2_avg)
        res = pd.DataFrame(
            {
                "distance_kb": bins[1:] / binsize,
                "r2_avg": r2_avg,
                "decay_rate": popt[0],
                "cohort": f"{locus.popu}_{locus.cohort}",
            }
        )
        decay_res.append(res)
    return pd.concat(decay_res, axis=0)


def cochran_q(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute Cochran-Q statistic for heterogeneity testing across cohorts.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple loci/cohorts.

    Returns
    -------
    pd.DataFrame
        DataFrame with SNPID as index and columns:
        - Q: Cochran-Q test statistic
        - Q_pvalue: p-value from chi-squared test
        - I_squared: I² heterogeneity statistic (percentage)

    Notes
    -----
    The Cochran-Q test assesses heterogeneity in effect sizes across studies:

    Q = Σ w_i(β_i - β_pooled)²

    where:
    - w_i = 1/SE_i² (inverse variance weights)
    - β_i = effect size in study i
    - β_pooled = weighted average effect size

    The I² statistic quantifies the proportion of total variation due to
    heterogeneity rather than chance:

    I² = max(0, (Q - df)/Q × 100%)

    Interpretation:
    - Q p-value < 0.05: significant heterogeneity
    - I² > 50%: substantial heterogeneity
    - I² > 75%: considerable heterogeneity

    High heterogeneity may indicate:
    - Population differences
    - Different LD patterns
    - Batch effects
    - Population stratification
    """
    merged_df = locus_set.loci[0].original_sumstats[[ColName.SNPID]].copy()
    for i, locus_obj in enumerate(locus_set.loci):
        locus_df = locus_obj.sumstats[
            [ColName.SNPID, ColName.BETA, ColName.SE, ColName.EAF]
        ].copy()
        locus_df.rename(
            columns={
                ColName.BETA: f"BETA_{i}",
                ColName.SE: f"SE_{i}",
                ColName.EAF: f"EAF_{i}",
            },
            inplace=True,
        )
        merged_df = pd.merge(
            merged_df, locus_df, on=ColName.SNPID, how="inner", suffixes=("", f"_{i}")
        )

    k: int = len(locus_set.loci)
    weights = []
    effects = []
    for i in range(k):
        weights.append((1 / (merged_df[f"SE_{i}"] ** 2)))
        effects.append(merged_df[f"BETA_{i}"])

    # Calculate weighted mean effect size
    weighted_mean = np.sum([w * e for w, e in zip(weights, effects)], axis=0) / np.sum(
        weights, axis=0
    )

    # Calculate Q statistic
    Q = np.sum([w * (e - weighted_mean) ** 2 for w, e in zip(weights, effects)], axis=0)

    # Calculate degrees of freedom
    df = k - 1

    # Calculate P-value
    p_value = stats.chi2.sf(Q, df)

    # Calculate I^2
    with np.errstate(invalid="ignore"):
        I_squared = np.maximum(0, (Q - df) / Q * 100)

    # Create output dataframe
    output_df = pd.DataFrame(
        {
            "SNPID": merged_df["SNPID"],
            "Q": Q,
            "Q_pvalue": p_value,
            "I_squared": I_squared,
        }
    )
    return output_df.set_index(ColName.SNPID)


def locus_qc(
    locus_set: LocusSet,
    r_tol: float = 1e-3,
    method: str = "null-mle",
    out_dir: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Perform comprehensive quality control analysis for a locus.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing loci to analyze.
    r_tol : float, optional
        Tolerance level for eigenvalue check of positive semidefinite matrix, by default 1e-3.
    method : str, optional
        Method to estimate s parameter, by default "null-mle".
        Options: "null-mle", "null-partialmle", or "null-pseudomle".
    out_dir : Optional[str], optional
        Output directory to save results, by default None.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary of quality control results with keys:
        - 'expected_z': kriging RSS results
        - 'dentist_s': Dentist-S test results
        - 'compare_maf': MAF comparison results
        - 'ld_4th_moment': 4th moment of LD matrix
        - 'ld_decay': LD decay analysis
        - 'cochran_q': heterogeneity test (if multiple cohorts)
        - 'snp_missingness': missingness analysis (if multiple cohorts)

    Notes
    -----
    This function performs comprehensive QC including:

    Single-locus analyses:
    - Kriging RSS for outlier detection
    - Dentist-S for outlier detection
    - MAF comparison between sumstats and LD reference
    - LD matrix 4th moment analysis
    - LD decay pattern analysis

    Multi-locus analyses (when applicable):
    - Cochran-Q heterogeneity testing
    - SNP missingness across cohorts

    TODO: Add LAVA (Local Analysis of Variant Associations) analysis.

    If out_dir is provided, results are saved as tab-separated files.
    """
    qc_metrics = {}
    all_expected_z = []
    all_dentist_s = []
    all_compare_maf = []
    for locus in locus_set.loci:
        lo = intersect_sumstat_ld(locus)
        eigens = get_eigen(lo.ld.r)
        lambda_s = estimate_s_rss(locus, r_tol, method, eigens)
        expected_z = kriging_rss(locus, r_tol, lambda_s, eigens)
        expected_z["lambda_s"] = lambda_s
        expected_z["cohort"] = f"{locus.popu}_{locus.cohort}"
        dentist_s = compute_dentist_s(locus)
        dentist_s["cohort"] = f"{locus.popu}_{locus.cohort}"
        compare_maf_res = compare_maf(locus)
        compare_maf_res["cohort"] = f"{locus.popu}_{locus.cohort}"
        all_expected_z.append(expected_z)
        all_dentist_s.append(dentist_s)
        all_compare_maf.append(compare_maf_res)
    all_expected_z = pd.concat(all_expected_z, axis=0)
    all_dentist_s = pd.concat(all_dentist_s, axis=0)
    all_compare_maf = pd.concat(all_compare_maf, axis=0)
    qc_metrics["expected_z"] = all_expected_z
    qc_metrics["dentist_s"] = all_dentist_s
    qc_metrics["compare_maf"] = all_compare_maf

    qc_metrics["ld_4th_moment"] = ld_4th_moment(locus_set)
    qc_metrics["ld_decay"] = ld_decay(locus_set)

    if len(locus_set.loci) > 1:
        qc_metrics["cochran_q"] = cochran_q(locus_set)
        qc_metrics["snp_missingness"] = snp_missingness(locus_set)

    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        for metric_name, metric_data in qc_metrics.items():
            metric_data.to_csv(f"{out_dir}/{metric_name}.txt", sep="\t", index=False)

    return qc_metrics


def qc_locus_cli(args: Tuple[str, pd.DataFrame, str]) -> str:
    """
    Quality control for a single locus (command-line interface wrapper).

    Parameters
    ----------
    args : Tuple[str, pd.DataFrame, str]
        Tuple containing:
        - locus_id : str
            Locus identifier
        - locus_info : pd.DataFrame
            DataFrame with locus information
        - base_out_dir : str
            Base output directory

    Returns
    -------
    str
        The locus_id that was processed.

    Notes
    -----
    This function is designed for multiprocessing and:

    1. Loads the locus set from the provided information
    2. Performs comprehensive QC analysis
    3. Creates locus-specific output directory
    4. Saves all QC results as compressed files
    5. Returns the processed locus_id for tracking

    Output files are saved as:
    {base_out_dir}/{locus_id}/{qc_metric}.txt.gz
    """
    locus_id, locus_info, base_out_dir = args
    locus_set = load_locus_set(locus_info)
    qc_metrics = locus_qc(locus_set)
    locus_out_dir = f"{base_out_dir}/{locus_id}"
    os.makedirs(locus_out_dir, exist_ok=True)
    for metric_name, metric_data in qc_metrics.items():
        metric_data.to_csv(
            f"{locus_out_dir}/{metric_name}.txt.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    return locus_id


def loci_qc(inputs: str, out_dir: str, threads: int = 1) -> None:
    """
    Perform quality control analysis on multiple loci in parallel.

    Parameters
    ----------
    inputs : str
        Path to input file containing locus information.
        Must be tab-separated with columns including 'locus_id'.
    out_dir : str
        Output directory path where results will be saved.
    threads : int, optional
        Number of parallel threads to use, by default 1.

    Returns
    -------
    None
        Results are saved to files in the output directory.

    Raises
    ------
    ValueError
        If the number of threads is less than 1.

    Notes
    -----
    This function processes multiple loci in parallel with the following workflow:

    1. Reads locus information from input file
    2. Groups loci by locus_id
    3. Processes each locus group using multiprocessing
    4. Displays progress bar for user feedback
    5. Saves results organized by locus_id

    The input file should contain columns: locus_id, prefix, popu, cohort, sample_size.

    Output structure:
    {out_dir}/{locus_id}/{qc_metric}.txt.gz

    Each locus gets its own subdirectory with compressed QC result files.
    """
    loci_info = pd.read_csv(inputs, sep="\t")

    # Create progress bar
    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
    )

    # Prepare arguments for multiprocessing
    locus_groups = [
        (locus_id, locus_info, out_dir)
        for locus_id, locus_info in loci_info.groupby("locus_id")
    ]

    with progress:
        task = progress.add_task("[cyan]Processing loci...", total=len(locus_groups))

        # Process loci in parallel with progress updates
        with Pool(threads) as pool:
            for _ in pool.imap_unordered(qc_locus_cli, locus_groups):  # type: ignore
                progress.update(task, advance=1)
