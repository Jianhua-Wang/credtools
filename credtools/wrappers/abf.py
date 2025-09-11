"""Wrapper of ABF fine-mapping method."""

import json
import logging
from typing import List

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet, combine_creds
from credtools.locus import Locus

logger = logging.getLogger("ABF")


def run_abf(
    locus: Locus, max_causal: int = 1, coverage: float = 0.95, var_prior: float = 0.2
) -> CredibleSet:
    """
    Run Approximate Bayes Factor (ABF) fine-mapping analysis.

    Calculate the approximate Bayes factor (ABF) from BETA and SE, using the
    formula:
    SNP_BF = sqrt(SE²/(SE² + W²)) * EXP(W²/(SE² + W²) * (BETA²/SE²)/2)
    where W is variance prior, usually set to 0.15 for quantitative traits
    and 0.2 for binary traits.

    The posterior probability of each variant being causal is calculated
    using the formula:
    PP(causal) = SNP_BF / sum(all_SNP_BFs)

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics for fine-mapping analysis.
    max_causal : int, optional
        Maximum number of causal variants, by default 1. ABF only supports
        single causal variant analysis, so this is always set to 1.
    coverage : float, optional
        Coverage probability for the credible set, by default 0.95.
        This determines the probability mass included in the credible set.
    var_prior : float, optional
        Variance prior parameter (W²), by default 0.2. This parameter controls
        the expected effect size:
        - 0.15 typically used for quantitative traits
        - 0.2 typically used for binary traits
        - Higher values assume larger effect sizes

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Posterior inclusion probabilities (PIPs) for all variants
        - Credible set variants that explain the specified coverage
        - Lead SNP with smallest p-value within credible set
        - Method-specific parameters and metadata

    Warnings
    --------
    If max_causal > 1, a warning is logged and max_causal is automatically set to 1,
    as ABF only supports single causal variant analysis.

    If no SNPs have p-value ≤ 1e-5, a warning is logged and an empty credible set
    is returned.

    Notes
    -----
    The ABF method assumes a single causal variant per locus and calculates
    Bayes factors for each variant independently. The method:

    1. Computes Bayes factors for each variant using effect size and standard error
    2. Normalizes Bayes factors to obtain posterior inclusion probabilities
    3. Selects variants for credible set based on ranked PIPs until coverage is reached
    4. Identifies lead SNP as the variant with smallest p-value in credible set

    The variance prior (W²) is a key parameter that affects the results:
    - Larger values favor variants with larger effect sizes
    - Smaller values are more conservative
    - Should be chosen based on trait type and expected effect sizes

    Reference:
    Asimit, J. L. et al. Eur J Hum Genet (2016).
    "Stochastic search and joint fine-mapping increases accuracy and identifies
    previously unreported associations in immune-mediated diseases"

    Examples
    --------
    >>> # Basic ABF analysis with default parameters
    >>> credible_set = run_abf(locus)
    >>> print(f"Found {credible_set.n_cs} credible set with {len(credible_set.snps[0])} variants")
    Found 1 credible set with 15 variants

    >>> # ABF analysis with custom variance prior for quantitative trait
    >>> credible_set = run_abf(locus, var_prior=0.15, coverage=0.99)
    >>> print(f"Coverage: {credible_set.coverage}")
    >>> print(f"Top PIP: {credible_set.pips.max():.4f}")
    Coverage: 0.99
    Top PIP: 0.6543

    >>> # Access posterior inclusion probabilities
    >>> pips_df = credible_set.pips.reset_index()
    >>> pips_df.columns = ['SNPID', 'PIP']
    >>> top_variants = pips_df.nlargest(5, 'PIP')
    >>> print(top_variants)
        SNPID           PIP
    0   rs123456    0.6543
    1   rs789012    0.1234
    2   rs345678    0.0987
    3   rs456789    0.0654
    4   rs567890    0.0321
    """
    if max_causal > 1:
        logger.warning(
            "ABF only support single causal variant. max_causal is set to 1."
        )
        max_causal = 1
    logger.info(f"Running ABF on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "var_prior": var_prior,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    df = locus.original_sumstats.copy()
    df["W2"] = var_prior**2
    df["SNP_BF"] = np.sqrt(
        (df[ColName.SE] ** 2 / (df[ColName.SE] ** 2 + df["W2"]))
    ) * np.exp(
        df["W2"]
        / (df[ColName.SE] ** 2 + df["W2"])
        * (df[ColName.BETA] ** 2 / df[ColName.SE] ** 2)
        / 2
    )
    df[ColName.PIP] = df["SNP_BF"] / df["SNP_BF"].sum()
    pips = pd.Series(
        data=df[ColName.PIP].values, index=df[ColName.SNPID].tolist(), name=ColName.ABF
    )
    if len(df[df[ColName.P] <= 1e-5]) > 0:
        ordering = np.argsort(pips.to_numpy())[::-1]
        idx = np.where(np.cumsum(pips.to_numpy()[ordering]) > coverage)[0][0]
        cs_snps = pips.index[ordering][: (idx + 1)].to_list()
        lead_snps = [
            str(
                df.loc[
                    df[df[ColName.SNPID].isin(cs_snps)][ColName.P].idxmin(),
                    ColName.SNPID,
                ]
            )
        ]
    else:
        logger.warning(
            "There are no SNPs with p-value <= 1e-5, output zero credible set"
        )
        cs_snps = []
        lead_snps = []
    logger.info(f"Finished ABF on {locus}")
    logger.info(f"N of credible set: {len(lead_snps)}")
    logger.info(f"Credible set size: [{len(cs_snps)}]")
    return CredibleSet(
        tool=Method.ABF,
        n_cs=1 if len(cs_snps) > 0 else 0,
        coverage=coverage,
        lead_snps=[lead_snps],  # type: ignore
        snps=[cs_snps] if len(cs_snps) > 0 else [],
        cs_sizes=[len(cs_snps)] if len(cs_snps) > 0 else [],
        pips=pips,
        parameters=parameters,
    )
