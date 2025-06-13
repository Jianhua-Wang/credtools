"""Wrapper for COJO."""

import logging
from typing import Optional

import pandas as pd
from cojopy.cojopy import COJO

from credtools.locus import Locus
from credtools.sumstats import ColName

logger = logging.getLogger("COJO")


def conditional_selection(
    locus: Locus,
    p_cutoff: float = 5e-8,
    collinear_cutoff: float = 0.9,
    window_size: int = 10000000,
    maf_cutoff: float = 0.01,
    diff_freq_cutoff: float = 0.2,
) -> pd.DataFrame:
    """
    Perform conditional selection on the locus using COJO method.

    Parameters
    ----------
    locus : Locus
        The locus to perform conditional selection on. Must contain summary statistics
        and LD matrix data.
    p_cutoff : float, optional
        The p-value cutoff for the conditional selection, by default 5e-8.
        If no SNPs pass this threshold, it will be relaxed to 1e-5.
    collinear_cutoff : float, optional
        The collinearity cutoff for the conditional selection, by default 0.9.
        SNPs with LD correlation above this threshold are considered collinear.
    window_size : int, optional
        The window size in base pairs for the conditional selection, by default 10000000.
        SNPs within this window are considered for conditional analysis.
    maf_cutoff : float, optional
        The minor allele frequency cutoff for the conditional selection, by default 0.01.
        SNPs with MAF below this threshold are excluded.
    diff_freq_cutoff : float, optional
        The difference in frequency cutoff between summary statistics and reference panel,
        by default 0.2. SNPs with frequency differences above this threshold are excluded.

    Returns
    -------
    pd.DataFrame
        The conditional selection results containing independently associated variants
        with columns including SNP identifiers, effect sizes, and conditional p-values.

    Warnings
    --------
    If no SNPs pass the initial p-value cutoff, the threshold is automatically
    relaxed to 1e-5 and a warning is logged.

    If AF2 (reference allele frequency) is not available in the LD matrix,
    a warning is logged and frequency checking is disabled.

    Notes
    -----
    COJO (Conditional and Joint analysis) performs stepwise conditional analysis
    to identify independently associated variants at a locus. The method:

    1. Identifies the most significant SNP
    2. Performs conditional analysis on remaining SNPs
    3. Iteratively adds independently associated SNPs
    4. Continues until no more SNPs meet significance criteria

    The algorithm accounts for linkage disequilibrium patterns and helps
    distinguish truly independent signals from those in LD with lead variants.

    Reference: Yang, J. et al. Conditional and joint multiple-SNP analysis of GWAS
    summary statistics identifies additional variants influencing complex traits.
    Nat Genet 44, 369-375 (2012).

    Examples
    --------
    >>> # Basic conditional selection
    >>> results = conditional_selection(locus)
    >>> print(f"Found {len(results)} independent signals")
    Found 3 independent signals

    >>> # With custom thresholds
    >>> results = conditional_selection(
    ...     locus,
    ...     p_cutoff=1e-6,
    ...     maf_cutoff=0.05
    ... )
    >>> print(results[['SNP', 'b', 'se', 'p']])
        SNP           b        se         p
    0   rs123456   0.15   0.025   1.2e-08
    1   rs789012  -0.08   0.020   4.5e-07
    """
    sumstats = locus.sumstats.copy()
    sumstats = sumstats[
        [
            ColName.SNPID,
            ColName.EA,
            ColName.NEA,
            ColName.BETA,
            ColName.SE,
            ColName.P,
            ColName.EAF,
        ]
    ]
    sumstats.columns = ["SNP", "A1", "A2", "b", "se", "p", "freq"]
    sumstats["N"] = locus.sample_size
    if p_cutoff < 1e-5 and len(sumstats[sumstats["p"] < p_cutoff]) == 0:
        logger.warning("No SNPs passed the p-value cutoff, using p_cutoff=1e-5")
        p_cutoff = 1e-5

    ld_matrix = locus.ld.r.copy()
    ld_freq: Optional[pd.DataFrame] = locus.ld.map.copy()
    if ld_freq is not None and "AF2" not in ld_freq.columns:
        logger.warning("AF2 is not in the LD matrix.")
        ld_freq = None
    elif ld_freq is not None:
        ld_freq = ld_freq[["SNPID", "AF2"]]
        ld_freq.columns = ["SNP", "freq"]
        ld_freq["freq"] = 1 - ld_freq["freq"]
    c = COJO(
        p_cutoff=p_cutoff,
        collinear_cutoff=collinear_cutoff,
        window_size=window_size,
        maf_cutoff=maf_cutoff,
        diff_freq_cutoff=diff_freq_cutoff,
    )
    c.load_sumstats(sumstats=sumstats, ld_matrix=ld_matrix, ld_freq=ld_freq)  # type: ignore
    cojo_result = c.conditional_selection()
    return cojo_result
