"""Wrapper for ABF+COJO combined fine-mapping method."""

import json
import logging
import tempfile
from typing import List, Optional

import numpy as np
import pandas as pd
from cojopy.cojopy import COJO

from credtools.cojo import conditional_selection
from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import Locus
from credtools.wrappers.abf import run_abf

logger = logging.getLogger("ABF_COJO")


def run_abf_cojo(
    locus: Locus,
    max_causal: int = 10,
    coverage: float = 0.95,
    var_prior: float = 0.2,
    p_cutoff: float = 5e-8,
    collinear_cutoff: float = 0.9,
    window_size: int = 10000000,
    maf_cutoff: float = 0.01,
    diff_freq_cutoff: float = 0.2,
) -> CredibleSet:
    """
    Run ABF fine-mapping with COJO conditional analysis for multi-signal loci.

    This method combines COJO stepwise selection with ABF fine-mapping to handle
    multi-signal loci. The approach:
    1. Uses COJO to detect independent signals at the locus
    2. If only one signal: runs standard ABF
    3. If multiple signals: performs conditional analysis for each signal and
       runs ABF on each conditional summary statistic

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix for analysis.
    max_causal : int, optional
        Maximum number of causal variants per signal, by default 10.
        This is also used as the limit for the number of independent signals
        to detect via COJO.
    coverage : float, optional
        Coverage probability for each credible set, by default 0.95.
    var_prior : float, optional
        Variance prior parameter (W²) for ABF calculation, by default 0.2.
        - 0.15 typically used for quantitative traits
        - 0.2 typically used for binary traits
    p_cutoff : float, optional
        P-value cutoff for COJO conditional selection, by default 5e-8.
    collinear_cutoff : float, optional
        Collinearity cutoff for COJO, by default 0.9.
    window_size : int, optional
        Window size for COJO conditional analysis, by default 10000000.
    maf_cutoff : float, optional
        Minor allele frequency cutoff for COJO, by default 0.01.
    diff_freq_cutoff : float, optional
        Allele frequency difference cutoff between summary stats and LD reference,
        by default 0.2.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Multiple credible sets (one per independent signal)
        - Posterior inclusion probabilities for all variants
        - Lead SNPs for each credible set
        - Combined PIPs across all signals

    Notes
    -----
    The method leverages COJO's ability to:
    1. Identify independent association signals through stepwise selection
    2. Compute conditional summary statistics for each signal
    3. Run ABF on each conditional dataset to get signal-specific credible sets

    For single-signal loci, results should be identical to standard ABF.
    For multi-signal loci, results should be similar to SuSiE but using ABF
    assumptions (single causal variant per signal).

    The final credible set object contains multiple credible sets, where each
    represents one independent signal. PIPs are computed per signal and then
    the maximum PIP across signals is taken for each SNP.

    Examples
    --------
    >>> # Basic ABF+COJO analysis
    >>> credible_set = run_abf_cojo(locus)
    >>> print(f"Found {credible_set.n_cs} independent signals")
    Found 3 independent signals

    >>> # Access credible sets for each signal
    >>> for i, cs_snps in enumerate(credible_set.snps):
    ...     print(f"Signal {i+1}: {len(cs_snps)} variants")
    Signal 1: 12 variants
    Signal 2: 8 variants
    Signal 3: 15 variants

    >>> # Get lead SNP for each signal
    >>> for i, lead_snp in enumerate(credible_set.lead_snps):
    ...     print(f"Signal {i+1} lead SNP: {lead_snp}")
    Signal 1 lead SNP: rs123456
    Signal 2 lead SNP: rs789012
    Signal 3 lead SNP: rs345678
    """
    logger.info(f"Running ABF+COJO on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "var_prior": var_prior,
        "p_cutoff": p_cutoff,
        "collinear_cutoff": collinear_cutoff,
        "window_size": window_size,
        "maf_cutoff": maf_cutoff,
        "diff_freq_cutoff": diff_freq_cutoff,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")

    # Step 1: Use COJO to detect independent signals
    logger.info(
        "Step 1: Running COJO conditional selection to detect independent signals"
    )
    cojo_results = conditional_selection(
        locus,
        p_cutoff=p_cutoff,
        collinear_cutoff=collinear_cutoff,
        window_size=window_size,
        maf_cutoff=maf_cutoff,
        diff_freq_cutoff=diff_freq_cutoff,
    )

    n_signals = len(cojo_results)
    logger.info(f"COJO detected {n_signals} independent signals")

    # Step 2: Handle based on number of signals detected
    if n_signals == 0:
        logger.warning("No independent signals detected by COJO")
        # Return empty credible set
        empty_pips = pd.Series(
            data=np.zeros(len(locus.sumstats)),
            index=locus.sumstats[ColName.SNPID].tolist(),
            name=f"{Method.ABF}_COJO",
        )
        return CredibleSet(
            tool=f"{Method.ABF}_COJO",
            n_cs=0,
            coverage=coverage,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=empty_pips,
            parameters=parameters,
        )

    elif n_signals == 1:
        logger.info("Single signal detected, running standard ABF")
        # For single signal, just run standard ABF
        abf_result = run_abf(
            locus, max_causal=1, coverage=coverage, var_prior=var_prior
        )
        # Update the method name and parameters
        abf_result._tool = f"{Method.ABF}_COJO"
        abf_result._parameters = parameters
        abf_result._pips.name = f"{Method.ABF}_COJO"
        return abf_result

    else:
        logger.info(
            f"Multiple signals detected ({n_signals}), running conditional ABF analysis"
        )
        # Step 3: For each signal, perform conditional analysis and run ABF
        return _run_conditional_abf_analysis(
            locus, cojo_results, coverage, var_prior, parameters
        )


def _run_conditional_abf_analysis(
    locus: Locus,
    cojo_results: pd.DataFrame,
    coverage: float,
    var_prior: float,
    parameters: dict,
) -> CredibleSet:
    """
    Run conditional ABF analysis for multiple independent signals.

    For each signal detected by COJO, this function:
    1. Creates conditional summary statistics (conditioning on other signals)
    2. Runs ABF on the conditional statistics
    3. Combines results into a multi-credible set structure
    """
    # Prepare COJO object for conditional analysis
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

    ld_matrix = locus.ld.r.copy()
    ld_freq: Optional[pd.DataFrame] = locus.ld.map.copy()
    if ld_freq is not None and "AF2" not in ld_freq.columns:
        logger.warning("AF2 is not in the LD matrix.")
        ld_freq = None
    elif ld_freq is not None:
        ld_freq = ld_freq[["SNPID", "AF2"]]
        ld_freq.columns = ["SNP", "freq"]
        ld_freq["freq"] = 1 - ld_freq["freq"]

    # Initialize COJO object
    c = COJO(
        p_cutoff=parameters["p_cutoff"],
        collinear_cutoff=parameters["collinear_cutoff"],
        window_size=parameters["window_size"],
        maf_cutoff=parameters["maf_cutoff"],
        diff_freq_cutoff=parameters["diff_freq_cutoff"],
    )
    c.load_sumstats(sumstats=sumstats, ld_matrix=ld_matrix, ld_freq=ld_freq)

    # Get list of independent SNPs from COJO results
    independent_snps = cojo_results["SNP"].tolist()
    n_signals = len(independent_snps)

    all_credible_sets = []
    all_lead_snps = []
    all_cs_sizes = []
    combined_pips = pd.Series(
        data=np.zeros(len(locus.sumstats)),
        index=locus.sumstats[ColName.SNPID].tolist(),
        name=f"{Method.ABF}_COJO",
    )

    logger.info(f"Running conditional ABF analysis for {n_signals} signals")

    for i, signal_snp in enumerate(independent_snps):
        logger.info(f"Processing signal {i+1}/{n_signals}: {signal_snp}")

        # Create list of conditioning SNPs (all other independent SNPs)
        conditioning_snps = [snp for snp in independent_snps if snp != signal_snp]

        if len(conditioning_snps) > 0:
            # Create temporary file with conditioning SNPs for COJO
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".txt"
            ) as f:
                for snp in conditioning_snps:
                    f.write(f"{snp}\n")
                cond_snps_file = f.name

            try:
                # Run conditional analysis
                logger.debug(
                    f"Running conditional analysis conditioning on: {conditioning_snps}"
                )
                conditional_results = c.run_conditional_analysis(
                    cond_snps_path=cond_snps_file
                )

                # Create modified locus with conditional summary statistics
                conditional_locus = _create_conditional_locus(
                    locus, conditional_results, signal_snp
                )

            finally:
                # Clean up temporary file
                import os

                os.unlink(cond_snps_file)
        else:
            # No conditioning needed (shouldn't happen with multiple signals, but safety check)
            conditional_locus = locus

        # Run ABF on conditional locus
        signal_abf_result = run_abf(
            conditional_locus, max_causal=1, coverage=coverage, var_prior=var_prior
        )

        # Extract results for this signal
        if signal_abf_result.n_cs > 0:
            all_credible_sets.extend(signal_abf_result.snps)
            all_lead_snps.extend(signal_abf_result.lead_snps)
            all_cs_sizes.extend(signal_abf_result.cs_sizes)

            # Update combined PIPs (take maximum across signals)
            signal_pips = signal_abf_result.pips
            for snp in signal_pips.index:
                if snp in combined_pips.index:
                    combined_pips[snp] = max(combined_pips[snp], signal_pips[snp])
        else:
            logger.warning(f"No credible set found for signal {i+1}: {signal_snp}")

    n_cs = len(all_credible_sets)
    logger.info(
        f"ABF+COJO analysis complete: {n_cs} credible sets from {n_signals} signals"
    )

    return CredibleSet(
        tool=f"{Method.ABF}_COJO",
        n_cs=n_cs,
        coverage=coverage,
        lead_snps=all_lead_snps,
        snps=all_credible_sets,
        cs_sizes=all_cs_sizes,
        pips=combined_pips,
        parameters=parameters,
    )


def _create_conditional_locus(
    original_locus: Locus, conditional_results: pd.DataFrame, signal_snp: str
) -> Locus:
    """
    Create a new Locus object with conditional summary statistics.

    This function updates the summary statistics with conditional effects
    while preserving the LD matrix and other locus properties.
    """
    # Create a copy of the original locus
    from copy import deepcopy

    conditional_locus = deepcopy(original_locus)

    # Update summary statistics with conditional results
    original_sumstats = conditional_locus.sumstats.copy()

    # Map conditional results back to original format
    # conditional_results has columns: SNP, original_beta, original_se, original_p, cond_beta, cond_se, cond_p
    conditional_map = conditional_results.set_index("SNP")[
        ["cond_beta", "cond_se", "cond_p"]
    ].to_dict()

    # Update the summary statistics
    for idx, row in original_sumstats.iterrows():
        snpid = row[ColName.SNPID]
        if snpid in conditional_map["cond_beta"]:
            original_sumstats.loc[idx, ColName.BETA] = conditional_map["cond_beta"][
                snpid
            ]
            original_sumstats.loc[idx, ColName.SE] = conditional_map["cond_se"][snpid]
            original_sumstats.loc[idx, ColName.P] = conditional_map["cond_p"][snpid]

    # Update the locus with new summary statistics
    conditional_locus._sumstats = original_sumstats

    # Also update original_sumstats if it exists
    if hasattr(conditional_locus, "_original_sumstats"):
        conditional_locus._original_sumstats = original_sumstats.copy()

    return conditional_locus
