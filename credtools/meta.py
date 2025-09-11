"""Meta analysis of multi-ancestry gwas data."""

import logging
import os
from multiprocessing import Pool
from typing import Any, Dict, List, Optional, Tuple

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

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet, intersect_sumstat_ld, load_locus_set
from credtools.sumstats import munge

logger = logging.getLogger("META")


def meta_sumstats(inputs: LocusSet) -> pd.DataFrame:
    """
    Perform fixed effect meta-analysis of summary statistics.

    Parameters
    ----------
    inputs : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    pd.DataFrame
        Meta-analysis summary statistics with columns: SNPID, BETA, SE, P, EAF, CHR, BP, EA, NEA.

    Notes
    -----
    This function performs inverse-variance weighted fixed-effects meta-analysis:

    1. Merges summary statistics from all studies on SNPID
    2. Calculates inverse-variance weights (1/SE²)
    3. Computes weighted average effect size
    4. Calculates meta-analysis standard error
    5. Computes Z-scores and p-values
    6. Performs sample-size weighted averaging of effect allele frequencies

    The meta-analysis formulas used:
    - Beta_meta = Σ(Beta_i * Weight_i) / Σ(Weight_i)
    - SE_meta = 1 / sqrt(Σ(Weight_i))
    - Weight_i = 1 / SE_i²
    """
    # Merge all dataframes on SNPID
    merged_df = inputs.loci[0].original_sumstats[[ColName.SNPID]].copy()
    n_sum = sum([input.sample_size for input in inputs.loci])
    eaf_weights = [input.sample_size / n_sum for input in inputs.loci]
    for i, df in enumerate(inputs.loci):
        df = df.sumstats[[ColName.SNPID, ColName.BETA, ColName.SE, ColName.EAF]].copy()
        df.rename(
            columns={
                ColName.BETA: f"BETA_{i}",
                ColName.SE: f"SE_{i}",
                ColName.EAF: f"EAF_{i}",
            },
            inplace=True,
        )
        merged_df = pd.merge(
            merged_df, df, on=ColName.SNPID, how="outer", suffixes=("", f"_{i}")
        )

    # Calculate weights (inverse of variance)
    for i in range(len(inputs.loci)):
        merged_df[f"weight_{i}"] = 1 / (merged_df[f"SE_{i}"] ** 2)
        # merged_df[f"EAF_{i}"] = merged_df[f"EAF_{i}"] * eaf_weights[i]

    merged_df.fillna(0, inplace=True)

    # Calculate meta-analysis beta
    beta_numerator = sum(
        merged_df[f"BETA_{i}"] * merged_df[f"weight_{i}"]
        for i in range(len(inputs.loci))
    )
    weight_sum = sum(merged_df[f"weight_{i}"] for i in range(len(inputs.loci)))
    meta_beta = beta_numerator / weight_sum

    # Calculate meta-analysis SE
    meta_se = np.sqrt(1 / weight_sum)

    # Calculate meta-analysis Z-score and p-value
    meta_z = meta_beta / meta_se
    meta_p = 2 * stats.norm.sf(abs(meta_z))

    # Calculate meta-analysis EAF
    meta_eaf = sum(
        merged_df[f"EAF_{i}"] * eaf_weights[i] for i in range(len(inputs.loci))
    )

    # Create output dataframe
    output_df = pd.DataFrame(
        {
            ColName.SNPID: merged_df[ColName.SNPID],
            ColName.BETA: meta_beta,
            ColName.SE: meta_se,
            ColName.P: meta_p,
            ColName.EAF: meta_eaf,
        }
    )
    output_df[[ColName.CHR, ColName.BP, ColName.EA, ColName.NEA]] = merged_df[
        ColName.SNPID
    ].str.split("-", expand=True)[[0, 1, 2, 3]]
    return munge(output_df)


def meta_lds(inputs: LocusSet) -> LDMatrix:
    """
    Perform meta-analysis of LD matrices using sample-size weighted averaging.

    Parameters
    ----------
    inputs : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    LDMatrix
        Meta-analyzed LD matrix with merged variant map.

    Notes
    -----
    This function performs the following operations:

    1. Identifies unique variants across all studies
    2. Creates a master variant list sorted by chromosome and position
    3. Performs sample-size weighted averaging of LD correlations
    4. Handles missing variants by setting weights to zero
    5. Optionally meta-analyzes allele frequencies if available

    The meta-analysis formula:
    LD_meta[i,j] = Σ(LD_k[i,j] * N_k) / Σ(N_k)

    where k indexes studies, N_k is sample size, and the sum is over studies
    that have both variants i and j.
    """
    # Get unique variants across all studies
    variant_dfs = [input.ld.map for input in inputs.loci]
    ld_matrices = [input.ld.r for input in inputs.loci]
    sample_sizes = [input.sample_size for input in inputs.loci]

    # Concatenate all variants
    merged_variants = pd.concat(variant_dfs, ignore_index=True)
    merged_variants.drop_duplicates(subset=[ColName.SNPID], inplace=True)
    merged_variants.sort_values([ColName.CHR, ColName.BP], inplace=True)
    merged_variants.reset_index(drop=True, inplace=True)
    # meta allele frequency of LD reference, if exists
    if all("AF2" in variant_df.columns for variant_df in variant_dfs):
        n_sum = sum([input.sample_size for input in inputs.loci])
        weights = [input.sample_size / n_sum for input in inputs.loci]
        af_df = merged_variants[[ColName.SNPID]].copy()
        af_df.set_index(ColName.SNPID, inplace=True)
        for i, variant_df in enumerate(variant_dfs):
            df = variant_df.copy()
            df.set_index(ColName.SNPID, inplace=True)
            af_df[f"AF2_{i}"] = df["AF2"] * weights[i]
        af_df.fillna(0, inplace=True)
        af_df["AF2_meta"] = af_df.sum(axis=1)
        merged_variants["AF2"] = merged_variants[ColName.SNPID].map(af_df["AF2_meta"])
    all_variants = merged_variants[ColName.SNPID].values
    variant_to_index = {snp: idx for idx, snp in enumerate(all_variants)}
    n_variants = len(all_variants)

    # Initialize arrays using numpy operations
    merged_ld = np.zeros((n_variants, n_variants))
    weight_matrix = np.zeros((n_variants, n_variants))

    # Process each study
    for ld_mat, variants_df, sample_size in zip(ld_matrices, variant_dfs, sample_sizes):
        # coverte float16 to float32, to avoid overflow
        # ld_mat = ld_mat.astype(np.float32)

        # Get indices in the master matrix
        study_snps = variants_df["SNPID"].values
        study_indices = np.array([variant_to_index[snp] for snp in study_snps])

        # Create index meshgrid for faster indexing
        idx_i, idx_j = np.meshgrid(study_indices, study_indices)

        # Update matrices using vectorized operations
        merged_ld[idx_i, idx_j] += ld_mat * sample_size
        weight_matrix[idx_i, idx_j] += sample_size

    # Compute weighted average
    mask = weight_matrix != 0
    merged_ld[mask] /= weight_matrix[mask]

    return LDMatrix(merged_variants, merged_ld.astype(np.float32))


def meta_all(inputs: LocusSet) -> Locus:
    """
    Perform comprehensive meta-analysis of both summary statistics and LD matrices.

    Parameters
    ----------
    inputs : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    Locus
        Meta-analyzed Locus object with combined population and cohort identifiers.

    Notes
    -----
    This function:

    1. Performs meta-analysis of summary statistics using inverse-variance weighting
    2. Performs meta-analysis of LD matrices using sample-size weighting
    3. Combines population and cohort names from all input studies
    4. Sums sample sizes across studies
    5. Intersects the meta-analyzed data to ensure consistency

    Population and cohort names are combined with "+" as separator and sorted alphabetically.
    """
    meta_sumstat = meta_sumstats(inputs)
    meta_ld = meta_lds(inputs)
    sample_size = sum([input.sample_size for input in inputs.loci])
    popu_set = set()
    for input in inputs.loci:
        for pop in input.popu.split(","):
            popu_set.add(pop)
    popu = "+".join(sorted(popu_set))
    cohort_set = set()
    for input in inputs.loci:
        for cohort_name in input.cohort.split(","):
            cohort_set.add(cohort_name)
    cohort = "+".join(sorted(cohort_set))

    return Locus(
        popu, cohort, sample_size, sumstats=meta_sumstat, ld=meta_ld, if_intersect=True
    )


def meta_by_population(inputs: LocusSet) -> Dict[str, Locus]:
    """
    Perform meta-analysis within each population separately.

    Parameters
    ----------
    inputs : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    Dict[str, Locus]
        Dictionary mapping population codes to meta-analyzed Locus objects.

    Notes
    -----
    This function:

    1. Groups studies by population code
    2. Performs meta-analysis within each population group
    3. For single-study populations, applies intersection without meta-analysis
    4. Returns a dictionary with population codes as keys

    This approach preserves population-specific LD patterns while still
    allowing meta-analysis of multiple cohorts within the same population.
    """
    meta_popu = {}
    for input in inputs.loci:
        popu = input.popu
        if popu not in meta_popu:
            meta_popu[popu] = [input]
        else:
            meta_popu[popu].append(input)

    result_dict = {}
    for popu in meta_popu:
        if len(meta_popu[popu]) > 1:
            result_dict[popu] = meta_all(LocusSet(meta_popu[popu]))
        else:
            result_dict[popu] = intersect_sumstat_ld(meta_popu[popu][0])
    return result_dict


def meta(inputs: LocusSet, meta_method: str = "meta_all") -> LocusSet:
    """
    Perform meta-analysis using the specified method.

    Parameters
    ----------
    inputs : LocusSet
        LocusSet containing input data from multiple studies.
    meta_method : str, optional
        Meta-analysis method to use, by default "meta_all".
        Options:
        - "meta_all": Meta-analyze all studies together
        - "meta_by_population": Meta-analyze within each population separately
        - "no_meta": No meta-analysis, just intersect individual studies

    Returns
    -------
    LocusSet
        LocusSet containing meta-analyzed results.

    Raises
    ------
    ValueError
        If an unsupported meta-analysis method is specified.

    Notes
    -----
    The different methods serve different purposes:

    - "meta_all": Maximizes power by combining all studies, but may be inappropriate
      if LD patterns differ substantially between populations
    - "meta_by_population": Preserves population-specific LD while allowing
      meta-analysis within populations
    - "no_meta": Keeps studies separate, useful for comparison or when
      meta-analysis is not appropriate
    """
    if meta_method == "meta_all":
        return LocusSet([meta_all(inputs)])
    elif meta_method == "meta_by_population":
        res = meta_by_population(inputs)
        return LocusSet([res[popu] for popu in res])
    elif meta_method == "no_meta":
        return LocusSet([intersect_sumstat_ld(i) for i in inputs.loci])
    else:
        raise ValueError(f"Unsupported meta-analysis method: {meta_method}")


def meta_locus(args: Tuple[str, pd.DataFrame, str, str, bool]) -> List[List[Any]]:
    """
    Process a single locus for meta-analysis.

    Parameters
    ----------
    args : Tuple[str, pd.DataFrame, str, str]
        A tuple containing:
        - locus_id : str
            The ID of the locus
        - locus_info : pd.DataFrame
            DataFrame containing locus information
        - outdir : str
            Output directory path
        - meta_method : str
            Method for meta-analysis

    Returns
    -------
    List[List[Any]]
        A list of results containing processed locus information.
        Each inner list contains: [chrom, start, end, popu, sample_size, cohort, out_prefix, locus_id]

    Notes
    -----
    This function is designed for parallel processing and:

    1. Loads the locus set from the provided information
    2. Performs meta-analysis using the specified method
    3. Creates output directory for the locus
    4. Saves results to compressed files (sumstats.gz, ld.npz, ldmap.gz)
    5. Returns metadata for each processed locus
    """
    locus_id, locus_info, outdir, meta_method, calculate_lambda_s = args
    results = []
    locus_set = load_locus_set(locus_info, calculate_lambda_s=calculate_lambda_s)
    locus_set = meta(locus_set, meta_method)
    out_dir = f"{outdir}/{locus_id}"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    for locus in locus_set.loci:
        out_prefix = f"{out_dir}/{locus.prefix}"
        locus.sumstats.to_csv(
            f"{out_prefix}.sumstats.gz", sep="\t", index=False, compression="gzip"
        )
        np.savez_compressed(f"{out_prefix}.ld.npz", ld=locus.ld.r.astype(np.float16))
        locus.ld.map.to_csv(
            f"{out_prefix}.ldmap.gz", sep="\t", index=False, compression="gzip"
        )
        chrom, start, end = locus.chrom, locus.start, locus.end
        results.append(
            [
                chrom,
                start,
                end,
                locus.popu,
                locus.sample_size,
                locus.cohort,
                out_prefix,
                f"chr{chrom}_{start}_{end}",
            ]
        )
    return results


def meta_loci(
    inputs: str,
    outdir: str,
    threads: int = 1,
    meta_method: str = "meta_all",
    calculate_lambda_s: bool = False,
) -> None:
    """
    Perform meta-analysis on multiple loci in parallel.

    Parameters
    ----------
    inputs : str
        Path to input file containing locus information.
        Must be a tab-separated file with columns including 'locus_id'.
    outdir : str
        Output directory path where results will be saved.
    threads : int, optional
        Number of parallel threads to use, by default 1.
    meta_method : str, optional
        Meta-analysis method to use, by default "meta_all".
    calculate_lambda_s : bool, optional
        Whether to calculate lambda_s parameter using estimate_s_rss function, by default False.
        See meta() function for available options.

    Returns
    -------
    None
        Results are saved to files in the output directory.

    Notes
    -----
    This function:

    1. Reads locus information from the input file
    2. Groups loci by locus_id for parallel processing
    3. Processes each locus group using the specified meta-analysis method
    4. Saves results with a progress bar for user feedback
    5. Creates a summary file (loci_info.txt) with all processed loci

    The input file should contain columns: locus_id, prefix, popu, cohort, sample_size.
    Each locus_id can have multiple rows representing different cohorts/populations.

    Output files are organized as:
    {outdir}/{locus_id}/{prefix}.{sumstats.gz,ld.npz,ldmap.gz}
    """
    loci_info = pd.read_csv(inputs, sep="\t")
    new_loci_info = pd.DataFrame(columns=loci_info.columns)

    # Group loci by locus_id
    grouped_loci = list(loci_info.groupby("locus_id"))
    total_loci = len(grouped_loci)
    os.makedirs(outdir, exist_ok=True)

    # Create process pool and process loci in parallel with progress bar
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task("[cyan]Meta-analysing...", total=total_loci)

        with Pool(threads) as pool:
            args = [
                (locus_id, locus_info, outdir, meta_method, calculate_lambda_s)
                for locus_id, locus_info in grouped_loci
            ]
            for result in pool.imap_unordered(meta_locus, args):  # type: ignore
                # Update results
                for i, res in enumerate(result):
                    new_loci_info.loc[len(new_loci_info)] = res
                # Update progress
                progress.advance(task)

    new_loci_info.to_csv(f"{outdir}/loci_info.txt", sep="\t", index=False)
