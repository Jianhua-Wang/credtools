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
from credtools.locus import (
    Locus,
    LocusSet,
    check_loci_info,
    intersect_sumstat_ld,
    load_locus_set,
)
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

    # Record which cohorts have EAF data before filling NaN
    eaf_present = pd.DataFrame(
        {i: merged_df[f"EAF_{i}"].notna() for i in range(len(inputs.loci))}
    )

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
    meta_p = np.maximum(2 * stats.norm.sf(abs(meta_z)), 1e-300)

    # Calculate meta-analysis EAF (normalize weights to present cohorts only)
    eaf_weight_sum = sum(
        eaf_present[i].astype(float) * eaf_weights[i] for i in range(len(inputs.loci))
    )
    meta_eaf = (
        sum(merged_df[f"EAF_{i}"] * eaf_weights[i] for i in range(len(inputs.loci)))
        / eaf_weight_sum
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
            af_df[f"AF2_{i}"] = df["AF2"]
        # Normalize weights to present cohorts only
        af_present = af_df.notna()
        af_df.fillna(0, inplace=True)
        weight_sum = sum(
            af_present[f"AF2_{i}"].astype(float) * weights[i]
            for i in range(len(variant_dfs))
        )
        af_df["AF2_meta"] = (
            sum(af_df[f"AF2_{i}"] * weights[i] for i in range(len(variant_dfs)))
            / weight_sum
        )
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

    # All input loci must have the same boundaries
    locus_starts = [locus._locus_start for locus in inputs.loci]
    locus_ends = [locus._locus_end for locus in inputs.loci]

    if not all(s == locus_starts[0] for s in locus_starts):
        raise ValueError("All input loci must have the same start position")
    if not all(e == locus_ends[0] for e in locus_ends):
        raise ValueError("All input loci must have the same end position")

    return Locus(
        popu,
        cohort,
        sample_size,
        meta_sumstat,
        locus_starts[0],
        locus_ends[0],
        ld=meta_ld,
        if_intersect=True,
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


def compute_heterogeneity(locus_set: LocusSet) -> Dict[str, pd.DataFrame]:
    """Compute heterogeneity metrics across cohorts before meta-analysis.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary of heterogeneity metrics with keys:
        - 'ld_4th_moment': 4th moment of LD matrix
        - 'ld_decay': LD decay analysis
        - 'cochran_q': heterogeneity test (if multiple cohorts)
        - 'snp_missingness': missingness analysis (if multiple cohorts)
    """
    from credtools.qc import cochran_q, ld_4th_moment, ld_decay, snp_missingness

    het_metrics: Dict[str, pd.DataFrame] = {}
    het_metrics["ld_4th_moment"] = ld_4th_moment(locus_set)
    het_metrics["ld_decay"] = ld_decay(locus_set)
    if len(locus_set.loci) > 1:
        het_metrics["cochran_q"] = cochran_q(locus_set)
        het_metrics["snp_missingness"] = snp_missingness(locus_set)
    return het_metrics


def compute_heterogeneity_by_population(
    locus_set: LocusSet,
) -> Dict[str, pd.DataFrame]:
    """Compute heterogeneity metrics with cochran_q and snp_missingness grouped by population.

    LD metrics (ld_4th_moment, ld_decay) are computed globally across all cohorts.
    Cochran-Q and SNP missingness are computed within each population separately,
    since heterogeneity should measure differences between cohorts of the same
    population, not across populations.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary of heterogeneity metrics with keys:
        - 'ld_4th_moment': 4th moment of LD matrix (global)
        - 'ld_decay': LD decay analysis (global)
        - 'cochran_q': heterogeneity test per population (with 'population' column)
        - 'snp_missingness': missingness per population (with 'population' column)
    """
    from credtools.qc import cochran_q, ld_4th_moment, ld_decay, snp_missingness

    het_metrics: Dict[str, pd.DataFrame] = {}
    het_metrics["ld_4th_moment"] = ld_4th_moment(locus_set)
    het_metrics["ld_decay"] = ld_decay(locus_set)

    # Group loci by population
    pop_groups: Dict[str, List[Locus]] = {}
    for locus in locus_set.loci:
        pop_groups.setdefault(locus.popu, []).append(locus)

    q_parts: List[pd.DataFrame] = []
    miss_parts: List[pd.DataFrame] = []
    for popu, loci in pop_groups.items():
        if len(loci) > 1:
            sub_set = LocusSet(loci)
            cq = cochran_q(sub_set)
            cq["population"] = popu
            q_parts.append(cq)
            ms = snp_missingness(sub_set)
            ms["population"] = popu
            miss_parts.append(ms)

    if q_parts:
        het_metrics["cochran_q"] = pd.concat(q_parts, ignore_index=True)
    if miss_parts:
        het_metrics["snp_missingness"] = pd.concat(miss_parts, ignore_index=True)

    return het_metrics


def heterogeneity_summary(
    het_metrics: Dict[str, pd.DataFrame],
    locus_set: LocusSet,
) -> pd.DataFrame:
    """Generate per-cohort summary of heterogeneity metrics.

    Parameters
    ----------
    het_metrics : Dict[str, pd.DataFrame]
        Dictionary of heterogeneity metrics from compute_heterogeneity().
    locus_set : LocusSet
        LocusSet containing input data from multiple studies.

    Returns
    -------
    pd.DataFrame
        Summary DataFrame with one row per cohort.
    """
    rows = []
    for locus in locus_set.loci:
        cohort_label = f"{locus.popu}_{locus.cohort}"
        row: Dict[str, Any] = {"popu": locus.popu, "cohort": locus.cohort}

        # LD 4th moment mean for this cohort
        ld_4th = het_metrics.get("ld_4th_moment")
        if ld_4th is not None and cohort_label in ld_4th.columns:
            row["ld_4th_moment_mean"] = float(ld_4th[cohort_label].mean())
        else:
            row["ld_4th_moment_mean"] = np.nan

        # LD decay rate for this cohort
        ld_dec = het_metrics.get("ld_decay")
        if ld_dec is not None:
            cohort_data = ld_dec[ld_dec["cohort"] == cohort_label]
            if not cohort_data.empty:
                row["ld_decay_rate"] = float(cohort_data["decay_rate"].iloc[0])
            else:
                row["ld_decay_rate"] = np.nan
        else:
            row["ld_decay_rate"] = np.nan

        # SNP missingness rate for this cohort
        miss = het_metrics.get("snp_missingness")
        if miss is not None and cohort_label in miss.columns:
            # Filter by population if column exists
            if "population" in miss.columns:
                miss_filtered = miss[miss["population"] == locus.popu]
            else:
                miss_filtered = miss
            if cohort_label in miss_filtered.columns and not miss_filtered.empty:
                row["missing_rate"] = float(
                    round(
                        1 - miss_filtered[cohort_label].sum() / miss_filtered.shape[0],
                        6,
                    )
                )
            else:
                row["missing_rate"] = np.nan
        else:
            row["missing_rate"] = np.nan

        # Cochran Q locus-level summary
        cq = het_metrics.get("cochran_q")
        if cq is not None and not cq.empty:
            # Filter by population if column exists
            if "population" in cq.columns:
                cq_filtered = cq[cq["population"] == locus.popu]
            else:
                cq_filtered = cq
            if not cq_filtered.empty:
                row["cochran_q_median"] = float(cq_filtered["Q"].median())
                row["i_squared_median"] = float(cq_filtered["I_squared"].median())
                row["n_het_snps"] = int((cq_filtered["Q_pvalue"] < 0.05).sum())
            else:
                row["cochran_q_median"] = np.nan
                row["i_squared_median"] = np.nan
                row["n_het_snps"] = np.nan
        else:
            row["cochran_q_median"] = np.nan
            row["i_squared_median"] = np.nan
            row["n_het_snps"] = np.nan

        rows.append(row)
    return pd.DataFrame(rows)


def save_heterogeneity(
    het_metrics: Dict[str, pd.DataFrame],
    out_dir: str,
    summary: Optional[pd.DataFrame] = None,
) -> None:
    """Save heterogeneity metrics to compressed TSV files.

    Parameters
    ----------
    het_metrics : Dict[str, pd.DataFrame]
        Dictionary of heterogeneity metrics from compute_heterogeneity().
    out_dir : str
        Output directory path.
    summary : pd.DataFrame, optional
        Per-cohort summary from heterogeneity_summary(). If provided,
        saved as heterogeneity.txt.gz.
    """
    os.makedirs(out_dir, exist_ok=True)
    for name, data in het_metrics.items():
        save_data = data
        # Drop population column from snp_missingness for file format compatibility
        if name == "snp_missingness" and "population" in data.columns:
            save_data = data.drop(columns=["population"])
        save_data.to_csv(
            f"{out_dir}/{name}.txt.gz",
            sep="\t",
            index=False,
            float_format="%.6f",
            compression="gzip",
        )
    if summary is not None and not summary.empty:
        summary.to_csv(
            f"{out_dir}/heterogeneity.txt.gz",
            sep="\t",
            index=False,
            float_format="%.6f",
            compression="gzip",
        )


def recover_completed_locus(
    locus_id: str,
    outdir: str,
    prev_loci_info: Optional[pd.DataFrame],
) -> Optional[Tuple[List[List[Any]], pd.DataFrame]]:
    """Recover a previously completed locus from existing output files.

    Parameters
    ----------
    locus_id : str
        The locus identifier.
    outdir : str
        Output directory path.
    prev_loci_info : Optional[pd.DataFrame]
        Previous loci_info DataFrame from a prior run, or None.

    Returns
    -------
    Optional[Tuple[List[List[Any]], pd.DataFrame]]
        (results, het_summary) if the locus is fully complete, else None.
    """
    if prev_loci_info is None:
        return None

    locus_dir = os.path.join(outdir, locus_id)
    if not os.path.isdir(locus_dir):
        return None

    rows = prev_loci_info[prev_loci_info["locus_id"] == locus_id]
    if rows.empty:
        return None

    # Check all 3 files exist for each prefix
    for _, row in rows.iterrows():
        prefix = row["prefix"]
        for ext in [".sumstats.gz", ".ld.npz", ".ldmap.gz"]:
            if not os.path.exists(f"{prefix}{ext}"):
                return None

    # Rebuild results list
    results = []
    for _, row in rows.iterrows():
        results.append(
            [
                row["chr"],
                row["start"],
                row["end"],
                row["popu"],
                row["sample_size"],
                row["cohort"],
                row["prefix"],
                row["locus_id"],
            ]
        )

    # Read heterogeneity summary if available
    het_path = os.path.join(locus_dir, "heterogeneity.txt.gz")
    if os.path.exists(het_path):
        het_summary = pd.read_csv(het_path, sep="\t", compression="gzip")
    else:
        het_summary = pd.DataFrame()

    return results, het_summary


def meta_locus(
    args: Tuple[str, pd.DataFrame, str, str, bool],
) -> Tuple[List[List[Any]], pd.DataFrame]:
    """Process a single locus for meta-analysis.

    Parameters
    ----------
    args : Tuple[str, pd.DataFrame, str, str, bool]
        A tuple containing:
        - locus_id : str
            The ID of the locus
        - locus_info : pd.DataFrame
            DataFrame containing locus information
        - outdir : str
            Output directory path
        - meta_method : str
            Method for meta-analysis
        - calculate_lambda_s : bool
            Whether to calculate lambda_s

    Returns
    -------
    Tuple[List[List[Any]], pd.DataFrame]
        A tuple of (results, het_summary) where results is a list of
        processed locus info and het_summary is the per-cohort heterogeneity
        summary DataFrame.

    Notes
    -----
    This function is designed for parallel processing and:

    1. Loads the locus set from the provided information
    2. Computes heterogeneity metrics before meta-analysis
    3. Performs meta-analysis using the specified method
    4. Creates output directory for the locus
    5. Saves results to compressed files (sumstats.gz, ld.npz, ldmap.gz)
    6. Returns metadata for each processed locus and heterogeneity summary
    """
    locus_id, locus_info, outdir, meta_method, calculate_lambda_s = args
    results = []
    locus_set = load_locus_set(locus_info, calculate_lambda_s=calculate_lambda_s)

    # Compute heterogeneity BEFORE meta combines data
    if meta_method == "meta_by_population":
        het_metrics = compute_heterogeneity_by_population(locus_set)
    else:
        het_metrics = compute_heterogeneity(locus_set)
    het_summary = heterogeneity_summary(het_metrics, locus_set)
    het_summary["locus_id"] = locus_id

    locus_set = meta(locus_set, meta_method)
    out_dir = f"{outdir}/{locus_id}"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Save heterogeneity alongside meta results
    save_heterogeneity(het_metrics, out_dir, summary=het_summary)

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
    return results, het_summary


def meta_loci(
    inputs: str,
    outdir: str,
    threads: int = 1,
    meta_method: str = "meta_all",
    calculate_lambda_s: bool = False,
    skip: bool = False,
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
    skip : bool, optional
        Skip loci already completed from a previous run, by default False.

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
    loci_info = check_loci_info(loci_info)  # Validate input data
    new_loci_info = pd.DataFrame(columns=loci_info.columns)
    all_het_summaries: List[pd.DataFrame] = []

    # Group loci by locus_id
    grouped_loci = list(loci_info.groupby("locus_id"))
    total_loci = len(grouped_loci)
    os.makedirs(outdir, exist_ok=True)

    # Try to recover completed loci when skip=True
    prev_loci_info: Optional[pd.DataFrame] = None
    pending_args: List[Tuple[str, pd.DataFrame, str, str, bool]] = []
    skipped_count = 0

    if skip:
        prev_path = os.path.join(outdir, "loci_info.txt")
        if os.path.exists(prev_path):
            prev_loci_info = pd.read_csv(prev_path, sep="\t")
            logger.info(f"Found previous loci_info.txt with {len(prev_loci_info)} rows")

    if prev_loci_info is not None:
        for locus_id, locus_info_group in grouped_loci:
            recovered = recover_completed_locus(locus_id, outdir, prev_loci_info)
            if recovered is not None:
                results, het_summary = recovered
                for res in results:
                    new_loci_info.loc[len(new_loci_info)] = res
                if het_summary is not None and not het_summary.empty:
                    all_het_summaries.append(het_summary)
                skipped_count += 1
                logger.info(f"Skipped completed locus: {locus_id}")
            else:
                pending_args.append(
                    (
                        locus_id,
                        locus_info_group,
                        outdir,
                        meta_method,
                        calculate_lambda_s,
                    )
                )
    else:
        pending_args = [
            (locus_id, locus_info_group, outdir, meta_method, calculate_lambda_s)
            for locus_id, locus_info_group in grouped_loci
        ]

    if skipped_count > 0:
        logger.info(
            f"Recovered {skipped_count}/{total_loci} loci, "
            f"processing {len(pending_args)} remaining"
        )

    # Create process pool and process loci in parallel with progress bar
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task(
            "[cyan]Meta-analysing...", total=total_loci, completed=skipped_count
        )

        if pending_args:
            with Pool(threads) as pool:
                for result, het_summary in pool.imap_unordered(meta_locus, pending_args):  # type: ignore
                    for i, res in enumerate(result):
                        new_loci_info.loc[len(new_loci_info)] = res
                    if het_summary is not None and not het_summary.empty:
                        all_het_summaries.append(het_summary)
                    progress.advance(task)

    new_loci_info.to_csv(f"{outdir}/loci_info.txt", sep="\t", index=False)

    # Save global heterogeneity summary
    if all_het_summaries:
        global_het = pd.concat(all_het_summaries, ignore_index=True)
        cols = ["locus_id"] + [c for c in global_het.columns if c != "locus_id"]
        global_het = global_het[cols]
        global_het.to_csv(
            f"{outdir}/heterogeneity.txt.gz",
            sep="\t",
            index=False,
            float_format="%.6f",
            compression="gzip",
        )
