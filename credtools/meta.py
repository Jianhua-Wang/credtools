"""Meta analysis of multi-ancestry gwas data."""

import json
import logging
import os
from multiprocessing import Pool
from typing import Any, Dict, List, Optional, Tuple, Union

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
from credtools.meta_weights import matched_meta, validate_ld_weighting

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
    This standalone function does not match LD. For paired fine-mapping
    statistics and LD, use meta_all(), which shares contribution eligibility.

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


def meta_lds(inputs: LocusSet, ld_weighting: str = "ess") -> LDMatrix:
    """Merge within-cohort GWAS/LD matches using geometric ESS or SE weights.

    ESS (default) uses the supplied scalar sample_size; SE uses 1/SE**2.
    Normalize at each SNP, then multiply the square-root coefficients at both
    ends of each LD correlation. Missing contributions are zero, not pairwise
    renormalized. Use meta_all() to obtain the paired IVW summary statistics.
    """
    return matched_meta(inputs, ld_weighting)[1]


def meta_all(inputs: LocusSet, ld_weighting: str = "ess") -> Locus:
    """
    Perform comprehensive meta-analysis of both summary statistics and LD matrices.

    Parameters
    ----------
    inputs : LocusSet
        LocusSet containing input data from multiple studies.
    ld_weighting : {"ess", "se"}, optional
        Geometric LD weights after within-cohort matching, by default "ess".
        Summary statistics always use matched inverse-variance weighting.

    Returns
    -------
    Locus
        Meta-analyzed Locus object with combined population and cohort identifiers.

    Notes
    -----
    This function:

    1. Performs meta-analysis of summary statistics using inverse-variance weighting
    2. Merges LD using geometric ESS (default) or SE weights on the same contributions
    3. Combines population and cohort names from all input studies
    4. Sums sample sizes across studies
    5. Intersects the meta-analyzed data to ensure consistency

    Population and cohort names are combined with "+" as separator and sorted alphabetically.
    """
    validate_ld_weighting(ld_weighting)
    if not inputs.loci:
        raise ValueError("No input cohorts for meta-analysis")
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

    meta_sumstat, meta_ld, audit, cohort_audit = matched_meta(inputs, ld_weighting)
    result = Locus(
        popu,
        cohort,
        sample_size,
        meta_sumstat,
        locus_starts[0],
        locus_ends[0],
        ld=meta_ld,
        if_intersect=True,
    )
    result.meta_audit = audit
    result.meta_cohort_audit = cohort_audit
    return result


def meta_by_population(inputs: LocusSet, ld_weighting: str = "ess") -> Dict[str, Locus]:
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
    3. For single-study populations, applies the same matching and validation
    4. Returns a dictionary with population codes as keys

    This approach preserves population-specific LD patterns while still
    allowing meta-analysis of multiple cohorts within the same population.
    """
    validate_ld_weighting(ld_weighting)
    meta_popu = {}
    for input in inputs.loci:
        popu = input.popu
        if popu not in meta_popu:
            meta_popu[popu] = [input]
        else:
            meta_popu[popu].append(input)

    result_dict = {}
    for popu in meta_popu:
        result_dict[popu] = meta_all(LocusSet(meta_popu[popu]), ld_weighting)
    return result_dict


def meta(inputs: LocusSet, meta_method: str = "meta_all", ld_weighting: str = "ess") -> LocusSet:
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
    validate_ld_weighting(ld_weighting)
    if meta_method == "meta_all":
        return LocusSet([meta_all(inputs, ld_weighting)])
    elif meta_method == "meta_by_population":
        res = meta_by_population(inputs, ld_weighting)
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


def meta_configuration(meta_method: str, ld_weighting: str) -> Dict[str, Any]:
    """Describe the input contract used to protect output reuse."""
    validate_ld_weighting(ld_weighting)
    if meta_method not in ("meta_all", "meta_by_population", "no_meta"):
        raise ValueError(f"Unsupported meta-analysis method: {meta_method}")
    return {
        "schema_version": 1,
        "algorithm": "cohort_matched_geometric_v1" if meta_method != "no_meta" else "no_meta_v1",
        "meta_method": meta_method,
        "ld_weighting": ld_weighting if meta_method != "no_meta" else None,
        "summary_weighting": "inverse_variance" if meta_method != "no_meta" else None,
    }


def ensure_meta_configuration(outdir: str, meta_method: str, ld_weighting: str) -> None:
    """Refuse legacy/mixed-mode outputs; write a manifest for new output directories."""
    expected = meta_configuration(meta_method, ld_weighting)
    path = os.path.join(outdir, "meta_config.json")
    os.makedirs(outdir, exist_ok=True)
    if os.path.exists(path):
        with open(path) as handle:
            existing = json.load(handle)
        if existing != expected:
            raise ValueError("Meta configuration differs from existing output; use a separate output directory")
        return
    # A pre-feature result has no evidence of matching or its LD weighting.
    # Do not silently reuse or overwrite it, even without --skip.
    if any(name == "loci_info.txt" or name.endswith((".ld.npz", ".sumstat", ".sumstats.gz"))
           for name in os.listdir(outdir)):
        raise ValueError("Existing meta output has no configuration manifest; use a new output directory")
    with open(path, "w") as handle:
        json.dump(expected, handle, indent=2)
        handle.write("\n")


def meta_output_prefix(locus: Locus, meta_method: str, ld_weighting: str) -> str:
    """Make matched LD weighting explicit in output filenames."""
    suffix = "no_meta" if meta_method == "no_meta" else f"matched_{ld_weighting}"
    return f"{locus.prefix}.{suffix}"


def save_meta_audit(locus: Locus, out_prefix: str) -> None:
    """Save SNP information retention and cohort contribution counts when available."""
    for attr, suffix in (("meta_audit", "meta_audit"), ("meta_cohort_audit", "meta_cohorts")):
        audit = getattr(locus, attr, None)
        if audit is not None:
            audit.to_csv(f"{out_prefix}.{suffix}.tsv.gz", sep="\t", index=False, compression="gzip")


def recover_completed_locus(
    locus_id: str,
    outdir: str,
    prev_loci_info: Optional[pd.DataFrame],
    configuration: Optional[Dict[str, Any]] = None,
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
    if configuration is not None:
        try:
            with open(os.path.join(locus_dir, "meta_config.json")) as handle:
                if json.load(handle) != configuration:
                    return None
        except (OSError, ValueError):
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
    args: Union[Tuple[str, pd.DataFrame, str, str, bool], Tuple[str, pd.DataFrame, str, str, bool, str]],
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
        - ld_weighting : str, optional sixth item
            Geometric LD weighting, "ess" (default) or "se".

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
    # Retain compatibility with the previous five-item worker argument tuple.
    locus_id, locus_info, outdir, meta_method, calculate_lambda_s = args[:5]
    ld_weighting = args[5] if len(args) > 5 else "ess"
    out_dir = os.path.abspath(f"{outdir}/{locus_id}")
    ensure_meta_configuration(out_dir, meta_method, ld_weighting)
    results = []
    locus_set = load_locus_set(locus_info, calculate_lambda_s=calculate_lambda_s)

    # Compute heterogeneity BEFORE meta combines data
    if meta_method == "meta_by_population":
        het_metrics = compute_heterogeneity_by_population(locus_set)
    else:
        het_metrics = compute_heterogeneity(locus_set)
    het_summary = heterogeneity_summary(het_metrics, locus_set)
    het_summary["locus_id"] = locus_id

    locus_set = meta(locus_set, meta_method, ld_weighting)
    out_dir = f"{outdir}/{locus_id}"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Save heterogeneity alongside meta results
    save_heterogeneity(het_metrics, out_dir, summary=het_summary)

    for locus in locus_set.loci:
        out_prefix = f"{out_dir}/{meta_output_prefix(locus, meta_method, ld_weighting)}"
        save_meta_audit(locus, out_prefix)
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
    ld_weighting: str = "ess",
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
    ld_weighting : {"ess", "se"}, optional
        Cohort-matched geometric LD weighting, by default "ess".

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
    configuration = meta_configuration(meta_method, ld_weighting)
    loci_info = pd.read_csv(inputs, sep="\t")
    loci_info = check_loci_info(loci_info)  # Validate input data
    # Worker records have this fixed schema; input TSV column order is arbitrary.
    # Do not relabel positional results using the input's column order.
    new_loci_info = pd.DataFrame(columns=[
        "chr", "start", "end", "popu", "sample_size", "cohort", "prefix", "locus_id"
    ])
    all_het_summaries: List[pd.DataFrame] = []

    # Group loci by locus_id
    grouped_loci = list(loci_info.groupby("locus_id"))
    total_loci = len(grouped_loci)
    ensure_meta_configuration(outdir, meta_method, ld_weighting)

    # Try to recover completed loci when skip=True
    prev_loci_info: Optional[pd.DataFrame] = None
    pending_args: List[Tuple[str, pd.DataFrame, str, str, bool, str]] = []
    skipped_count = 0

    if skip:
        prev_path = os.path.join(outdir, "loci_info.txt")
        if os.path.exists(prev_path):
            prev_loci_info = pd.read_csv(prev_path, sep="\t")
            logger.info(f"Found previous loci_info.txt with {len(prev_loci_info)} rows")

    if prev_loci_info is not None:
        for locus_id, locus_info_group in grouped_loci:
            recovered = recover_completed_locus(locus_id, outdir, prev_loci_info, configuration)
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
                        ld_weighting,
                    )
                )
    else:
        pending_args = [
            (locus_id, locus_info_group, outdir, meta_method, calculate_lambda_s, ld_weighting)
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
