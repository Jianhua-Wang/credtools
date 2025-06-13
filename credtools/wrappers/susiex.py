"""Wrapper for SuSiEx multi-ancestry fine-mapping."""

import json
import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import Locus, LocusSet, intersect_sumstat_ld
from credtools.utils import io_in_tempdir, tool_manager

logger = logging.getLogger("SuSiEx")


@io_in_tempdir("./tmp/SuSiEx")
def run_susiex(
    locus_set: LocusSet,
    max_causal: int = 1,
    coverage: float = 0.95,
    pval_thresh: float = 1,
    maf_thresh: float = 0.005,
    mult_step: bool = True,
    keep_ambig: bool = True,
    n_threads: int = 1,
    min_purity: float = 0,
    max_iter: int = 100,
    tol: float = 1e-3,
    temp_dir: Optional[str] = None,
) -> CredibleSet:
    """
    Run SuSiEx multi-ancestry fine-mapping analysis on a LocusSet.

    SuSiEx (SuSiE for Cross-ancestry) extends the SuSiE framework for
    multi-ancestry fine-mapping, allowing joint analysis of GWAS summary
    statistics from multiple populations. It accounts for population-specific
    LD structures while identifying shared causal variants across ancestries.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple Locus objects from different
        populations/ancestries covering the same genomic region.
        Each locus should have summary statistics and LD matrix data.
    max_causal : int, optional
        Maximum number of causal variants to detect, by default 1.
        This sets the upper limit on the number of independent signals
        that can be identified in the analysis.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
        Determines the cumulative posterior probability mass required
        for credible set inclusion.
    pval_thresh : float, optional
        P-value threshold for variant inclusion, by default 1.
        Variants with p-values above this threshold are excluded
        from the analysis to reduce computational burden.
    maf_thresh : float, optional
        Minor allele frequency threshold for variant filtering, by default 0.005.
        Variants with MAF below this threshold in any population are
        excluded to avoid spurious associations from rare variants.
    mult_step : bool, optional
        Whether to use multi-step refinement procedure, by default True.
        Multi-step refinement can improve fine-mapping resolution by
        iteratively updating the analysis with refined credible sets.
    keep_ambig : bool, optional
        Whether to retain ambiguous variants in the analysis, by default True.
        Ambiguous variants are those with unclear strand orientation
        or allele assignments across populations.
    n_threads : int, optional
        Number of parallel threads for computation, by default 1.
        Higher values can speed up analysis but require more memory.
    min_purity : float, optional
        Minimum purity threshold for credible sets, by default 0.
        Credible sets with purity below this threshold may be filtered
        out as potentially unreliable.
    max_iter : int, optional
        Maximum number of iterations for the optimization algorithm, by default 100.
        More iterations may improve convergence but increase runtime.
    tol : float, optional
        Convergence tolerance for the optimization algorithm, by default 1e-3.
        Smaller values require tighter convergence but may improve accuracy.
    temp_dir : Optional[str], optional
        Temporary directory for intermediate files, by default None.
        If None, a temporary directory is automatically created.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Multi-ancestry posterior inclusion probabilities
        - Credible sets identified across populations
        - Lead SNPs for each credible set
        - Coverage and purity statistics
        - SuSiEx-specific parameters and results

    Notes
    -----
    SuSiEx implements a multi-ancestry extension of the SuSiE model:

    y_k = Σ(l=1 to L) X_k * γ_l * β_l + ε_k

    where:
    - y_k is the phenotype vector for population k
    - X_k is the genotype matrix for population k
    - γ_l is the binary indicator vector for causal variants in component l
    - β_l is the effect size vector for component l
    - ε_k is the residual error for population k

    Key features:

    1. **Cross-ancestry information sharing**: Leverages evidence from
       multiple populations to improve fine-mapping power and resolution

    2. **Population-specific LD modeling**: Accounts for different LD
       patterns across populations while assuming shared causal variants

    3. **Adaptive variant filtering**: Applies population-specific quality
       control filters while maintaining variant overlap for joint analysis

    4. **Multi-step refinement**: Optionally refines results through
       iterative analysis with updated variant sets

    The algorithm workflow:
    1. Harmonize summary statistics and LD data across populations
    2. Apply quality control filters (MAF, p-value thresholds)
    3. Convert data to SuSiEx input format
    4. Run multi-ancestry SuSiE analysis
    5. Extract credible sets with cross-population evidence
    6. Apply purity filters and post-processing

    File format requirements:
    - Summary statistics: tab-separated with standard column names
    - LD matrices: binary format with corresponding variant maps
    - Population-specific allele frequency data

    Advantages over single-population methods:
    - Increased statistical power through meta-analysis
    - Better fine-mapping resolution via diverse LD patterns
    - Robustness to population-specific confounding
    - Natural framework for trans-ethnic studies

    The method automatically handles:
    - Variant harmonization across populations
    - Population-specific quality control
    - LD matrix format conversion
    - Multi-threading for computational efficiency

    Examples
    --------
    >>> # Basic multi-ancestry fine-mapping
    >>> credible_set = run_susiex(locus_set)
    >>> print(f"Found {credible_set.n_cs} credible sets")
    >>> print(f"Populations: {len(locus_set.loci)}")
    Found 1 credible sets
    Populations: 3

    >>> # SuSiEx with multiple signals and strict quality control
    >>> credible_set = run_susiex(
    ...     locus_set,
    ...     max_causal=5,
    ...     coverage=0.99,
    ...     maf_thresh=0.01,    # Stricter MAF filter
    ...     min_purity=0.5,     # Require high purity
    ...     mult_step=True      # Multi-step refinement
    ... )
    >>> print(f"Detected {credible_set.n_cs} high-confidence signals")
    >>> print(f"Credible set sizes: {credible_set.cs_sizes}")
    Detected 2 high-confidence signals
    Credible set sizes: [8, 12]

    >>> # High-performance analysis with parallel processing
    >>> credible_set = run_susiex(
    ...     locus_set,
    ...     n_threads=8,        # Parallel processing
    ...     max_iter=200,       # More iterations
    ...     tol=1e-6           # Tight convergence
    ... )
    >>> print(f"Analysis completed with {credible_set.n_cs} signals")
    Analysis completed with 3 signals

    >>> # Access cross-population posterior inclusion probabilities
    >>> top_variants = credible_set.pips.nlargest(10)
    >>> print("Top 10 variants by cross-population PIP:")
    >>> print(top_variants)
    Top 10 variants by cross-population PIP:
    rs123456    0.8934
    rs789012    0.7234
    rs345678    0.6456
    ...

    >>> # Examine credible sets from multi-ancestry analysis
    >>> for i, snps in enumerate(credible_set.snps):
    ...     lead_snp = credible_set.lead_snps[i]
    ...     pip = credible_set.pips[lead_snp]
    ...     size = len(snps)
    ...     print(f"Signal {i+1}: {lead_snp} (PIP={pip:.4f}, size={size})")
    Signal 1: rs123456 (PIP=0.8934, size=8)
    Signal 2: rs789012 (PIP=0.7234, size=12)
    Signal 3: rs345678 (PIP=0.6456, size=15)
    """
    logger.info(f"Running SuSiEx on {locus_set}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "pval_thresh": pval_thresh,
        "maf_thresh": maf_thresh,
        "mult_step": mult_step,
        "keep_ambig": keep_ambig,
        "n_threads": n_threads,
        "min_purity": min_purity,
        "max_iter": max_iter,
        "tol": tol,
    }
    logger.info(f"Parameters: {parameters}")

    input_prefix_list = []
    for locus in locus_set.loci:
        locus = intersect_sumstat_ld(locus)
        input_prefix = f"{temp_dir}/{locus.popu}.{locus.cohort}"
        logger.debug(f"Writing {input_prefix}.sumstats")
        locus.sumstats.to_csv(f"{input_prefix}.sumstats", sep="\t", index=False)
        ldmap = locus.ld.map.copy()
        ldmap["cm"] = 0
        logger.debug(f"Writing {input_prefix}_ref.bim")
        ldmap[
            [ColName.CHR, ColName.SNPID, "cm", ColName.BP, ColName.A1, ColName.A2]
        ].to_csv(f"{input_prefix}_ref.bim", sep="\t", index=False, header=False)
        ldmap["MAF"] = ldmap["AF2"].apply(lambda x: min(x, 1 - x))
        ldmap["NCHROBS"] = locus.sample_size * 2
        ldmap.rename(columns={"SNPID": "SNP"}, inplace=True)
        logger.debug(f"Writing {input_prefix}_frq.frq")
        ldmap[["CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"]].to_csv(
            f"{input_prefix}_frq.frq", sep="\t", index=False
        )
        logger.debug(f"Writing {input_prefix}.ld.bin")
        ld = locus.ld.r**2
        ld.astype(np.float32).tofile(f"{input_prefix}.ld.bin")
        input_prefix_list.append(input_prefix)

    sst_file = ",".join([i + ".sumstats" for i in input_prefix_list])
    n_gwas = ",".join([str(locus.sample_size) for locus in locus_set.loci])
    ld_file = ",".join(input_prefix_list)
    chrom, start, end = locus_set.chrom, locus_set.start, locus_set.end
    cmd = [
        f"--sst_file={sst_file}",
        f"--n_gwas={n_gwas}",
        f"--ld_file={ld_file}",
        f"--out_dir={temp_dir}",
        f"--out_name=chr{chrom}_{start}_{end}",
        f"--level={coverage}",
        f"--pval_thresh={pval_thresh}",
        f"--maf={maf_thresh}",
        f"--chr={chrom}",
        f"--bp={start},{end}",
        f'--snp_col={",".join(["1"]*locus_set.n_loci)}',
        f'--chr_col={",".join(["2"]*locus_set.n_loci)}',
        f'--bp_col={",".join(["3"]*locus_set.n_loci)}',
        f'--a1_col={",".join(["5"]*locus_set.n_loci)}',
        f'--a2_col={",".join(["6"]*locus_set.n_loci)}',
        f'--eff_col={",".join(["9"]*locus_set.n_loci)}',
        f'--se_col={",".join(["10"]*locus_set.n_loci)}',
        f'--pval_col={",".join(["11"]*locus_set.n_loci)}',
        "--plink=../utilities/plink",
        f"--n_sig={max_causal}",
        f"--mult-step={mult_step}",
        f"--keep-ambig={keep_ambig}",
        f"--threads={n_threads}",
        f"--min_purity={min_purity}",
        f"--max_iter={max_iter}",
        f"--tol={tol}",
    ]
    required_output_files = [
        f"{temp_dir}/chr{chrom}_{start}_{end}.snp",
        f"{temp_dir}/chr{chrom}_{start}_{end}.cs",
    ]
    logger.info(f"Running SuSiEx with command: {' '.join(cmd)}.")
    tool_manager.run_tool("SuSiEx", cmd, f"{temp_dir}/run.log", required_output_files)

    pip_df = pd.read_csv(f"{temp_dir}/chr{chrom}_{start}_{end}.snp", sep="\t")
    cs_snp: List[List[str]] = []
    if len(pip_df.columns) == 2:
        logger.warning("No credible set found, please try other parameters.")
        pip = pd.Series(index=pip_df["SNP"].values.tolist())
    else:
        cs_df = pd.read_csv(f"{temp_dir}/chr{chrom}_{start}_{end}.cs", sep="\t")
        for _, sub_df in cs_df.groupby("CS_ID"):
            cs_snp.append(sub_df["SNP"].values.tolist())
        pip_cols = [col for col in pip_df.columns if col.startswith("PIP")]
        pip_df = pip_df[pip_cols + ["SNP"]].copy()
        pip_df.set_index("SNP", inplace=True)
        pip_df["PIP"] = pip_df[pip_cols].max(axis=1)
        pip = pd.Series(
            index=pip_df.index.values.tolist(), data=pip_df["PIP"].values.tolist()
        )

    logger.info(f"Finished SuSiEx on {locus_set}")
    logger.info(f"N of credible set: {len(cs_snp)}")
    logger.info(f"Credible set size: {[len(i) for i in cs_snp]}")
    return CredibleSet(
        tool=Method.SUSIEX,
        n_cs=len(cs_snp),
        coverage=coverage,
        lead_snps=[str(pip[pip.index.isin(i)].idxmax()) for i in cs_snp],
        snps=cs_snp,
        cs_sizes=[len(i) for i in cs_snp],
        pips=pip,
        parameters=parameters,
    )
