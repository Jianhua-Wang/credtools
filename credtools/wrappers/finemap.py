"""Wrapper for FINEMAP."""

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import Locus, intersect_sumstat_ld
from credtools.utils import io_in_tempdir, tool_manager

logger = logging.getLogger("FINEMAP")


@io_in_tempdir("./tmp/FINEMAP")
def run_finemap(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    n_iter: int = 100000,
    n_threads: int = 1,
    temp_dir: Optional[str] = None,
) -> CredibleSet:
    """
    Run FINEMAP with shotgun stochastic search for fine-mapping analysis.

    FINEMAP uses a shotgun stochastic search algorithm to efficiently explore
    the space of possible causal configurations and identify credible sets
    of variants that jointly explain the observed association signal.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix data.
        Must include MAF (minor allele frequency) information.
    max_causal : int, optional
        Maximum number of causal variants to consider, by default 1.
        Higher values allow for more complex causal models but increase
        computational time exponentially.
    coverage : float, optional
        Coverage probability of the credible set, by default 0.95.
        This determines the cumulative posterior probability mass
        included in the credible set.
    n_iter : int, optional
        Number of iterations for the stochastic search, by default 100000.
        More iterations generally improve accuracy but increase runtime.
    n_threads : int, optional
        Number of threads to use for parallel computation, by default 1.
        Multiple threads can significantly speed up analysis for large datasets.
    temp_dir : Optional[str], optional
        Temporary directory for intermediate files, by default None.
        Automatically provided by the @io_in_tempdir decorator.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Number of causal variants identified
        - Posterior inclusion probabilities for all variants
        - Credible set membership for each causal configuration
        - Lead SNPs (most significant in each credible set)
        - Coverage probability and method parameters

    Raises
    ------
    ValueError
        If MAF column is not present in the locus summary statistics.

    Warnings
    --------
    If the summary statistics and LD matrix are not matched (different variant orders),
    they will be automatically intersected and reordered.

    Notes
    -----
    FINEMAP implements a shotgun stochastic search algorithm that:

    1. Samples causal configurations from the posterior distribution
    2. Uses Bayes factors to evaluate each configuration
    3. Accounts for linkage disequilibrium through the LD matrix
    4. Provides posterior probabilities for different numbers of causal variants
    5. Constructs credible sets based on posterior inclusion probabilities

    The algorithm can handle multiple causal variants and provides:
    - Posterior probabilities for 0, 1, 2, ... max_causal variants
    - Configuration-specific credible sets
    - Model selection based on posterior probabilities

    Input file requirements:
    - Summary statistics with SNPID, CHR, BP, EA, NEA, MAF, BETA, SE
    - LD matrix in lower-triangular format
    - Master file linking all input components

    Output interpretation:
    - Higher posterior inclusion probabilities indicate stronger evidence for causality
    - Multiple credible sets may be identified for different causal models
    - The optimal number of causal variants is selected based on posterior probabilities

    Reference:
    Benner, C. et al. FINEMAP: efficient variable selection using summary data from
    genome-wide association studies. Bioinformatics 32, 1493-1501 (2016).

    Examples
    --------
    >>> # Basic FINEMAP analysis
    >>> credible_set = run_finemap(locus)
    >>> print(f"Found {credible_set.n_cs} causal variants")
    >>> print(f"Credible set sizes: {credible_set.cs_sizes}")
    Found 2 causal variants
    Credible set sizes: [5, 8]

    >>> # FINEMAP with multiple causal variants and high coverage
    >>> credible_set = run_finemap(
    ...     locus,
    ...     max_causal=3,
    ...     coverage=0.99,
    ...     n_iter=500000
    ... )
    >>> print(f"Coverage: {credible_set.coverage}")
    >>> print(f"Lead SNPs: {credible_set.lead_snps}")
    Coverage: 0.99
    Lead SNPs: ['rs123456', 'rs789012', 'rs345678']

    >>> # Access posterior inclusion probabilities
    >>> top_pips = credible_set.pips.nlargest(10)
    >>> print("Top 10 variants by PIP:")
    >>> print(top_pips)
    Top 10 variants by PIP:
    rs123456    0.85
    rs789012    0.72
    rs345678    0.68
    ...
    """
    logger.info(f"Running FINEMAP on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "n_iter": n_iter,
        "n_threads": n_threads,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    if not locus.is_matched:
        logger.warning(
            "The sumstat and LD are not matched, will match them in same order."
        )
        locus = intersect_sumstat_ld(locus)
    # write z file
    if ColName.MAF not in locus.sumstats.columns:
        raise ValueError(f"{ColName.MAF} is required for FINEMAP.")
    finemap_input = locus.sumstats.copy()
    finemap_input[ColName.MAF] = finemap_input[ColName.MAF].replace(0, 0.00001)
    finemap_input = finemap_input[
        [
            ColName.SNPID,
            ColName.CHR,
            ColName.BP,
            ColName.EA,
            ColName.NEA,
            ColName.MAF,
            ColName.BETA,
            ColName.SE,
        ]
    ]
    finemap_input.rename(
        columns={
            ColName.SNPID: "rsid",
            ColName.CHR: "chromosome",
            ColName.BP: "position",
            ColName.MAF: "maf",
            ColName.BETA: "beta",
            ColName.SE: "se",
            ColName.EA: "allele1",
            ColName.NEA: "allele2",
        },
        inplace=True,
    )
    # change maf to 0.00001 if maf is 0
    finemap_input.loc[finemap_input["maf"] <= 0.00001, "maf"] = 0.00001
    logger.info(f"Writing FINEMAP input to {temp_dir}/finemap.z")
    finemap_input.to_csv(
        f"{temp_dir}/finemap.z", sep=" ", index=False, float_format="%0.5f"
    )

    # write ld file
    logger.info(f"Writing FINEMAP LD file to {temp_dir}/finemap.ld")
    np.savetxt(f"{temp_dir}/finemap.ld", locus.ld.r, delimiter=" ", fmt="%0.4f")
    # TODO: write ld file only once for multiple tools
    # TODO: use BCOR file for LD

    # write master file
    logger.info(f"Writing FINEMAP master file to {temp_dir}/finemap.master")
    with open(f"{temp_dir}/finemap.master", "w") as f:
        master_content = [
            f"{temp_dir}/finemap.z",
            f"{temp_dir}/finemap.ld",
            f"{temp_dir}/finemap.snp",
            f"{temp_dir}/finemap.config",
            f"{temp_dir}/finemap.cred",
            f"{temp_dir}/finemap.log",
            str(locus.sample_size),
        ]
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(";".join(master_content))

    # run finemap
    cmd = [
        "--sss",
        "--in-files",
        f"{temp_dir}/finemap.master",
        "--n-causal-snps",
        str(max_causal),
        "--n-iter",
        str(n_iter),
        "--n-threads",
        str(n_threads),
        "--prob-cred-set",
        str(coverage),
    ]
    required_output_files = [f"{temp_dir}/finemap.snp", f"{temp_dir}/finemap.config"]
    logger.info(f"Running FINEMAP with command: {' '.join(cmd)}.")
    tool_manager.run_tool("finemap", cmd, f"{temp_dir}/run.log", required_output_files)

    # get credible set
    cred_file_list = Path(f"{temp_dir}/").glob("finemap.cred*")
    cred_prob = {}
    pip = pd.Series(index=finemap_input["rsid"].values.tolist(), data=0.0)
    n_cs = 0
    cs_snps = []
    lead_snps = []
    cs_sizes = []
    for cred_file in cred_file_list:
        with open(cred_file, "r") as f:
            n_causal = int(cred_file.name[-1])
            first_line = f.readline()
            cred_prob[n_causal] = float(first_line.split()[-1])
    if len(cred_prob) == 0:
        logger.warning("FINEMAP output is empty.")
    else:
        # get the credible set with the highest posterior probability
        n_cs = max(cred_prob, key=lambda k: float(cred_prob[k]))
        logger.info(
            f"FINEMAP found {n_cs} causal SNPs with the post-prob {cred_prob[n_cs]}."
        )
        cred_set = pd.read_csv(f"{temp_dir}/finemap.cred{n_cs}", sep=" ", comment="#")
        for cred_idx in range(1, n_cs + 1):
            cred_df = cred_set[[f"cred{cred_idx}", f"prob{cred_idx}"]].copy()
            cred_df.rename(
                columns={f"cred{cred_idx}": "snp", f"prob{cred_idx}": "pip"},
                inplace=True,
            )
            cred_df.dropna(inplace=True)
            cs_snps.append(cred_df["snp"].values.tolist())
            cs_sizes.append(len(cred_df["snp"].values.tolist()))
            pip[cred_df["snp"].values.tolist()] = cred_df["pip"].values.tolist()
            lead_snps.append(
                str(
                    locus.sumstats.loc[
                        locus.sumstats[
                            locus.sumstats[ColName.SNPID].isin(
                                cred_df["snp"].values.tolist()
                            )
                        ][ColName.P].idxmin(),
                        ColName.SNPID,
                    ]
                )
            )

    # output
    logger.info(f"Finished FINEMAP on {locus}")
    logger.info(f"N of credible set: {n_cs}")
    logger.info(f"Credible set size: {cs_sizes}")
    return CredibleSet(
        tool=Method.FINEMAP,
        n_cs=n_cs,
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=cs_sizes,
        pips=pip,
        parameters=parameters,
    )
