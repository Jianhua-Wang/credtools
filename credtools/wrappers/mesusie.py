"""Wrapper for MESuSiE multi-ancestry fine-mapping."""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import LocusSet, intersect_sumstat_ld
from credtools.utils import io_in_tempdir

logger = logging.getLogger("MESUSIE")

# Path to the R wrapper script bundled with this package
_R_SCRIPT_PATH = str(Path(__file__).parent / "mesusie_wrapper.R")


def _check_r_and_mesusie():
    """Check that Rscript and the MESuSiE R package are available.

    Raises
    ------
    FileNotFoundError
        If Rscript is not on PATH or MESuSiE is not installed.
    """
    if shutil.which("Rscript") is None:
        raise FileNotFoundError(
            "Rscript not found on PATH. Please install R "
            "(https://cran.r-project.org/) and ensure Rscript is available."
        )

    result = subprocess.run(
        ["Rscript", "-e", "library(MESuSiE)"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise FileNotFoundError(
            "MESuSiE R package is not installed. Please install it with:\n"
            "  R -e 'devtools::install_github(\"borangao/MESuSiE\")'"
        )


@io_in_tempdir("./tmp/MESuSiE")
def run_mesusie(
    locus_set: LocusSet,
    max_causal: int = 1,
    coverage: float = 0.95,
    max_iter: int = 100,
    tol: float = 1e-3,
    purity: float = 0.0,
    estimate_residual_variance: bool = False,
    empty_on_nonconvergence: bool = False,
    temp_dir: Optional[str] = None,
) -> CredibleSet:
    """
    Run MESuSiE multi-ancestry fine-mapping analysis on a LocusSet.

    MESuSiE (Multiple Ancestry Sum of Single Effects) extends the SuSiE
    framework for multi-ancestry fine-mapping, identifying both shared and
    ancestry-specific causal variants across populations.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple Locus objects from different
        populations/ancestries covering the same genomic region.
    max_causal : int, optional
        Maximum number of causal variants (L parameter), by default 1.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
    max_iter : int, optional
        Maximum number of iterations, by default 100.
    tol : float, optional
        Convergence tolerance, by default 1e-3.
    purity : float, optional
        Minimum purity threshold for credible sets, by default 0.0.
    estimate_residual_variance : bool, optional
        Whether to estimate residual variance, by default False.
    empty_on_nonconvergence : bool, optional
        When True and the IBSS algorithm did not converge, return an empty
        credible set (``n_cs=0``, all PIPs zeroed) instead of the unreliable
        partially-fit results. Defaults to False (preserve legacy behavior).
    temp_dir : Optional[str], optional
        Temporary directory for intermediate files, by default None.

    Returns
    -------
    CredibleSet
        Credible set object containing PIPs, credible sets, and MESuSiE-specific
        parameters including cs_types (shared/ancestry-specific classification).

    Raises
    ------
    FileNotFoundError
        If Rscript or MESuSiE R package is not available.
    RuntimeError
        If the R script execution fails.
    """
    logger.info(f"Running MESuSiE on {locus_set}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "max_iter": max_iter,
        "tol": tol,
        "purity": purity,
        "estimate_residual_variance": estimate_residual_variance,
        "empty_on_nonconvergence": empty_on_nonconvergence,
    }
    logger.info(f"Parameters: {parameters}")

    # Check R and MESuSiE availability
    _check_r_and_mesusie()

    # Step 1: Intersect sumstats and LD within each population
    loci = []
    for locus in locus_set.loci:
        loci.append(intersect_sumstat_ld(locus))

    # Step 2: Find common SNPs across all populations
    common_snps = None
    for locus in loci:
        snp_set = set(locus.sumstats[ColName.SNPID].values)
        if common_snps is None:
            common_snps = snp_set
        else:
            common_snps = common_snps.intersection(snp_set)

    common_snps_list = sorted(common_snps)  # Sort for deterministic order
    logger.info(
        f"Found {len(common_snps_list)} common SNPs across {len(loci)} populations"
    )

    if len(common_snps_list) == 0:
        raise ValueError("No common SNPs found across populations")

    # Step 3: Write input files for each population (subset to common SNPs)
    for i, locus in enumerate(loci):
        sumstat = locus.sumstats.copy()
        sumstat = sumstat[sumstat[ColName.SNPID].isin(common_snps)].copy()
        # Sort by common_snps_list order for consistency
        sumstat = sumstat.set_index(ColName.SNPID).loc[common_snps_list].reset_index()

        sumstat[ColName.Z] = sumstat[ColName.BETA] / sumstat[ColName.SE]
        ss_out = pd.DataFrame(
            {
                "SNP": sumstat[ColName.SNPID].values,
                "Beta": sumstat[ColName.BETA].values,
                "Se": sumstat[ColName.SE].values,
                "Z": sumstat[ColName.Z].values,
                "N": [locus.sample_size] * len(sumstat),
            }
        )
        ss_out.to_csv(f"{temp_dir}/pop_{i}_sumstats.csv", index=False)

        # Write SNP IDs
        snpids = sumstat[ColName.SNPID].values.tolist()
        with open(f"{temp_dir}/pop_{i}_snpids.txt", "w") as f:
            f.write("\n".join(snpids) + "\n")

        # Subset LD matrix to common SNPs
        ldmap = locus.ld.map.copy()
        ld_r = locus.ld.r.copy()
        # Get indices of common SNPs in LD map
        ldmap_snps = ldmap[ColName.SNPID].values.tolist()
        common_indices = [
            ldmap_snps.index(s) for s in common_snps_list if s in ldmap_snps
        ]
        ld_subset = ld_r[np.ix_(common_indices, common_indices)]

        # Write LD matrix as float64 binary
        ld_subset.astype(np.float64).tofile(f"{temp_dir}/pop_{i}_ld.bin")

        # Write LD matrix dimension
        with open(f"{temp_dir}/pop_{i}_ld_dim.txt", "w") as f:
            f.write(str(ld_subset.shape[0]) + "\n")

        logger.debug(
            f"Wrote input files for population {i} ({locus.popu}): {len(snpids)} SNPs"
        )

    # Write population names file
    pop_names = [locus.popu for locus in locus_set.loci]
    with open(f"{temp_dir}/pop_names.txt", "w") as f:
        f.write("\n".join(pop_names) + "\n")

    # Build Rscript command
    cmd = [
        "Rscript",
        _R_SCRIPT_PATH,
        "--temp_dir",
        temp_dir,
        "--n_pop",
        str(locus_set.n_loci),
        "--L",
        str(max_causal),
        "--coverage",
        str(coverage),
        "--max_iter",
        str(max_iter),
        "--purity",
        str(purity),
        "--estimate_residual_variance",
        str(estimate_residual_variance).upper(),
    ]

    logger.info(f"Running MESuSiE R script: {' '.join(cmd)}")

    # Run R script
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"MESuSiE R script stderr: {result.stderr}")
        logger.error(f"MESuSiE R script stdout: {result.stdout}")
        raise RuntimeError(
            f"MESuSiE R script failed (return code {result.returncode}).\n"
            f"stderr: {result.stderr}\n"
            f"stdout: {result.stdout}"
        )

    logger.debug(f"MESuSiE R script stdout: {result.stdout}")

    # Parse output files
    pip_df = pd.read_csv(f"{temp_dir}/mesusie_pips.csv")
    cs_df = pd.read_csv(f"{temp_dir}/mesusie_cs.csv")
    purity_df = pd.read_csv(f"{temp_dir}/mesusie_purity.csv")

    # Read convergence status
    converged_file = f"{temp_dir}/mesusie_converged.txt"
    converged: Optional[bool] = None
    if os.path.exists(converged_file):
        with open(converged_file) as f:
            converged = f.read().strip() == "TRUE"
        if not converged:
            logger.warning("MESuSiE did not converge. Results may be unreliable.")

    if converged is False and empty_on_nonconvergence:
        logger.error(
            "MESuSiE did not converge in %d iterations; returning empty "
            "credible set because empty_on_nonconvergence=True.",
            max_iter,
        )
        zero_pips = pd.Series(
            data=np.zeros(len(common_snps_list), dtype=float),
            index=common_snps_list,
        )
        return CredibleSet(
            tool=Method.MESUSIE,
            n_cs=0,
            coverage=coverage,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=zero_pips,
            parameters=parameters,
            converged=False,
        )

    # Build PIP series
    pip = pd.Series(
        index=pip_df["SNP"].values.tolist(), data=pip_df["PIP"].values.tolist()
    )

    # Build credible sets
    cs_snp: List[List[str]] = []
    purity_list: List[Optional[float]] = []
    cs_types: List[str] = []

    if len(cs_df) > 0:
        for cs_id, sub_df in cs_df.groupby("CS_ID"):
            cs_snp.append(sub_df["SNP"].values.tolist())

            # Get purity for this CS
            cs_purity_row = purity_df[purity_df["CS_ID"] == cs_id]
            if len(cs_purity_row) > 0 and not pd.isna(cs_purity_row.iloc[0]["PURITY"]):
                purity_list.append(float(cs_purity_row.iloc[0]["PURITY"]))
            else:
                purity_list.append(None)

            # Get CS type
            if len(cs_purity_row) > 0 and "CS_TYPE" in cs_purity_row.columns:
                cs_types.append(str(cs_purity_row.iloc[0]["CS_TYPE"]))
            else:
                cs_types.append("unknown")
    else:
        logger.warning("No credible set found, please try other parameters.")

    # Store cs_types in parameters
    if cs_types:
        parameters["cs_types"] = cs_types

    logger.info(f"Finished MESuSiE on {locus_set}")
    logger.info(f"N of credible set: {len(cs_snp)}")
    logger.info(f"Credible set size: {[len(i) for i in cs_snp]}")
    logger.info(f"Credible set purity: {purity_list}")
    if cs_types:
        logger.info(f"Credible set types: {cs_types}")

    return CredibleSet(
        tool=Method.MESUSIE,
        n_cs=len(cs_snp),
        coverage=coverage,
        lead_snps=[str(pip[pip.index.isin(i)].idxmax()) for i in cs_snp],
        snps=cs_snp,
        cs_sizes=[len(i) for i in cs_snp],
        pips=pip,
        parameters=parameters,
        purity=purity_list if len(purity_list) > 0 else None,
        converged=converged,
    )
