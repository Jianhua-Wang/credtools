"""Wrapper for CARMA fine-mapping (CAusal Robust Mapping Method with Annotations)."""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet, calculate_cs_purity
from credtools.locus import Locus, intersect_sumstat_ld

logger = logging.getLogger("CARMA")

_R_SCRIPT_PATH = str(Path(__file__).parent / "carma_wrapper.R")


def _check_r_and_carma() -> None:
    """Check that Rscript and the CARMA R package are available."""
    if shutil.which("Rscript") is None:
        raise FileNotFoundError(
            "Rscript not found on PATH. Please install R "
            "(https://cran.r-project.org/) and ensure Rscript is available."
        )

    result = subprocess.run(
        ["Rscript", "-e", "library(CARMA)"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise FileNotFoundError(
            "CARMA R package is not installed. Please install it with:\n"
            "  R -e 'devtools::install_github(\"ZikunY/CARMA\")'"
        )


def run_carma(
    locus: Locus,
    max_causal: int = 10,
    coverage: float = 0.99,
    effect_size_prior: str = "Spike-slab",
    bf_threshold: float = 10.0,
    outlier_switch: bool = True,
    outlier_bf_threshold: float = 1.0 / 3.2,
    em_dist: str = "Logistic",
    max_model_dim: int = 200000,
    all_iter: int = 3,
    all_inner_iter: int = 10,
    input_alpha: float = 0.0,
    epsilon_threshold: float = 1e-5,
    tau: float = 0.04,
    y_var: float = 1.0,
    significant_threshold: float = 5e-8,
) -> CredibleSet:
    """
    Run CARMA fine-mapping with summary statistics and an LD matrix.

    CARMA (CAusal Robust Mapping Method with Annotations) is a Bayesian
    fine-mapping framework that explicitly models discrepancies between the
    summary statistics and the reference LD as outliers, providing robustness
    when the in-sample LD is unavailable.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix data.
    max_causal : int, optional
        Maximum number of causal variants (CARMA's ``num.causal``), by
        default 10.
    coverage : float, optional
        Cumulative posterior probability for credible sets (CARMA's
        ``rho.index``), by default 0.99.
    effect_size_prior : {"Spike-slab", "Cauchy"}, optional
        Prior on the causal effect size, by default "Spike-slab".
    bf_threshold : float, optional
        Bayes-factor threshold to keep models in the credible model set
        (CARMA's ``BF.index``), by default 10.0.
    outlier_switch : bool, optional
        Enable CARMA's outlier-detection step that identifies SNPs whose
        marginal effect is inconsistent with the LD-implied joint model.
        Defaults to True.
    outlier_bf_threshold : float, optional
        Bayes-factor threshold for outlier inclusion, by default 1/3.2.
    em_dist : str, optional
        EM distribution for the annotation prior, by default "Logistic".
    max_model_dim : int, optional
        Maximum number of candidate causal models, by default 2e5.
    all_iter : int, optional
        Number of outer EM iterations, by default 3.
    all_inner_iter : int, optional
        Number of inner shotgun-stochastic-search iterations, by default 10.
    input_alpha : float, optional
        Annotation regression intercept prior, by default 0.0.
    epsilon_threshold : float, optional
        Convergence tolerance, by default 1e-5.
    tau : float, optional
        Effect-size prior variance, by default 0.04.
    y_var : float, optional
        Phenotype variance assumption, by default 1.0.
    significant_threshold : float, optional
        Minimum p-value for the locus to be fine-mapped. If no variants
        cross this threshold, an empty credible set is returned without
        invoking R. Defaults to 5e-8.

    Returns
    -------
    CredibleSet
        Credible set object with PIPs, credible sets, lead SNPs, purity,
        and CARMA-specific metadata (e.g., detected outlier SNPs).

    Raises
    ------
    FileNotFoundError
        If Rscript or the CARMA R package is not available.
    RuntimeError
        If the CARMA R script execution fails.

    Notes
    -----
    The CARMA R wrapper is bundled with credtools and is invoked via
    ``Rscript`` in a temporary directory. Inputs are written as CSV/binary
    files; outputs are read back and translated into a :class:`CredibleSet`.

    Reference
    ---------
    Yang, Z. et al. CARMA is a new Bayesian model for fine-mapping in
    genome-wide association meta-analyses. *Nat. Genet.* **55**,
    1057-1065 (2023).
    """
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "effect_size_prior": effect_size_prior,
        "bf_threshold": bf_threshold,
        "outlier_switch": outlier_switch,
        "outlier_bf_threshold": outlier_bf_threshold,
        "em_dist": em_dist,
        "max_model_dim": max_model_dim,
        "all_iter": all_iter,
        "all_inner_iter": all_inner_iter,
        "input_alpha": input_alpha,
        "epsilon_threshold": epsilon_threshold,
        "tau": tau,
        "y_var": y_var,
        "significant_threshold": significant_threshold,
    }
    logger.info(f"Running CARMA on {locus}")
    logger.info(f"Parameters: {parameters}")

    if not (locus.sumstats[ColName.P] <= significant_threshold).any():
        logger.warning(
            "No variants pass the significance threshold %.2e. "
            "Returning empty result.",
            significant_threshold,
        )
        zero_pips = pd.Series(
            data=np.zeros(len(locus.sumstats), dtype=float),
            index=locus.sumstats[ColName.SNPID].tolist(),
        )
        return CredibleSet(
            tool=Method.CARMA,
            n_cs=0,
            coverage=coverage,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=zero_pips,
            parameters=parameters,
        )

    _check_r_and_carma()

    if not locus.is_matched:
        logger.warning(
            "The sumstat and LD are not matched, will match them in same order."
        )
        locus = intersect_sumstat_ld(locus)

    import tempfile

    with tempfile.TemporaryDirectory(prefix="carma_") as temp_dir:
        sumstat = locus.sumstats.copy()
        sumstat[ColName.Z] = sumstat[ColName.BETA] / sumstat[ColName.SE]
        ss_out = pd.DataFrame(
            {
                "SNP": sumstat[ColName.SNPID].values,
                "Z": sumstat[ColName.Z].values,
            }
        )
        ss_out.to_csv(f"{temp_dir}/sumstats.csv", index=False)

        snpids = sumstat[ColName.SNPID].tolist()
        with open(f"{temp_dir}/snpids.txt", "w") as f:
            f.write("\n".join(snpids) + "\n")

        ld_r = locus.ld.r.astype(np.float64)
        ld_r.tofile(f"{temp_dir}/ld.bin")
        with open(f"{temp_dir}/ld_dim.txt", "w") as f:
            f.write(str(ld_r.shape[0]) + "\n")

        cmd = [
            "Rscript",
            _R_SCRIPT_PATH,
            "--temp_dir",
            temp_dir,
            "--num_causal",
            str(max_causal),
            "--rho_index",
            str(coverage),
            "--bf_index",
            str(bf_threshold),
            "--outlier_switch",
            str(outlier_switch).upper(),
            "--outlier_bf_index",
            str(outlier_bf_threshold),
            "--effect_size_prior",
            effect_size_prior,
            "--em_dist",
            em_dist,
            "--max_model_dim",
            str(max_model_dim),
            "--all_iter",
            str(all_iter),
            "--all_inner_iter",
            str(all_inner_iter),
            "--input_alpha",
            str(input_alpha),
            "--epsilon_threshold",
            str(epsilon_threshold),
            "--tau",
            str(tau),
            "--y_var",
            str(y_var),
        ]
        logger.info(f"Running CARMA R script: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"CARMA R script stderr: {result.stderr}")
            logger.error(f"CARMA R script stdout: {result.stdout}")
            raise RuntimeError(
                f"CARMA R script failed (return code {result.returncode}).\n"
                f"stderr: {result.stderr}\n"
                f"stdout: {result.stdout}"
            )
        logger.debug(f"CARMA R script stdout: {result.stdout}")

        pip_df = pd.read_csv(f"{temp_dir}/carma_pips.csv")
        cs_df = pd.read_csv(f"{temp_dir}/carma_cs.csv")
        outliers_path = f"{temp_dir}/carma_outliers.csv"
        outliers: List[str] = []
        if os.path.exists(outliers_path):
            outliers_df = pd.read_csv(outliers_path)
            if len(outliers_df) > 0 and "SNP" in outliers_df.columns:
                outliers = outliers_df["SNP"].astype(str).tolist()

    pips = pd.Series(
        index=pip_df["SNP"].astype(str).tolist(),
        data=pip_df["PIP"].astype(float).tolist(),
    )

    cs_snps: List[List[str]] = []
    if len(cs_df) > 0:
        for _, sub_df in cs_df.groupby("CS_ID"):
            cs_snps.append(sub_df["SNP"].astype(str).tolist())
    else:
        logger.warning("CARMA found no credible sets.")

    lead_snps = [str(pips[pips.index.isin(s)].idxmax()) for s in cs_snps]
    cs_sizes = [len(s) for s in cs_snps]

    purity: Optional[List[float]] = None
    if cs_snps and locus.ld is not None:
        purity = [calculate_cs_purity(locus.ld, s) for s in cs_snps]

    parameters["outliers"] = outliers

    logger.info(f"Finished CARMA on {locus}")
    logger.info(f"N of credible set: {len(cs_snps)}")
    logger.info(f"Credible set size: {cs_sizes}")
    if outliers:
        logger.info(f"Detected {len(outliers)} outlier SNPs")

    return CredibleSet(
        tool=Method.CARMA,
        n_cs=len(cs_snps),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=cs_sizes,
        pips=pips,
        parameters=parameters,
        purity=purity,
    )
