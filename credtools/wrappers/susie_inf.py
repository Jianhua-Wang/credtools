"""Wrapper for SuSiE-inf fine-mapping (SuSiE 2.0 with infinitesimal background)."""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet, calculate_cs_purity
from credtools.locus import Locus, intersect_sumstat_ld

logger = logging.getLogger("SUSIE_INF")

_R_SCRIPT_PATH = str(Path(__file__).parent / "susie_inf_wrapper.R")


def _check_r_and_susie_inf() -> None:
    """Check that Rscript and a SuSiE-2.0-capable susieR are available."""
    if shutil.which("Rscript") is None:
        raise FileNotFoundError(
            "Rscript not found on PATH. Please install R "
            "(https://cran.r-project.org/) and ensure Rscript is available."
        )

    result = subprocess.run(
        ["Rscript", "-e", "library(susieR)"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise FileNotFoundError(
            "susieR R package is not installed. Please install it with:\n"
            "  R -e 'install.packages(\"susieR\")'\n"
            "SuSiE-inf requires susieR >= 0.16.1 (the SuSiE 2.0 release)."
        )


def _read_status(temp_dir: str) -> dict:
    """Parse the status file written by the R wrapper."""
    path = Path(temp_dir) / "susie_inf_status.txt"
    status: dict = {}
    if not path.exists():
        return status
    for raw in path.read_text().splitlines():
        if "\t" not in raw:
            continue
        key, value = raw.split("\t", 1)
        status[key.strip()] = value.strip()
    return status


def run_susie_inf(
    locus: Locus,
    max_causal: int = 10,
    coverage: float = 0.95,
    max_iter: int = 100,
    estimate_residual_variance: bool = False,
    purity: float = 0.0,
    convergence_tol: float = 1e-3,
    significant_threshold: float = 5e-8,
    empty_on_nonconvergence: bool = False,
) -> CredibleSet:
    """
    Run SuSiE-inf fine-mapping with summary statistics and an LD matrix.

    SuSiE-inf extends SuSiE with a single-Gaussian "infinitesimal" prior
    over the unmappable / background effects, capturing a broad polygenic
    background where many variants contribute small effects of comparable
    magnitude. It is well suited to loci that look like a strong peak
    rising out of a high plateau of background signal.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix data.
    max_causal : int, optional
        Maximum number of single-effect components (SuSiE's L), by default 10.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
    max_iter : int, optional
        Maximum number of IBSS iterations, by default 100.
    estimate_residual_variance : bool, optional
        Whether to estimate residual variance from data, by default False.
    purity : float, optional
        Minimum absolute correlation (susieR's ``min_abs_corr``) for credible
        sets, by default 0.0.
    convergence_tol : float, optional
        Convergence tolerance for the ELBO, by default 1e-3.
    significant_threshold : float, optional
        Minimum p-value required for the locus to be fine-mapped. If no
        variants cross this threshold, an empty credible set is returned
        without invoking R. Defaults to 5e-8.
    empty_on_nonconvergence : bool, optional
        When True and SuSiE-inf did not converge, return an empty credible
        set instead of the unreliable partially-fit results. Defaults to
        False (preserve the SuSiE wrapper's legacy behaviour).

    Returns
    -------
    CredibleSet
        Credible set object with PIPs, credible sets, lead SNPs, purity,
        and SuSiE-inf convergence metadata.

    Raises
    ------
    FileNotFoundError
        If Rscript or the susieR R package is not available.
    RuntimeError
        If the SuSiE-inf R script execution fails.

    Notes
    -----
    The R wrapper invokes ``susieR::susie_rss(..., unmappable_effects="inf")``
    in a temporary directory. Inputs are written as CSV/binary files; outputs
    are read back and translated into a :class:`CredibleSet`.

    Reference
    ---------
    McCreight, J. et al. SuSiE 2.0: a flexible Bayesian framework for
    fine-mapping with adaptive shrinkage and infinitesimal background
    effects. *bioRxiv* (2025), DOI: 10.1101/2025.11.25.690514.
    """
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "max_iter": max_iter,
        "estimate_residual_variance": estimate_residual_variance,
        "min_abs_corr": purity,
        "convergence_tol": convergence_tol,
        "estimate_prior_method": "optim",
        "unmappable_effects": "inf",
        "significant_threshold": significant_threshold,
        "empty_on_nonconvergence": empty_on_nonconvergence,
    }
    logger.info(f"Running SuSiE-inf on {locus}")
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
            tool=Method.SUSIE_INF,
            n_cs=0,
            coverage=coverage,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=zero_pips,
            parameters=parameters,
        )

    _check_r_and_susie_inf()

    if not locus.is_matched:
        logger.warning(
            "The sumstat and LD are not matched, will match them in same order."
        )
        locus = intersect_sumstat_ld(locus)

    import tempfile

    with tempfile.TemporaryDirectory(prefix="susie_inf_") as temp_dir:
        sumstat = locus.sumstats
        ss_out = pd.DataFrame(
            {
                "SNP": sumstat[ColName.SNPID].values,
                "BHAT": sumstat[ColName.BETA].values.astype(float),
                "SHAT": sumstat[ColName.SE].values.astype(float),
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
            "--n",
            str(int(locus.sample_size)),
            "--L",
            str(max_causal),
            "--coverage",
            str(coverage),
            "--max_iter",
            str(max_iter),
            "--estimate_residual_variance",
            str(estimate_residual_variance).upper(),
            "--min_abs_corr",
            str(purity),
            "--tol",
            str(convergence_tol),
        ]
        logger.info(f"Running SuSiE-inf R script: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"SuSiE-inf R script stderr: {result.stderr}")
            logger.error(f"SuSiE-inf R script stdout: {result.stdout}")
            raise RuntimeError(
                f"SuSiE-inf R script failed (return code {result.returncode}).\n"
                f"stderr: {result.stderr}\n"
                f"stdout: {result.stdout}"
            )
        logger.debug(f"SuSiE-inf R script stdout: {result.stdout}")

        pip_df = pd.read_csv(f"{temp_dir}/susie_inf_pips.csv")
        cs_df = pd.read_csv(f"{temp_dir}/susie_inf_cs.csv")
        status = _read_status(temp_dir)

    pips = pd.Series(
        index=pip_df["SNP"].astype(str).tolist(),
        data=pip_df["PIP"].astype(float).tolist(),
    )

    converged_str = status.get("converged", "").upper()
    converged: Optional[bool]
    if converged_str == "TRUE":
        converged = True
    elif converged_str == "FALSE":
        converged = False
    else:
        converged = None

    n_iter_str = status.get("niter", "")
    n_iter: Optional[int]
    try:
        n_iter = int(n_iter_str) if n_iter_str else None
    except ValueError:
        n_iter = None

    if converged is False and empty_on_nonconvergence:
        logger.error(
            "SuSiE-inf did not converge in %d iterations; returning empty "
            "credible set because empty_on_nonconvergence=True.",
            max_iter,
        )
        zero_pips = pd.Series(
            data=np.zeros(len(pips), dtype=float),
            index=pips.index.tolist(),
        )
        return CredibleSet(
            tool=Method.SUSIE_INF,
            n_cs=0,
            coverage=coverage,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=zero_pips,
            parameters=parameters,
            converged=False,
            n_iter=n_iter,
        )

    cs_snps: List[List[str]] = []
    if len(cs_df) > 0:
        for _, sub_df in cs_df.groupby("CS_ID"):
            cs_snps.append(sub_df["SNP"].astype(str).tolist())
    else:
        logger.warning("SuSiE-inf found no credible sets.")

    lead_snps = [str(pips[pips.index.isin(s)].idxmax()) for s in cs_snps]
    cs_sizes = [len(s) for s in cs_snps]

    purity_list: Optional[List[float]] = None
    if cs_snps and locus.ld is not None:
        purity_list = [calculate_cs_purity(locus.ld, s) for s in cs_snps]

    logger.info(f"Finished SuSiE-inf on {locus}")
    logger.info(f"N of credible set: {len(cs_snps)}")
    logger.info(f"Credible set size: {cs_sizes}")
    logger.info(f"converged={converged}, n_iter={n_iter}")

    return CredibleSet(
        tool=Method.SUSIE_INF,
        n_cs=len(cs_snps),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=cs_sizes,
        pips=pips,
        parameters=parameters,
        purity=purity_list,
        converged=converged,
        n_iter=n_iter,
    )
