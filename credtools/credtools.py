"""Main module."""

import inspect
import json
import logging
import os
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import numpy as np
import pandas as pd
import toml  # type: ignore

from credtools.cojo import conditional_selection
from credtools.credibleset import CredibleSet, combine_creds
from credtools.locus import LocusSet, load_locus_set
from credtools.meta import meta
from credtools.qc import locus_qc
from credtools.wrappers import (
    run_abf,
    run_carma,
    run_finemap,
    run_multisusie,
    run_rsparsepro,
    run_susie,
    run_susiex,
)

logger = logging.getLogger("CREDTOOLS")


def fine_map(
    locus_set: LocusSet,
    strategy: str = "single_input",
    tool: str = "susie",
    max_causal: int = 1,
    set_L_by_cojo: bool = True,
    p_cutoff: float = 5e-8,
    collinear_cutoff: float = 0.9,
    window_size: int = 10000000,
    maf_cutoff: float = 0.01,
    diff_freq_cutoff: float = 0.2,
    combine_cred: str = "union",
    combine_pip: str = "max",
    jaccard_threshold: float = 0.1,
    **kwargs,
) -> CredibleSet:
    """
    Perform fine-mapping on a locus set.

    Parameters
    ----------
    locus_set : LocusSet
        Locus set to fine-mapping.
    strategy : str
        Fine-mapping strategy. Choose from ["single_input", "multi_input", "post_hoc_combine"]
        single_input: use fine-mapping tools which take a single locus as input, these tools are:
            abf, carma, finemap, rsparsepro, susie
        multi_input: use fine-mapping tools which take multiple loci as input, these tools are:
            multisusie, susiex
        post_hoc_combine: use fine-mapping tools which take single loci as input (see single_input options), and then combine the results,
    tool : str
        Fine-mapping tool. Choose from ["abf", "carma", "finemap", "rsparsepro", "susie", "multisusie", "susiex"]
    combine_cred : str, optional
        Method to combine credible sets, by default "union".
        Options: "union", "intersection", "cluster".
        "union":        Union of all credible sets to form a merged credible set.
        "intersection": Frist merge the credible sets from the same tool,
                        then take the intersection of all merged credible sets.
                        no credible set will be returned if no common SNPs found.
        "cluster":      Merge credible sets with Jaccard index > 0.1.
    combine_pip : str, optional
        Method to combine PIPs, by default "max".
        Options: "max", "min", "mean", "meta".
        "meta": PIP_meta = 1 - prod(1 - PIP_i), where i is the index of tools,
                PIP_i = 0 when the SNP is not in the credible set of the tool.
        "max":  Maximum PIP value for each SNP across all tools.
        "min":  Minimum PIP value for each SNP across all tools.
        "mean": Mean PIP value for each SNP across all tools.
    jaccard_threshold : float, optional
        Jaccard index threshold for the "cluster" method, by default 0.1.
    max_causal : int, optional
        Maximum number of causal variants, by default 1.
    """
    tool_func_dict: Dict[str, Callable[..., Any]] = {
        "abf": run_abf,
        "carma": run_carma,
        "finemap": run_finemap,
        "rsparsepro": run_rsparsepro,
        "susie": run_susie,
        "multisusie": run_multisusie,
        "susiex": run_susiex,
    }
    inspect_dict = {
        "abf": set(inspect.signature(run_abf).parameters),
        "carma": set(inspect.signature(run_carma).parameters),
        "finemap": set(inspect.signature(run_finemap).parameters),
        "rsparsepro": set(inspect.signature(run_rsparsepro).parameters),
        "susie": set(inspect.signature(run_susie).parameters),
        "multisusie": set(inspect.signature(run_multisusie).parameters),
        "susiex": set(inspect.signature(run_susiex).parameters),
    }
    params_dict = {}
    for t, args in inspect_dict.items():
        params_dict[t] = {k: v for k, v in kwargs.items() if k in args}
    if strategy == "single_input":
        if locus_set.n_loci > 1:
            raise ValueError(
                "Locus set must contain only one locus for single-input strategy"
            )
        if tool in ["abf", "carma", "finemap", "rsparsepro", "susie"]:
            if set_L_by_cojo:
                max_causal = len(
                    conditional_selection(
                        locus_set.loci[0],
                        p_cutoff=p_cutoff,
                        collinear_cutoff=collinear_cutoff,
                        window_size=window_size,
                        maf_cutoff=maf_cutoff,
                        diff_freq_cutoff=diff_freq_cutoff,
                    )
                )
                if max_causal == 0:
                    logger.warning(
                        "No significant SNPs found by COJO, using max_causal=1"
                    )
                    max_causal = 1
            return tool_func_dict[tool](
                locus_set.loci[0], max_causal=max_causal, **params_dict[tool]
            )
        else:
            raise ValueError(f"Tool {tool} not supported for single-input strategy")
    elif strategy == "multi_input":
        # if locus_set.n_loci < 2:
        #     raise ValueError("Locus set must contain at least two loci for multi-input strategy")
        if tool in ["multisusie", "susiex"]:
            return tool_func_dict[tool](
                locus_set, max_causal=max_causal, **params_dict[tool]
            )
        else:
            raise ValueError(f"Tool {tool} not supported for multi-input strategy")
    elif strategy == "post_hoc_combine":
        # if locus_set.n_loci < 2:
        #     raise ValueError("Locus set must contain at least two loci for post-hoc combine strategy")
        if tool in ["abf", "carma", "finemap", "rsparsepro", "susie"]:
            all_creds = []
            for locus in locus_set.loci:
                creds = tool_func_dict[tool](
                    locus, max_causal=max_causal, **params_dict[tool]
                )
                all_creds.append(creds)
            return combine_creds(
                all_creds,
                combine_cred=combine_cred,
                combine_pip=combine_pip,
                jaccard_threshold=jaccard_threshold,
            )
        else:
            raise ValueError(f"Tool {tool} not supported for post-hoc combine strategy")
    else:
        raise ValueError(f"Strategy {strategy} not supported")


def pipeline(
    loci_df: pd.DataFrame,
    meta_method: str = "meta_all",
    skip_qc: bool = False,
    strategy: str = "single_input",
    tool: str = "susie",
    outdir: str = ".",
    **kwargs,
):
    """
    Run whole fine-mapping pipeline on a list of loci.

    Parameters
    ----------
    loci_df : pd.DataFrame
        Dataframe containing the locus information.
    meta_method : str, optional
        Meta-analysis method, by default "meta_all"
        Options: "meta_all", "meta_by_population", "no_meta".
    skip_qc : bool, optional
        Skip QC, by default False.
    strategy : str, optional
        Fine-mapping strategy, by default "single_input".
    tool : str, optional
        Fine-mapping tool, by default "susie".
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    locus_set = load_locus_set(loci_df)
    # meta-analysis
    locus_set = meta(locus_set, meta_method=meta_method)
    logger.info(f"Meta-analysis complete, {locus_set.n_loci} loci loaded.")
    logger.info(f"Save meta-analysis results to {outdir}.")
    for locus in locus_set.loci:
        out_prefix = f"{outdir}/{locus.prefix}"
        locus.sumstats.to_csv(f"{out_prefix}.sumstat", sep="\t", index=False)
        np.savez_compressed(f"{out_prefix}.ld.npz", ld=locus.ld.r.astype(np.float16))
        locus.ld.map.to_csv(f"{out_prefix}.ldmap", sep="\t", index=False)
    # QC
    if not skip_qc:
        qc_metrics = locus_qc(locus_set)
        logger.info(f"QC complete, {qc_metrics.keys()} metrics saved.")
        for k, v in qc_metrics.items():
            v.to_csv(f"{outdir}/{k}.txt", sep="\t", index=False)
    # fine-mapping
    creds = fine_map(locus_set, strategy=strategy, tool=tool, **kwargs)
    creds.pips.to_csv(f"{outdir}/pips.txt", sep="\t", header=False, index=True)
    with open(f"{outdir}/creds.json", "w") as f:
        json.dump(creds.to_dict(), f, indent=4)
    logger.info(f"Fine-mapping complete, {creds.n_cs} credible sets saved.")
    return
