"""Export finemapping results for web visualization."""

import json
import os
from functools import partial
from multiprocessing import Pool
from subprocess import call
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

from credtools.locus import load_locus


def process_locus(
    locus_id: str,
    loci_info: pd.DataFrame,
    data_dir: str,
    tool_list: List[str],
    webdata_dir: str,
) -> None:
    """Process a single locus for web visualization.

    Args
    ----
    locus_id : str
        The locus identifier
    loci_info : pd.DataFrame
        DataFrame containing locus information
    data_dir : str
        Directory containing the data
    tool_list : List[str]
        List of fine-mapping tools to process
    webdata_dir : str
        Output directory for web data
    """
    locus_dir = f"{webdata_dir}/{data_dir}/{locus_id}"
    os.makedirs(locus_dir, exist_ok=True)
    locus_info_subset = loci_info[loci_info["locus_id"] == locus_id]

    for i, row in locus_info_subset.iterrows():
        prefix, popu, cohort, sample_size = (
            row["prefix"],
            row["popu"],
            row["cohort"],
            row["sample_size"],
        )
        orig_sumstats = pd.read_csv(f"{prefix}.sumstats.gz", sep="\t")
        locus = load_locus(prefix, popu, cohort, sample_size, if_intersect=True)
        lead_idx = locus.sumstats["P"].idxmin()
        locus.sumstats["r2"] = locus.ld.r[lead_idx, :]  # type: ignore
        locus.sumstats["r2"] = locus.sumstats["r2"] ** 2
        orig_sumstats["r2"] = orig_sumstats["SNPID"].map(  # type: ignore
            pd.Series(locus.sumstats["r2"].values, index=locus.sumstats["SNPID"].values)  # type: ignore
        )

        for tool in tool_list:
            cred_dir = f"data/real/credset/{data_dir}/{tool}/{locus_id}"
            orig_sumstats[f"{tool}_cred"] = 0
            if os.path.exists(f"{cred_dir}/creds.json"):
                creds = json.load(open(f"{cred_dir}/creds.json"))
                cred_snps = creds["snps"]
                if len(cred_snps) > 0:
                    for cred_idx, snps in enumerate(cred_snps):
                        orig_sumstats.loc[
                            orig_sumstats["SNPID"].isin(snps), f"{tool}_cred"
                        ] = (cred_idx + 1)
            if os.path.exists(f"{cred_dir}/pips.txt"):
                pip_df = pd.read_csv(
                    f"{cred_dir}/pips.txt",
                    sep="\t",
                    names=["SNPID", "pip"],
                    index_col=0,
                )
                orig_sumstats[f"{tool}_pip"] = orig_sumstats["SNPID"].map(pip_df["pip"])
            else:
                orig_sumstats[f"{tool}_pip"] = 0
        orig_sumstats.to_csv(
            f"{locus_dir}/{popu}_{cohort}.res.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )


def process_qc(
    locus_id: str,
    loci_info: pd.DataFrame,
    data_dir: str,
    tool_list: List[str],
    webdata_dir: str,
) -> None:
    """Process QC data for a single locus.

    Args
    ----
    locus_id : str
        The locus identifier
    loci_info : pd.DataFrame
        DataFrame containing locus information
    data_dir : str
        Directory containing the data
    tool_list : List[str]
        List of fine-mapping tools to process
    webdata_dir : str
        Output directory for web data
    """
    locus_dir = f"{webdata_dir}/{data_dir}/{locus_id}"
    os.makedirs(locus_dir, exist_ok=True)
    locus_info_subset = loci_info[loci_info["locus_id"] == locus_id]
    qc_df = pd.DataFrame(
        columns=[
            "meta_type",
            "tool",
            "locus_id",
            "chr",
            "start",
            "end",
            "n_credsets",
            "whole_credsize",
            "n_PIP_gt_0.1",
            "popu",
            "cohort",
            "nsnp",
            "nsnp_1e-5",
            "nsnp_5e-8",
            "lambda",
            "dentist-s",
            "maf_corr",
            "ld_decay_rate",
        ]
    )
    ith = 0
    for i, row in locus_info_subset.iterrows():
        prefix, popu, cohort, sample_size = (
            row["prefix"],
            row["popu"],
            row["cohort"],
            row["sample_size"],
        )
        cred_res = pd.read_csv(f"{locus_dir}/{popu}_{cohort}.res.gz", sep="\t")
        for tool in tool_list:
            n_credsets = len(cred_res[f"{tool}_cred"].unique()) - 1
            whole_credsize = cred_res[cred_res[f"{tool}_cred"] > 0].shape[0]
            n_PIP_gt_01 = cred_res[cred_res[f"{tool}_pip"] > 0.1].shape[0]
            qc_res_dir = f"data/real/qc/{data_dir}/{locus_id}"
            lambda_res = pd.read_csv(f"{qc_res_dir}/expected_z.txt.gz", sep="\t")
            lambda_s = lambda_res[lambda_res["cohort"] == f"{popu}_{cohort}"][
                "lambda_s"
            ].max()
            dentist_s = pd.read_csv(f"{qc_res_dir}/dentist_s.txt.gz", sep="\t")
            dentist_s["p_dentist_s"] = stats.chi2.logsf(
                dentist_s["t_dentist_s"], df=1
            ) / -np.log(10)
            sig_dentist_s = dentist_s[dentist_s["p_dentist_s"] > 4]["SNPID"].unique()
            sig_dentist_s = len(
                cred_res[
                    (cred_res["r2"] > 0.6) & (cred_res["SNPID"].isin(sig_dentist_s))
                ]
            )
            maf_df = pd.read_csv(f"{qc_res_dir}/compare_maf.txt.gz", sep="\t")
            maf_df = maf_df[maf_df["cohort"] == f"{popu}_{cohort}"]
            if maf_df.shape[0] > 0:
                maf_corr = stats.pearsonr(maf_df["MAF_sumstats"], maf_df["MAF_ld"])[0]
            else:
                maf_corr = 1
            ld_decay_rate = pd.read_csv(f"{qc_res_dir}/ld_decay.txt.gz", sep="\t")
            ld_decay_rate = ld_decay_rate[
                ld_decay_rate["cohort"] == f"{popu}_{cohort}"
            ]["decay_rate"].max()
            qc_df.loc[ith] = [
                data_dir,
                tool,
                locus_id,
                row["chr"],
                row["start"],
                row["end"],
                n_credsets,
                whole_credsize,
                n_PIP_gt_01,
                popu,
                cohort,
                cred_res.shape[0],
                cred_res[cred_res["P"] < 1e-5].shape[0],
                cred_res[cred_res["P"] < 5e-8].shape[0],
                lambda_s,
                sig_dentist_s,
                maf_corr,
                ld_decay_rate,
            ]
            ith += 1
    call(f"cp {qc_res_dir}/expected_z.txt.gz {locus_dir}/expected_z.txt.gz", shell=True)
    call(
        f"cp {qc_res_dir}/ld_4th_moment.txt.gz {locus_dir}/ld_4th_moment.txt.gz",
        shell=True,
    )
    if os.path.exists(f"{qc_res_dir}/snp_missingness.txt.gz"):
        call(
            f"cp {qc_res_dir}/snp_missingness.txt.gz {locus_dir}/snp_missingness.txt.gz",
            shell=True,
        )
    call(f"cp {qc_res_dir}/ld_decay.txt.gz {locus_dir}/ld_decay.txt.gz", shell=True)
    qc_df.to_csv(f"{locus_dir}/qc.txt.gz", sep="\t", index=False, compression="gzip")


def export_for_web(
    data_base_dir: str,
    webdata_dir: str = "webdata",
    allmeta_loci_file: Optional[str] = None,
    popumeta_loci_file: Optional[str] = None,
    nometa_loci_file: Optional[str] = None,
    threads: int = 10,
) -> None:
    """Export finemapping results for web visualization.

    Args
    ----
    data_base_dir : str
        Base directory containing the finemapping data
    webdata_dir : str
        Output directory for processed web data
    allmeta_loci_file : Optional[str]
        Path to allmeta loci info file
    popumeta_loci_file : Optional[str]
        Path to popumeta loci info file
    nometa_loci_file : Optional[str]
        Path to nometa loci info file
    threads : int
        Number of parallel threads to use
    """
    # Change to the base directory
    os.chdir(data_base_dir)
    os.makedirs(webdata_dir, exist_ok=True)

    all_loci_info = []

    # Process allmeta if provided
    if allmeta_loci_file and os.path.exists(allmeta_loci_file):
        allmeta_loci = pd.read_csv(allmeta_loci_file, sep="\t")
        tool_list = ["abf", "susie", "finemap", "rsparsepro", "carma"]

        # Process loci
        partial_func = partial(
            process_locus,
            loci_info=allmeta_loci,
            data_dir="allmeta",
            tool_list=tool_list,
            webdata_dir=webdata_dir,
        )
        with Pool(threads) as p:
            p.map(partial_func, allmeta_loci["locus_id"].unique())

        # Process QC
        partial_func = partial(
            process_qc,
            loci_info=allmeta_loci,
            data_dir="allmeta",
            tool_list=tool_list,
            webdata_dir=webdata_dir,
        )
        with Pool(threads) as p:
            p.map(partial_func, allmeta_loci["locus_id"].unique())

        # Collect loci info
        for locus_id in allmeta_loci["locus_id"].unique():
            locus_info = pd.read_csv(
                f"{webdata_dir}/allmeta/{locus_id}/qc.txt.gz",
                sep="\t",
                compression="gzip",
            )
            all_loci_info.append(locus_info)

    # Process popumeta if provided
    if popumeta_loci_file and os.path.exists(popumeta_loci_file):
        popumeta_loci = pd.read_csv(popumeta_loci_file, sep="\t")
        tool_list = [
            "abf",
            "susie",
            "finemap",
            "rsparsepro",
            "carma",
            "susiex",
            "multisusie",
        ]

        # Process loci
        partial_func = partial(
            process_locus,
            loci_info=popumeta_loci,
            data_dir="popumeta",
            tool_list=tool_list,
            webdata_dir=webdata_dir,
        )
        with Pool(threads) as p:
            p.map(partial_func, popumeta_loci["locus_id"].unique())

        # Process QC
        partial_func = partial(
            process_qc,
            loci_info=popumeta_loci,
            data_dir="popumeta",
            tool_list=tool_list,
            webdata_dir=webdata_dir,
        )
        with Pool(threads) as p:
            p.map(partial_func, popumeta_loci["locus_id"].unique())

        # Collect loci info
        for locus_id in popumeta_loci["locus_id"].unique():
            locus_info = pd.read_csv(
                f"{webdata_dir}/popumeta/{locus_id}/qc.txt.gz",
                sep="\t",
                compression="gzip",
            )
            all_loci_info.append(locus_info)

    # Process nometa if provided
    if nometa_loci_file and os.path.exists(nometa_loci_file):
        nometa_loci = pd.read_csv(nometa_loci_file, sep="\t")
        tool_list = ["susiex", "multisusie"]

        # Process loci
        partial_func = partial(
            process_locus,
            loci_info=nometa_loci,
            data_dir="nometa",
            tool_list=tool_list,
            webdata_dir=webdata_dir,
        )
        with Pool(threads) as p:
            p.map(partial_func, nometa_loci["locus_id"].unique())

        # Process QC
        partial_func = partial(
            process_qc,
            loci_info=nometa_loci,
            data_dir="nometa",
            tool_list=tool_list,
            webdata_dir=webdata_dir,
        )
        with Pool(threads) as p:
            p.map(partial_func, nometa_loci["locus_id"].unique())

        # Collect loci info
        for locus_id in nometa_loci["locus_id"].unique():
            locus_info = pd.read_csv(
                f"{webdata_dir}/nometa/{locus_id}/qc.txt.gz",
                sep="\t",
                compression="gzip",
            )
            all_loci_info.append(locus_info)

    # Combine all loci info
    if all_loci_info:
        all_loci_info_df = pd.concat(all_loci_info, ignore_index=True)
        all_loci_info_df.to_csv(
            f"{webdata_dir}/all_loci_info.txt", sep="\t", index=False
        )

    print(f"Web data exported to {webdata_dir}")
