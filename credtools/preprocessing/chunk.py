"""Chunk whole genome summary statistics into independent loci.

This module provides functionality to identify independent lead SNPs and create
regional chunks suitable for fine-mapping analysis across multiple ancestries.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from tqdm import tqdm

logger = logging.getLogger("Preprocessing")


def identify_independent_loci(
    sumstats_files: Union[Dict[str, str], str],
    output_dir: str,
    distance_threshold: int = 500000,
    pvalue_threshold: float = 5e-8,
    merge_overlapping: bool = True,
    use_most_sig_if_no_sig: bool = True,
    min_variants_per_locus: int = 10,
    **kwargs,
) -> pd.DataFrame:
    """
    Identify independent loci across multiple ancestries.

    Parameters
    ----------
    sumstats_files : Union[Dict[str, str], str]
        Dictionary mapping ancestry/cohort names to munged sumstats files,
        or single file path.
    output_dir : str
        Output directory for results.
    distance_threshold : int, optional
        Distance threshold in base pairs for independence, by default 500000.
    pvalue_threshold : float, optional
        P-value threshold for significance, by default 5e-8.
    merge_overlapping : bool, optional
        Whether to merge overlapping loci across ancestries, by default True.
    use_most_sig_if_no_sig : bool, optional
        Whether to use most significant SNP if no significant SNPs found, by default True.
    min_variants_per_locus : int, optional
        Minimum number of variants required per locus, by default 10.
    **kwargs
        Additional parameters.

    Returns
    -------
    pd.DataFrame
        DataFrame with identified loci coordinates and lead SNPs.

    Examples
    --------
    >>> files = {"EUR": "eur.munged.txt.gz", "ASN": "asn.munged.txt.gz"}
    >>> loci_df = identify_independent_loci(files, "output/")
    """
    os.makedirs(output_dir, exist_ok=True)

    # Normalize input to dictionary
    if isinstance(sumstats_files, str):
        ancestry_key = Path(sumstats_files).stem.replace(".munged", "")
        sumstats_files = {ancestry_key: sumstats_files}

    all_loci = []

    # Process each ancestry file
    for ancestry, file_path in tqdm(
        sumstats_files.items(), desc="Processing ancestries"
    ):
        logger.info(f"Processing {ancestry}: {file_path}")

        # Load sumstats
        sumstats = pd.read_csv(file_path, sep="\t", compression="gzip")

        # Identify independent SNPs for this ancestry
        ancestry_loci = _identify_independent_snps_by_distance(
            sumstats=sumstats,
            ancestry=ancestry,
            distance_threshold=distance_threshold,
            pvalue_threshold=pvalue_threshold,
            use_most_sig_if_no_sig=use_most_sig_if_no_sig,
            min_variants_per_locus=min_variants_per_locus,
        )

        all_loci.extend(ancestry_loci)

    # Convert to DataFrame
    loci_df = pd.DataFrame(all_loci)

    if len(loci_df) == 0:
        logger.warning("No loci identified")
        return loci_df

    # Merge overlapping loci across ancestries if requested
    if merge_overlapping and len(sumstats_files) > 1:
        logger.info("Merging overlapping loci across ancestries")
        loci_df = _merge_overlapping_loci(loci_df)

    # Sort by chromosome and position
    loci_df = loci_df.sort_values(["chr", "start"]).reset_index(drop=True)

    # Add locus IDs
    loci_df["locus_id"] = [
        f"chr{row['chr']}_{row['start']}_{row['end']}" for _, row in loci_df.iterrows()
    ]

    # Save results
    output_file = os.path.join(output_dir, "identified_loci.txt")
    loci_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Identified {len(loci_df)} loci, saved to {output_file}")

    return loci_df


def _identify_independent_snps_by_distance(
    sumstats: pd.DataFrame,
    ancestry: str,
    distance_threshold: int,
    pvalue_threshold: float,
    use_most_sig_if_no_sig: bool,
    min_variants_per_locus: int,
) -> List[Dict]:
    """
    Identify independent SNPs by distance for a single ancestry.

    Adapted from easyfinemap loci identification approach.
    """
    # Get significant SNPs
    sig_snps = sumstats[sumstats["P"] <= pvalue_threshold].copy()

    if len(sig_snps) == 0 and use_most_sig_if_no_sig:
        # Use most significant SNP from each chromosome
        logger.info(
            f"{ancestry}: No significant SNPs, using most significant per chromosome"
        )
        sig_snps = sumstats.loc[sumstats.groupby("CHR")["P"].idxmin()].copy()
    elif len(sig_snps) == 0:
        logger.warning(f"{ancestry}: No significant SNPs found")
        return []

    # Sort by p-value
    sig_snps = sig_snps.sort_values("P").reset_index(drop=True)

    independent_loci = []
    used_positions = []

    # Process each chromosome separately
    for chrom in sig_snps["CHR"].unique():
        chrom_snps = sig_snps[sig_snps["CHR"] == chrom].copy()
        chrom_used = []

        for _, snp in chrom_snps.iterrows():
            pos = snp["BP"]

            # Check if this position conflicts with already selected SNPs
            conflict = False
            for used_pos in chrom_used:
                if abs(pos - used_pos) < distance_threshold:
                    conflict = True
                    break

            if not conflict:
                # Define locus boundaries
                start_pos = max(1, pos - distance_threshold // 2)
                end_pos = pos + distance_threshold // 2

                # Count variants in this region
                region_variants = sumstats[
                    (sumstats["CHR"] == chrom)
                    & (sumstats["BP"] >= start_pos)
                    & (sumstats["BP"] <= end_pos)
                ]

                if len(region_variants) >= min_variants_per_locus:
                    independent_loci.append(
                        {
                            "chr": chrom,
                            "start": start_pos,
                            "end": end_pos,
                            "lead_snp": snp["SNPID"],
                            "lead_bp": pos,
                            "lead_p": snp["P"],
                            "ancestry": ancestry,
                            "n_variants": len(region_variants),
                        }
                    )
                    chrom_used.append(pos)
                    used_positions.append((chrom, pos))

    logger.info(f"{ancestry}: Identified {len(independent_loci)} independent loci")
    return independent_loci


def _merge_overlapping_loci(loci_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping loci across ancestries.

    Adapted from easyfinemap loci merging approach.
    """
    merged_loci = []
    processed_indices = set()

    for i, locus in loci_df.iterrows():
        if i in processed_indices:
            continue

        # Find overlapping loci
        overlapping = loci_df[
            (loci_df["chr"] == locus["chr"])
            & (loci_df.index != i)
            & ~((loci_df["end"] < locus["start"]) | (loci_df["start"] > locus["end"]))
        ]

        if len(overlapping) == 0:
            # No overlap, keep as is
            merged_loci.append(locus.to_dict())
            processed_indices.add(i)
        else:
            # Merge overlapping loci
            all_loci = pd.concat([locus.to_frame().T, overlapping])

            merged_locus = {
                "chr": locus["chr"],
                "start": all_loci["start"].min(),
                "end": all_loci["end"].max(),
                "lead_snp": all_loci.loc[all_loci["lead_p"].idxmin(), "lead_snp"],
                "lead_bp": all_loci.loc[all_loci["lead_p"].idxmin(), "lead_bp"],
                "lead_p": all_loci["lead_p"].min(),
                "ancestry": ",".join(sorted(all_loci["ancestry"].unique())),
                "n_variants": all_loci["n_variants"].max(),
            }

            merged_loci.append(merged_locus)
            processed_indices.add(i)
            processed_indices.update(overlapping.index)

    return pd.DataFrame(merged_loci)


def chunk_sumstats(
    loci_df: pd.DataFrame,
    sumstats_files: Dict[str, str],
    output_dir: str,
    threads: int = 1,
    compress: bool = True,
) -> pd.DataFrame:
    """
    Chunk summary statistics files into loci-specific files.

    Parameters
    ----------
    loci_df : pd.DataFrame
        DataFrame with loci coordinates from identify_independent_loci.
    sumstats_files : Dict[str, str]
        Dictionary mapping ancestry names to sumstats file paths.
    output_dir : str
        Output directory for chunked files.
    threads : int, optional
        Number of threads for parallel processing, by default 1.
    compress : bool, optional
        Whether to compress output files, by default True.

    Returns
    -------
    pd.DataFrame
        DataFrame with information about generated files.

    Examples
    --------
    >>> loci_df = identify_independent_loci(files, "output/")
    >>> file_info = chunk_sumstats(loci_df, files, "output/chunks/")
    """
    os.makedirs(output_dir, exist_ok=True)

    file_info_list = []

    # Load all sumstats files
    sumstats_data = {}
    for ancestry, file_path in sumstats_files.items():
        logger.info(f"Loading {ancestry} sumstats: {file_path}")
        sumstats_data[ancestry] = pd.read_csv(file_path, sep="\t", compression="gzip")

    # Process each locus
    for _, locus in tqdm(loci_df.iterrows(), total=len(loci_df), desc="Chunking loci"):
        locus_id = locus["locus_id"]
        chrom = locus["chr"]
        start = locus["start"]
        end = locus["end"]

        # Get ancestries for this locus
        if "," in locus["ancestry"]:
            ancestries = locus["ancestry"].split(",")
        else:
            ancestries = [locus["ancestry"]]

        # Extract data for each ancestry
        for ancestry in ancestries:
            if ancestry not in sumstats_data:
                logger.warning(f"No data available for ancestry: {ancestry}")
                continue

            # Extract locus data
            locus_data = sumstats_data[ancestry][
                (sumstats_data[ancestry]["CHR"] == chrom)
                & (sumstats_data[ancestry]["BP"] >= start)
                & (sumstats_data[ancestry]["BP"] <= end)
            ].copy()

            if len(locus_data) == 0:
                logger.warning(f"No variants found for {ancestry} {locus_id}")
                continue

            # Define output file
            suffix = ".gz" if compress else ""
            output_file = os.path.join(
                output_dir, f"{ancestry}.{locus_id}.sumstats{suffix}"
            )

            # Save chunked data
            if compress:
                locus_data.to_csv(
                    output_file, sep="\t", index=False, compression="gzip"
                )
            else:
                locus_data.to_csv(output_file, sep="\t", index=False)

            # Record file info
            file_info_list.append(
                {
                    "locus_id": locus_id,
                    "ancestry": ancestry,
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "n_variants": len(locus_data),
                    "sumstats_file": output_file,
                }
            )

    # Create file info DataFrame
    file_info_df = pd.DataFrame(file_info_list)

    # Save file info
    info_file = os.path.join(output_dir, "chunk_info.txt")
    file_info_df.to_csv(info_file, sep="\t", index=False)

    logger.info(f"Generated {len(file_info_list)} chunked files")
    logger.info(f"File information saved to {info_file}")

    return file_info_df


def create_loci_list_for_credtools(
    chunk_info_df: pd.DataFrame,
    ld_info_df: Optional[pd.DataFrame] = None,
    output_file: str = "loci_list.txt",
) -> pd.DataFrame:
    """
    Create loci list file compatible with credtools format.

    Parameters
    ----------
    chunk_info_df : pd.DataFrame
        DataFrame from chunk_sumstats with file information.
    ld_info_df : Optional[pd.DataFrame], optional
        DataFrame with LD file information, by default None.
    output_file : str, optional
        Output file path, by default "loci_list.txt".

    Returns
    -------
    pd.DataFrame
        DataFrame in credtools loci list format.
    """
    # Group by locus_id to create credtools format
    loci_list = []

    for locus_id, group in chunk_info_df.groupby("locus_id"):
        for _, row in group.iterrows():
            # Extract prefix from sumstats file
            sumstats_file = row["sumstats_file"]
            prefix = str(Path(sumstats_file).with_suffix("")).replace(".sumstats", "")

            locus_entry = {
                "locus_id": locus_id,
                "chr": row["chr"],
                "start": row["start"],
                "end": row["end"],
                "popu": row["ancestry"],
                "cohort": row["ancestry"],  # Use ancestry as cohort for now
                "sample_size": 50000,  # Placeholder - should be provided by user
                "prefix": prefix,
            }

            # Add LD file info if available
            if ld_info_df is not None:
                ld_match = ld_info_df[
                    (ld_info_df["locus_id"] == locus_id)
                    & (ld_info_df["ancestry"] == row["ancestry"])
                ]
                if len(ld_match) > 0:
                    locus_entry.update(ld_match.iloc[0].to_dict())

            loci_list.append(locus_entry)

    loci_df = pd.DataFrame(loci_list)

    # Save to file
    loci_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Created credtools loci list: {output_file}")

    return loci_df
