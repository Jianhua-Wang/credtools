"""Prepare LD matrices and final fine-mapping inputs.

This module provides functionality to extract LD matrices from genotype data
and create final input files compatible with credtools fine-mapping pipeline.
"""

import logging
import os
import subprocess
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from tqdm import tqdm

from ..constants import ColName
from ..ldmatrix import LDMatrix
from ..sumstats import make_SNPID_unique, munge

logger = logging.getLogger("Preprocessing")


def prepare_finemap_inputs(
    chunk_info_df: pd.DataFrame,
    genotype_files: Dict[str, str],
    output_dir: str,
    threads: int = 1,
    ld_format: str = "plink",
    keep_intermediate: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """
    Prepare final fine-mapping input files from chunked sumstats and genotype data.

    Parameters
    ----------
    chunk_info_df : pd.DataFrame
        DataFrame from chunk_sumstats with chunked file information.
    genotype_files : Dict[str, str]
        Dictionary mapping ancestry names to genotype file prefixes.
        Supports PLINK format (.bed/.bim/.fam) and VCF format.
    output_dir : str
        Output directory for prepared files.
    threads : int, optional
        Number of threads for parallel processing, by default 1.
    ld_format : str, optional
        Format for LD computation ("plink", "vcf"), by default "plink".
    keep_intermediate : bool, optional
        Whether to keep intermediate files, by default False.
    **kwargs
        Additional parameters.

    Returns
    -------
    pd.DataFrame
        DataFrame with information about prepared files.

    Examples
    --------
    >>> genotype_files = {"EUR": "eur_genotypes", "ASN": "asn_genotypes"}
    >>> prepared_df = prepare_finemap_inputs(chunk_info_df, genotype_files, "output/")
    """
    os.makedirs(output_dir, exist_ok=True)

    # Group by ancestry for parallel processing
    ancestry_groups = chunk_info_df.groupby("ancestry")

    # Prepare arguments for parallel processing
    prepare_args = []
    for ancestry, group in ancestry_groups:
        if ancestry not in genotype_files:
            logger.warning(f"No genotype file specified for ancestry: {ancestry}")
            continue

        prepare_args.append(
            (
                ancestry,
                genotype_files[ancestry],
                group,
                output_dir,
                ld_format,
                keep_intermediate,
                kwargs,
            )
        )

    # Process ancestries in parallel
    if threads > 1:
        logger.info(
            f"Processing {len(prepare_args)} ancestries using {threads} threads"
        )
        with Pool(threads) as pool:
            results = list(
                tqdm(
                    pool.starmap(_prepare_ancestry_files, prepare_args),
                    total=len(prepare_args),
                    desc="Preparing ancestries",
                )
            )
    else:
        logger.info(f"Processing {len(prepare_args)} ancestries sequentially")
        results = [
            _prepare_ancestry_files(*args)
            for args in tqdm(prepare_args, desc="Preparing ancestries")
        ]

    # Combine results
    all_prepared_files = []
    for result in results:
        all_prepared_files.extend(result)

    prepared_df = pd.DataFrame(all_prepared_files)

    # Save results
    output_file = os.path.join(output_dir, "prepared_files.txt")
    prepared_df.to_csv(output_file, sep="\t", index=False)

    logger.info(f"Prepared {len(all_prepared_files)} locus files")
    logger.info(f"Results saved to {output_file}")

    return prepared_df


def _prepare_ancestry_files(
    ancestry: str,
    genotype_prefix: str,
    chunk_group: pd.DataFrame,
    output_dir: str,
    ld_format: str,
    keep_intermediate: bool,
    kwargs: dict,
) -> List[Dict]:
    """Prepare files for a single ancestry."""
    logger.info(f"Processing ancestry: {ancestry}")

    prepared_files = []

    for _, chunk_info in chunk_group.iterrows():
        try:
            result = _prepare_single_locus(
                chunk_info,
                ancestry,
                genotype_prefix,
                output_dir,
                ld_format,
                keep_intermediate,
                **kwargs,
            )
            if result:
                prepared_files.append(result)

        except Exception as e:
            logger.error(
                f"Failed to process {chunk_info['locus_id']} for {ancestry}: {str(e)}"
            )
            continue

    return prepared_files


def _prepare_single_locus(
    chunk_info: pd.Series,
    ancestry: str,
    genotype_prefix: str,
    output_dir: str,
    ld_format: str,
    keep_intermediate: bool,
    **kwargs,
) -> Optional[Dict]:
    """
    Prepare files for a single locus.

    Adapted from preprocessing.py workflow.
    """
    locus_id = chunk_info["locus_id"]
    chrom = chunk_info["chr"]
    start = chunk_info["start"]
    end = chunk_info["end"]
    sumstats_file = chunk_info["sumstats_file"]

    # Define output prefix
    output_prefix = os.path.join(output_dir, f"{ancestry}.{locus_id}")

    # Check if output files already exist
    expected_files = [
        f"{output_prefix}.sumstats.gz",
        f"{output_prefix}.ld.npz",
        f"{output_prefix}.ldmap.gz",
    ]

    if all(os.path.exists(f) for f in expected_files):
        logger.info(f"Output files exist for {ancestry} {locus_id}, skipping")
        return {
            "locus_id": locus_id,
            "ancestry": ancestry,
            "chr": chrom,
            "start": start,
            "end": end,
            "prefix": output_prefix,
            "status": "existed",
        }

    try:
        # Load and process sumstats
        logger.debug(f"Processing sumstats for {ancestry} {locus_id}")
        sumstats = pd.read_csv(sumstats_file, sep="\t", compression="gzip")

        # Munge sumstats to ensure proper format
        sumstats = munge(sumstats)

        # Extract and process LD matrix
        logger.debug(f"Extracting LD matrix for {ancestry} {locus_id}")
        ld_result = _extract_ld_matrix(
            genotype_prefix=genotype_prefix,
            chrom=chrom,
            start=start,
            end=end,
            output_prefix=output_prefix,
            ld_format=ld_format,
            keep_intermediate=keep_intermediate,
        )

        if not ld_result:
            logger.warning(f"Failed to extract LD matrix for {ancestry} {locus_id}")
            return None

        ldmap, ld_matrix = ld_result

        # Intersect sumstats and LD data
        logger.debug(f"Intersecting data for {ancestry} {locus_id}")
        intersected_sumstats, intersected_ld, intersected_ldmap = (
            _intersect_sumstats_ld(sumstats, ld_matrix, ldmap)
        )

        if len(intersected_sumstats) == 0:
            logger.warning(f"No common variants found for {ancestry} {locus_id}")
            return None

        # Save final files
        logger.debug(f"Saving final files for {ancestry} {locus_id}")

        # Save sumstats
        intersected_sumstats.to_csv(f"{output_prefix}.sumstats", sep="\t", index=False)
        subprocess.run(f"gzip -f {output_prefix}.sumstats", shell=True, check=True)

        # Save LD matrix
        np.savez_compressed(
            f"{output_prefix}.ld.npz", intersected_ld.astype(np.float16)
        )

        # Save LD map
        intersected_ldmap.to_csv(f"{output_prefix}.ldmap", sep="\t", index=False)
        subprocess.run(f"gzip -f {output_prefix}.ldmap", shell=True, check=True)

        return {
            "locus_id": locus_id,
            "ancestry": ancestry,
            "chr": chrom,
            "start": start,
            "end": end,
            "prefix": output_prefix,
            "n_variants": len(intersected_sumstats),
            "status": "created",
        }

    except Exception as e:
        logger.error(f"Error processing {ancestry} {locus_id}: {str(e)}")
        # Clean up partial files
        for expected_file in expected_files:
            if os.path.exists(expected_file):
                os.remove(expected_file)
        return None


def _extract_ld_matrix(
    genotype_prefix: str,
    chrom: int,
    start: int,
    end: int,
    output_prefix: str,
    ld_format: str,
    keep_intermediate: bool,
) -> Optional[Tuple[pd.DataFrame, np.ndarray]]:
    """Extract LD matrix from genotype data for a genomic region."""
    if ld_format.lower() == "plink":
        return _extract_ld_plink(
            genotype_prefix, chrom, start, end, output_prefix, keep_intermediate
        )
    elif ld_format.lower() == "vcf":
        return _extract_ld_vcf(
            genotype_prefix, chrom, start, end, output_prefix, keep_intermediate
        )
    else:
        raise ValueError(f"Unsupported LD format: {ld_format}")


def _extract_ld_plink(
    genotype_prefix: str,
    chrom: int,
    start: int,
    end: int,
    output_prefix: str,
    keep_intermediate: bool,
) -> Optional[Tuple[pd.DataFrame, np.ndarray]]:
    """Extract LD matrix using PLINK format files."""
    try:
        # Use PLINK to compute LD matrix for the region
        temp_prefix = f"{output_prefix}_temp"

        # Extract variants in region
        plink_cmd = [
            "plink",
            "--bfile",
            genotype_prefix,
            "--chr",
            str(chrom),
            "--from-bp",
            str(start),
            "--to-bp",
            str(end),
            "--make-bed",
            "--out",
            temp_prefix,
            "--silent",
        ]

        result = subprocess.run(plink_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"PLINK extraction failed: {result.stderr}")
            return None

        # Check if any variants were extracted
        if not os.path.exists(f"{temp_prefix}.bim"):
            logger.warning(f"No variants found in region chr{chrom}:{start}-{end}")
            return None

        # Read variant information
        bim_df = pd.read_csv(
            f"{temp_prefix}.bim",
            sep="\t",
            header=None,
            names=["CHR", "SNPID", "CM", "BP", "A1", "A2"],
        )

        if len(bim_df) < 2:
            logger.warning(f"Insufficient variants for LD computation: {len(bim_df)}")
            return None

        # Compute LD matrix
        ld_cmd = [
            "plink",
            "--bfile",
            temp_prefix,
            "--r",
            "square",
            "--out",
            temp_prefix,
            "--silent",
        ]

        result = subprocess.run(ld_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"PLINK LD computation failed: {result.stderr}")
            return None

        # Load LD matrix
        ld_file = f"{temp_prefix}.ld"
        if not os.path.exists(ld_file):
            logger.error(f"LD matrix file not found: {ld_file}")
            return None

        ld_matrix = np.loadtxt(ld_file)

        # Prepare LD map
        ldmap = bim_df[["CHR", "BP", "A1", "A2"]].copy()
        ldmap["AF2"] = 0.5  # Placeholder - should compute from genotypes

        # Clean up intermediate files
        if not keep_intermediate:
            temp_files = [
                f"{temp_prefix}.bed",
                f"{temp_prefix}.bim",
                f"{temp_prefix}.fam",
                f"{temp_prefix}.ld",
                f"{temp_prefix}.log",
            ]
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    os.remove(temp_file)

        return ldmap, ld_matrix

    except Exception as e:
        logger.error(f"Error extracting LD matrix: {str(e)}")
        return None


def _extract_ld_vcf(
    genotype_prefix: str,
    chrom: int,
    start: int,
    end: int,
    output_prefix: str,
    keep_intermediate: bool,
) -> Optional[Tuple[pd.DataFrame, np.ndarray]]:
    """Extract LD matrix from VCF format files."""
    # This would require vcftools or similar
    # For now, raise NotImplementedError
    raise NotImplementedError("VCF format LD extraction not yet implemented")


def _intersect_sumstats_ld(
    sumstats: pd.DataFrame, ld_matrix: np.ndarray, ldmap: pd.DataFrame
) -> Tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    """
    Intersect summary statistics and LD matrix data.

    Adapted from existing credtools intersection logic.
    """
    # Make SNPID unique for both datasets
    sumstats = make_SNPID_unique(sumstats)
    ldmap = make_SNPID_unique(ldmap, col_ea="A1", col_nea="A2")

    # Find common variants
    common_snpids = set(sumstats[ColName.SNPID]) & set(ldmap[ColName.SNPID])

    if len(common_snpids) == 0:
        logger.error("No common variants between sumstats and LD data")
        return pd.DataFrame(), np.array([]), pd.DataFrame()

    # Subset sumstats
    intersected_sumstats = (
        sumstats[sumstats[ColName.SNPID].isin(common_snpids)]
        .copy()
        .sort_values([ColName.CHR, ColName.BP])
        .reset_index(drop=True)
    )

    # Subset LD data maintaining order
    ldmap_indexed = ldmap.copy()
    ldmap_indexed["original_index"] = ldmap_indexed.index
    ldmap_indexed = ldmap_indexed.set_index(ColName.SNPID)

    # Get ordered common SNPs
    ordered_snpids = intersected_sumstats[ColName.SNPID].tolist()
    intersected_ldmap = ldmap_indexed.loc[ordered_snpids].copy()
    original_indices = intersected_ldmap["original_index"].values

    # Subset LD matrix
    intersected_ld = ld_matrix[np.ix_(original_indices, original_indices)]  # type: ignore

    # Clean up ldmap
    intersected_ldmap = intersected_ldmap.drop("original_index", axis=1).reset_index(
        drop=True
    )

    # Ensure diagonal is 1
    np.fill_diagonal(intersected_ld, 1.0)

    # Handle allele flipping if needed
    intersected_ldmap, intersected_ld = _handle_allele_flipping(
        intersected_ldmap, intersected_ld
    )

    logger.info(f"Intersected {len(intersected_sumstats)} variants")

    return intersected_sumstats, intersected_ld, intersected_ldmap


def _handle_allele_flipping(
    ldmap: pd.DataFrame, ld_matrix: np.ndarray
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Handle allele flipping to ensure consistent allele ordering.

    Adapted from preprocessing.py logic.
    """
    # Sort alleles alphabetically
    ldmap[["sort_a1", "sort_a2"]] = np.sort(ldmap[["A1", "A2"]], axis=1)

    # Find indices where alleles were swapped
    swapped_indices = ldmap[ldmap["A1"] != ldmap["sort_a1"]].index

    # Flip LD matrix values for swapped alleles
    if len(swapped_indices) > 0:
        ld_matrix[swapped_indices] *= -1
        ld_matrix[:, swapped_indices] *= -1

    # Update allele columns
    ldmap["A1"] = ldmap["sort_a1"]
    ldmap["A2"] = ldmap["sort_a2"]
    ldmap = ldmap.drop(columns=["sort_a1", "sort_a2"])

    # Ensure diagonal is still 1
    np.fill_diagonal(ld_matrix, 1.0)

    return ldmap, ld_matrix
