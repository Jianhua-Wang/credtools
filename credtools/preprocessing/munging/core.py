"""
Core munging functions for GWAS summary statistics.

Adapted from smunger (https://github.com/Jianhua-Wang/smunger)
Original author: Jianhua Wang
License: MIT

This module provides the main data cleaning and standardization functions
for processing GWAS summary statistics.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

from .constants import ColAllowNA, ColName, ColRange, ColType
from .validation import check_mandatory_cols, validate_and_clean_column

logger = logging.getLogger("Munging")


def munge(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and standardize GWAS summary statistics.

    Adapted from smunger.munge() function.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with GWAS summary statistics.

    Returns
    -------
    pd.DataFrame
        Cleaned and standardized DataFrame.

    Notes
    -----
    This function performs comprehensive data cleaning including:
    1. Removing columns with all NA values
    2. Cleaning and validating core columns (CHR, BP, alleles)
    3. Creating unique SNP identifiers
    4. Processing p-values, effect sizes, and other statistics
    5. Sorting by chromosome and position
    """
    logger.info("Starting munging process")
    original_rows = len(df)

    # Check mandatory columns
    check_mandatory_cols(df)

    # Make a copy to avoid modifying input
    outdf = df.copy()

    # Remove columns that are all NA
    outdf = _remove_all_na_columns(outdf)

    # Clean core columns
    outdf = _munge_chr(outdf)
    outdf = _munge_bp(outdf)
    outdf = _munge_alleles(outdf)

    # Create unique SNP identifiers
    outdf = make_SNPID_unique(outdf)

    # Process statistical columns
    outdf = _munge_pvalue(outdf)
    outdf = _munge_beta(outdf)
    outdf = _munge_se(outdf)

    # Process frequency columns if present
    if ColName.EAF in outdf.columns:
        outdf = _munge_eaf(outdf)
        outdf[ColName.MAF] = outdf[ColName.EAF].apply(
            lambda x: min(x, 1 - x) if pd.notna(x) else x
        )
        outdf = _munge_maf(outdf)

    # Sort by chromosome and position
    outdf = outdf.sort_values([ColName.CHR, ColName.BP])
    outdf.reset_index(drop=True, inplace=True)

    # Ensure correct column order and types
    outdf = _finalize_columns(outdf)

    final_rows = len(outdf)
    logger.info(
        f"Munging complete: {original_rows} -> {final_rows} rows ({original_rows - final_rows} removed)"
    )

    return outdf


def make_SNPID_unique(
    df: pd.DataFrame,
    remove_duplicates: bool = True,
    col_chr: str = ColName.CHR,
    col_bp: str = ColName.BP,
    col_ea: str = ColName.EA,
    col_nea: str = ColName.NEA,
    col_p: str = ColName.P,
) -> pd.DataFrame:
    """
    Generate unique SNP identifiers and optionally remove duplicates.

    Adapted from smunger.make_SNPID_unique() function.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    remove_duplicates : bool, optional
        Whether to remove duplicated SNPs, keeping the one with smallest p-value.
    col_chr : str, optional
        Column name for chromosome.
    col_bp : str, optional
        Column name for base pair position.
    col_ea : str, optional
        Column name for effect allele.
    col_nea : str, optional
        Column name for non-effect allele.
    col_p : str, optional
        Column name for p-value.

    Returns
    -------
    pd.DataFrame
        DataFrame with unique SNPID column.
    """
    outdf = df.copy()

    # Sort alleles alphabetically to ensure consistent SNP IDs
    allele_df = outdf[[col_ea, col_nea]].apply(
        lambda row: sorted([str(row[col_ea]), str(row[col_nea])]),
        axis=1,
        result_type="expand",
    )
    allele_df.columns = ["allele1", "allele2"]

    # Create unique SNPID: chr-bp-allele1-allele2
    outdf[ColName.SNPID] = (
        outdf[col_chr].astype(str)
        + "-"
        + outdf[col_bp].astype(str)
        + "-"
        + allele_df["allele1"]
        + "-"
        + allele_df["allele2"]
    )

    # Move SNPID to first column
    cols = outdf.columns.tolist()
    if ColName.SNPID in cols:
        cols.insert(0, cols.pop(cols.index(ColName.SNPID)))
        outdf = outdf[cols]

    # Handle duplicates
    n_duplicated = outdf.duplicated(subset=[ColName.SNPID]).sum()

    if remove_duplicates and n_duplicated > 0:
        logger.info(f"Removing {n_duplicated} duplicate SNPs")
        if col_p in outdf.columns:
            # Sort by p-value to keep most significant
            outdf = outdf.sort_values(col_p)
        outdf = outdf.drop_duplicates(subset=[ColName.SNPID], keep="first")
        outdf = outdf.sort_values([col_chr, col_bp])
        outdf.reset_index(drop=True, inplace=True)
    elif n_duplicated > 0:
        logger.warning(f"Found {n_duplicated} duplicate SNPs, keeping all")
        # Add suffix to make unique
        dup_suffix = "-" + outdf.groupby(ColName.SNPID).cumcount().astype(str)
        dup_suffix = dup_suffix.str.replace("-0", "")
        outdf[ColName.SNPID] = outdf[ColName.SNPID] + dup_suffix

    logger.debug(f"Created unique SNPIDs: {len(outdf)} variants")

    return outdf


def _remove_all_na_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Remove columns that are entirely NA."""
    outdf = df.copy()
    outdf = outdf.replace("", np.nan)  # Convert empty strings to NaN

    for col in outdf.columns:
        if outdf[col].isna().all():
            logger.debug(f"Removing all-NA column: {col}")
            outdf = outdf.drop(columns=[col])

    return outdf


def _munge_chr(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate chromosome column."""
    return validate_and_clean_column(
        df=df,
        col_name=ColName.CHR,
        col_type=ColType.CHR,
        min_val=ColRange.CHR_MIN,
        max_val=ColRange.CHR_MAX,
        allow_na=ColAllowNA.CHR,
        transform_func=_transform_chr,
    )


def _transform_chr(series: pd.Series) -> pd.Series:
    """Transform chromosome values to standard format."""
    # Convert to string first
    result = series.astype(str)

    # Remove 'chr' prefix
    result = result.str.replace("chr", "", case=False)

    # Convert X to 23
    result = result.replace(["X", "x"], "23")

    # Convert to numeric
    result = pd.to_numeric(result, errors="coerce")

    return result


def _munge_bp(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate base pair position column."""
    return validate_and_clean_column(
        df=df,
        col_name=ColName.BP,
        col_type=ColType.BP,
        min_val=ColRange.BP_MIN,
        max_val=ColRange.BP_MAX,
        allow_na=ColAllowNA.BP,
    )


def _munge_alleles(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate allele columns."""
    outdf = df.copy()

    # Clean effect allele
    outdf = validate_and_clean_column(
        df=outdf,
        col_name=ColName.EA,
        col_type=ColType.EA,
        allow_na=ColAllowNA.EA,
        transform_func=_transform_allele,
    )

    # Clean non-effect allele
    outdf = validate_and_clean_column(
        df=outdf,
        col_name=ColName.NEA,
        col_type=ColType.NEA,
        allow_na=ColAllowNA.NEA,
        transform_func=_transform_allele,
    )

    # Remove rows where EA == NEA
    before_count = len(outdf)
    outdf = outdf[outdf[ColName.EA] != outdf[ColName.NEA]]
    after_count = len(outdf)

    if before_count > after_count:
        logger.debug(
            f"Removed {before_count - after_count} rows with identical alleles"
        )

    return outdf


def _transform_allele(series: pd.Series) -> pd.Series:
    """Transform allele values to standard format."""
    result = series.astype(str).str.upper()

    # Only keep valid DNA bases
    valid_pattern = r"^[ACGT]+$"
    mask = result.str.match(valid_pattern, na=False)
    result.loc[~mask] = np.nan

    return result


def _munge_pvalue(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate p-value column."""
    return validate_and_clean_column(
        df=df,
        col_name=ColName.P,
        col_type=ColType.P,
        min_val=ColRange.P_MIN,
        max_val=ColRange.P_MAX,
        allow_na=ColAllowNA.P,
        exclude_min=True,  # P-value cannot be exactly 0
    )


def _munge_beta(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate beta (effect size) column."""
    return validate_and_clean_column(
        df=df, col_name=ColName.BETA, col_type=ColType.BETA, allow_na=ColAllowNA.BETA
    )


def _munge_se(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate standard error column."""
    return validate_and_clean_column(
        df=df,
        col_name=ColName.SE,
        col_type=ColType.SE,
        min_val=ColRange.SE_MIN,
        allow_na=ColAllowNA.SE,
        exclude_min=True,  # SE cannot be exactly 0
    )


def _munge_eaf(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate effect allele frequency column."""
    return validate_and_clean_column(
        df=df,
        col_name=ColName.EAF,
        col_type=ColType.EAF,
        min_val=ColRange.EAF_MIN,
        max_val=ColRange.EAF_MAX,
        allow_na=ColAllowNA.EAF,
    )


def _munge_maf(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and validate minor allele frequency column."""
    return validate_and_clean_column(
        df=df,
        col_name=ColName.MAF,
        col_type=ColType.MAF,
        min_val=ColRange.MAF_MIN,
        max_val=ColRange.MAF_MAX,
        allow_na=ColAllowNA.MAF,
    )


def _finalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure correct column order and add missing columns."""
    outdf = df.copy()

    # Add missing optional columns with None values
    for col in ColName.sumstat_cols:
        if col not in outdf.columns:
            outdf[col] = None

    # Reorder columns to standard order
    outdf = outdf[ColName.sumstat_cols]

    return outdf
