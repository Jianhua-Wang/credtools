"""
Validation functions for GWAS summary statistics columns.

Adapted from smunger (https://github.com/Jianhua-Wang/smunger)
Original author: Jianhua Wang
License: MIT

This module provides validation and cleaning functions for individual columns
in GWAS summary statistics data.
"""

import logging
from typing import Any, Callable, Optional

import numpy as np
import pandas as pd

from .constants import ColName

logger = logging.getLogger("Munging")


def check_mandatory_cols(df: pd.DataFrame) -> None:
    """
    Check if DataFrame contains all mandatory columns.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to validate.

    Raises
    ------
    ValueError
        If any mandatory columns are missing.
    """
    missing_cols = set(ColName.mandatory_cols) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing mandatory columns: {missing_cols}")


def validate_and_clean_column(
    df: pd.DataFrame,
    col_name: str,
    col_type: Any,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    allow_na: bool = True,
    exclude_min: bool = False,
    exclude_max: bool = False,
    transform_func: Optional[Callable] = None,
) -> pd.DataFrame:
    """
    Validate and clean a single column in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    col_name : str
        Name of column to validate.
    col_type : type
        Target data type for the column.
    min_val : float, optional
        Minimum allowed value.
    max_val : float, optional
        Maximum allowed value.
    allow_na : bool, optional
        Whether NA values are allowed.
    exclude_min : bool, optional
        Whether to exclude the minimum value itself.
    exclude_max : bool, optional
        Whether to exclude the maximum value itself.
    transform_func : callable, optional
        Function to transform values before validation.

    Returns
    -------
    pd.DataFrame
        DataFrame with cleaned column.
    """
    if col_name not in df.columns:
        if not allow_na:
            raise ValueError(f"Required column '{col_name}' not found")
        return df

    outdf = df.copy()
    original_count = len(outdf)

    # Apply transformation if provided
    if transform_func is not None:
        outdf[col_name] = transform_func(outdf[col_name])

    # Remove rows with NA values if not allowed
    if not allow_na:
        outdf = outdf[outdf[col_name].notna()]

    # Convert to target type for numeric columns
    if col_type in [
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        np.float16,
        np.float32,
        np.float64,
    ]:
        outdf[col_name] = pd.to_numeric(outdf[col_name], errors="coerce")

        # Remove rows that couldn't be converted
        if not allow_na:
            outdf = outdf[outdf[col_name].notna()]

    # Apply range validation
    if min_val is not None:
        if exclude_min:
            mask = outdf[col_name] > min_val
        else:
            mask = outdf[col_name] >= min_val
        outdf = outdf[mask | outdf[col_name].isna()]

    if max_val is not None:
        if exclude_max:
            mask = outdf[col_name] < max_val
        else:
            mask = outdf[col_name] <= max_val
        outdf = outdf[mask | outdf[col_name].isna()]

    # Convert to final type
    if col_type == str:
        outdf[col_name] = outdf[col_name].astype(str)
        outdf.loc[outdf[col_name] == "nan", col_name] = np.nan
    else:
        outdf[col_name] = outdf[col_name].astype(col_type)

    final_count = len(outdf)
    if original_count > final_count:
        logger.debug(
            f"Column {col_name}: removed {original_count - final_count} invalid rows"
        )

    return outdf


def validate_allele_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate that alleles are consistent and biallelic.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with EA and NEA columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with consistent alleles.
    """
    outdf = df.copy()

    if ColName.EA not in outdf.columns or ColName.NEA not in outdf.columns:
        return outdf

    original_count = len(outdf)

    # Remove rows where alleles are identical
    outdf = outdf[outdf[ColName.EA] != outdf[ColName.NEA]]

    # Remove rows where either allele is missing (if required)
    outdf = outdf[outdf[ColName.EA].notna() & outdf[ColName.NEA].notna()]

    final_count = len(outdf)
    if original_count > final_count:
        logger.debug(
            f"Removed {original_count - final_count} rows with invalid alleles"
        )

    return outdf


def validate_statistical_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate statistical consistency (e.g., beta and SE relationship).

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with statistical columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with statistically consistent data.
    """
    outdf = df.copy()

    # Check if we can compute Z-score from BETA and SE
    if all(col in outdf.columns for col in [ColName.BETA, ColName.SE]):
        # Remove rows where SE is 0 or negative
        mask = (outdf[ColName.SE] > 0) | outdf[ColName.SE].isna()
        removed = len(outdf) - mask.sum()
        if removed > 0:
            logger.debug(f"Removed {removed} rows with invalid SE values")
        outdf = outdf[mask]

    return outdf


def validate_frequency_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate frequency column consistency.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with frequency columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with consistent frequency data.
    """
    outdf = df.copy()

    # Ensure MAF is actually minor (≤ 0.5)
    if ColName.MAF in outdf.columns:
        # Convert frequencies > 0.5 to 1 - frequency
        mask = outdf[ColName.MAF] > 0.5
        outdf.loc[mask, ColName.MAF] = 1 - outdf.loc[mask, ColName.MAF]

    # Ensure EAF and MAF are consistent
    if all(col in outdf.columns for col in [ColName.EAF, ColName.MAF]):
        # MAF should be min(EAF, 1-EAF)
        expected_maf = outdf[ColName.EAF].apply(
            lambda x: min(x, 1 - x) if pd.notna(x) else np.nan
        )

        # Check for large discrepancies
        discrepancy = abs(outdf[ColName.MAF] - expected_maf)
        large_discrepancy = discrepancy > 0.01  # 1% tolerance

        if large_discrepancy.any():
            n_discrepant = large_discrepancy.sum()
            logger.warning(f"Found {n_discrepant} rows with EAF/MAF discrepancies")
            # Use computed MAF
            outdf.loc[large_discrepancy, ColName.MAF] = expected_maf.loc[
                large_discrepancy
            ]

    return outdf


def validate_pvalue_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate p-value consistency with other statistics.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with p-value and other statistical columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with consistent p-values.
    """
    outdf = df.copy()

    # If we have BETA, SE, and P, check consistency
    if all(col in outdf.columns for col in [ColName.BETA, ColName.SE, ColName.P]):
        # Compute Z-score
        z_score = outdf[ColName.BETA] / outdf[ColName.SE]

        # Compute expected p-value (two-tailed test)
        try:
            from scipy.stats import norm

            expected_p = 2 * (1 - norm.cdf(abs(z_score)))
        except ImportError:
            # Fallback if scipy not available
            logger.warning("scipy not available, skipping p-value consistency check")
            return outdf

        # Check for large discrepancies (order of magnitude)
        log_ratio = np.log10(outdf[ColName.P]) - np.log10(expected_p)
        large_discrepancy = abs(log_ratio) > 1  # 10-fold difference

        if large_discrepancy.any():
            n_discrepant = large_discrepancy.sum()
            logger.warning(f"Found {n_discrepant} rows with P-value discrepancies")

    return outdf
