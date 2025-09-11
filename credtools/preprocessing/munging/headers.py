"""
Header detection and mapping utilities for GWAS summary statistics.

Adapted from smunger (https://github.com/Jianhua-Wang/smunger)
Original author: Jianhua Wang
License: MIT

This module provides functionality to detect, map, and standardize column
headers from various GWAS summary statistics file formats.
"""

import logging
from typing import Any, Dict, List, Optional

import pandas as pd

from .constants import COMMON_COLNAMES, ColName, suggest_column_mapping

logger = logging.getLogger("Munging")


def inspect_headers(
    file_path: str, sep: Optional[str] = None, nrows: int = 5
) -> List[str]:
    """
    Inspect file headers and return column names.

    Parameters
    ----------
    file_path : str
        Path to the input file.
    sep : str, optional
        Column separator. If None, will try to auto-detect.
    nrows : int, optional
        Number of rows to read for inspection.

    Returns
    -------
    List[str]
        List of column headers.
    """
    # Auto-detect separator if not provided
    if sep is None:
        sep = _detect_separator(file_path)

    try:
        # Read just the header and a few rows
        df = pd.read_csv(file_path, sep=sep, nrows=nrows)
        headers = df.columns.tolist()

        logger.info(f"Detected {len(headers)} columns in {file_path}")
        logger.debug(f"Headers: {headers}")

        return headers

    except Exception as e:
        logger.error(f"Failed to inspect headers in {file_path}: {str(e)}")
        raise


def _detect_separator(file_path: str) -> str:
    """
    Auto-detect column separator in file.

    Parameters
    ----------
    file_path : str
        Path to input file.

    Returns
    -------
    str
        Detected separator character.
    """
    import gzip

    # Determine if file is gzipped
    is_gzipped = file_path.endswith(".gz")

    try:
        # Read first line
        if is_gzipped:
            with gzip.open(file_path, "rt") as f:
                first_line = f.readline().strip()
        else:
            with open(file_path, "r") as f:
                first_line = f.readline().strip()

        # Count separators
        tab_count = first_line.count("\t")
        comma_count = first_line.count(",")
        space_count = len(first_line.split()) - 1

        # Choose most frequent separator
        if tab_count > comma_count and tab_count > space_count:
            sep = "\t"
        elif comma_count > space_count:
            sep = ","
        else:
            sep = r"\s+"  # Multiple whitespace

        logger.debug(f"Auto-detected separator: '{sep}'")
        return sep

    except Exception as e:
        logger.warning(f"Failed to auto-detect separator: {e}. Using tab.")
        return "\t"


def map_headers_automatic(headers: List[str]) -> Dict[str, str]:
    """
    Automatically map headers to standard column names.

    Parameters
    ----------
    headers : List[str]
        List of column headers from input file.

    Returns
    -------
    Dict[str, str]
        Mapping from original headers to standard column names.
    """
    mapping = {}

    for header in headers:
        if header in COMMON_COLNAMES:
            mapping[header] = COMMON_COLNAMES[header]
        else:
            # Try fuzzy matching
            mapped = _fuzzy_match_header(header)
            mapping[header] = mapped or header

    logger.info(
        f"Automatically mapped {sum(1 for v in mapping.values() if v in ColName.sumstat_cols)} columns"
    )

    return mapping


def _fuzzy_match_header(header: str) -> Optional[str]:
    """
    Try to match header using fuzzy string matching.

    Parameters
    ----------
    header : str
        Header to match.

    Returns
    -------
    Optional[str]
        Matched standard column name, or None if no match.
    """
    header_lower = header.lower()

    # Pattern matching for common variations
    patterns = {
        ColName.CHR: ["chr", "chrom", "chromosome"],
        ColName.BP: ["bp", "pos", "position", "coordinate"],
        ColName.SNPID: ["snp", "variant", "marker", "id"],
        ColName.EA: ["a1", "ea", "effect", "alt", "allele1"],
        ColName.NEA: ["a2", "nea", "other", "ref", "allele2"],
        ColName.BETA: ["beta", "effect", "coef"],
        ColName.SE: ["se", "stderr", "error"],
        ColName.P: ["p", "pval", "pvalue", "p_val"],
        ColName.EAF: ["eaf", "freq", "maf", "af"],
        ColName.OR: ["or", "odds"],
        ColName.N: ["n", "sample", "size"],
        ColName.INFO: ["info", "imputation", "quality"],
        ColName.Z: ["z", "zscore", "stat"],
        ColName.RSID: ["rs", "rsid", "dbsnp"],
    }

    for standard_col, pattern_list in patterns.items():
        for pattern in pattern_list:
            if pattern in header_lower:
                logger.debug(
                    f"Fuzzy matched '{header}' -> '{standard_col}' (pattern: '{pattern}')"
                )
                return standard_col

    return None


def create_config_template(
    headers: List[str], interactive: bool = False
) -> Dict[str, Any]:
    """
    Create configuration template for column mapping.

    Parameters
    ----------
    headers : List[str]
        List of column headers.
    interactive : bool, optional
        Whether to use interactive mode for mapping.

    Returns
    -------
    Dict[str, Any]
        Configuration dictionary with column mappings.
    """
    if interactive:
        return _create_interactive_config(headers)
    else:
        return _create_automatic_config(headers)


def _create_automatic_config(headers: List[str]) -> Dict[str, Any]:
    """
    Create configuration using automatic header mapping.

    Parameters
    ----------
    headers : List[str]
        List of column headers.

    Returns
    -------
    Dict[str, Any]
        Configuration dictionary.
    """
    mapping = map_headers_automatic(headers)
    suggestions = suggest_column_mapping(headers)

    config = {"column_mapping": mapping, "suggestions": suggestions, "headers": headers}

    return config


def _create_interactive_config(headers: List[str]) -> Dict[str, Any]:
    """
    Create configuration using interactive prompts.

    Note: This is a placeholder for interactive functionality.
    In a real CLI environment, this would prompt the user for input.

    Parameters
    ----------
    headers : List[str]
        List of column headers.

    Returns
    -------
    Dict[str, Any]
        Configuration dictionary.
    """
    logger.info("Interactive configuration not fully implemented")
    logger.info("Falling back to automatic configuration")

    return _create_automatic_config(headers)


def apply_header_mapping(df: pd.DataFrame, mapping: Dict[str, str]) -> pd.DataFrame:
    """
    Apply header mapping to DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    mapping : Dict[str, str]
        Mapping from current column names to new names.

    Returns
    -------
    pd.DataFrame
        DataFrame with renamed columns.
    """
    outdf = df.copy()

    # Only rename columns that exist in the DataFrame
    existing_mapping = {k: v for k, v in mapping.items() if k in outdf.columns}

    if existing_mapping:
        outdf = outdf.rename(columns=existing_mapping)
        logger.info(f"Applied mapping to {len(existing_mapping)} columns")

    return outdf


def validate_required_columns(
    df: pd.DataFrame, required: Optional[List[str]] = None
) -> bool:
    """
    Validate that required columns are present.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to validate.
    required : List[str], optional
        List of required column names. Uses mandatory columns if None.

    Returns
    -------
    bool
        True if all required columns are present.
    """
    if required is None:
        required = ColName.mandatory_cols

    missing = set(required) - set(df.columns)

    if missing:
        logger.error(f"Missing required columns: {missing}")
        return False

    logger.info("All required columns are present")
    return True


def suggest_missing_mappings(
    headers: List[str], mapped_headers: Dict[str, str]
) -> Dict[str, str]:
    """
    Suggest mappings for unmapped headers.

    Parameters
    ----------
    headers : List[str]
        Original headers.
    mapped_headers : Dict[str, str]
        Already mapped headers.

    Returns
    -------
    Dict[str, str]
        Suggestions for unmapped headers.
    """
    mapped_standards = set(mapped_headers.values())
    required_standards = set(ColName.mandatory_cols)

    missing_standards = required_standards - mapped_standards
    unmapped_headers = [h for h in headers if h not in mapped_headers]

    suggestions = {}

    # Try to map missing required columns
    for missing_std in missing_standards:
        for unmapped in unmapped_headers:
            fuzzy_match = _fuzzy_match_header(unmapped)
            if fuzzy_match == missing_std:
                suggestions[unmapped] = missing_std
                unmapped_headers.remove(unmapped)
                break

    logger.info(f"Generated {len(suggestions)} additional mapping suggestions")

    return suggestions
