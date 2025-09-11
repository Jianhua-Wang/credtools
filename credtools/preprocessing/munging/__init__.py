"""
Internal munging module for credtools.

Adapted from smunger (https://github.com/Jianhua-Wang/smunger)
Original author: Jianhua Wang
License: MIT

This module provides GWAS summary statistics munging functionality
without requiring external dependencies.

Main Functions
--------------
munge : Main data munging and cleaning function
read_config : Read configuration from file
inspect_headers : Examine file headers
create_config_template : Generate column mapping configuration
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .constants import ColName
from .core import make_SNPID_unique, munge
from .headers import (
    apply_header_mapping,
    create_config_template,
    inspect_headers,
    map_headers_automatic,
    validate_required_columns,
)

logger = logging.getLogger("Munging")

# Main API functions (compatible with smunger)
__all__ = [
    "munge",
    "read_config",
    "inspect_headers",
    "create_config_template",
    "make_SNPID_unique",
    "load_and_munge",
]


def read_config(config_path: str) -> Dict[str, Any]:
    """
    Read configuration from JSON file.

    Parameters
    ----------
    config_path : str
        Path to JSON configuration file.

    Returns
    -------
    Dict[str, Any]
        Configuration dictionary.

    Raises
    ------
    FileNotFoundError
        If configuration file doesn't exist.
    ValueError
        If configuration file is invalid JSON.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    try:
        with open(config_path, "r") as f:
            config = json.load(f)

        logger.info(f"Loaded configuration from {config_path}")
        return config

    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in configuration file: {e}")


def load_and_munge(
    file_path: str,
    config: Optional[Dict[str, Any]] = None,
    sep: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Load file and perform munging in one step.

    Parameters
    ----------
    file_path : str
        Path to input file.
    config : Dict[str, Any], optional
        Configuration dictionary with column mappings.
    sep : str, optional
        Column separator. Auto-detected if None.
    **kwargs
        Additional arguments for pandas.read_csv.

    Returns
    -------
    pd.DataFrame
        Munged DataFrame.
    """
    logger.info(f"Loading and munging file: {file_path}")

    # Auto-detect separator if not provided
    if sep is None:
        from .headers import _detect_separator

        sep = _detect_separator(file_path)

    # Load data
    compression = "gzip" if file_path.endswith(".gz") else None
    df = pd.read_csv(file_path, sep=sep, compression=compression, **kwargs)

    logger.info(f"Loaded {len(df)} rows with {len(df.columns)} columns")

    # Apply column mapping if provided in config
    if config and "column_mapping" in config:
        df = apply_header_mapping(df, config["column_mapping"])
    else:
        # Auto-map headers
        headers = df.columns.tolist()
        mapping = map_headers_automatic(headers)
        df = apply_header_mapping(df, mapping)

    # Validate required columns are present
    if not validate_required_columns(df):
        missing = set(ColName.mandatory_cols) - set(df.columns)
        raise ValueError(f"Cannot proceed: missing required columns {missing}")

    # Perform munging
    munged_df = munge(df)

    logger.info(f"Munging complete: {len(munged_df)} rows retained")

    return munged_df


def create_example_config(
    headers: List[str], output_path: str = "munge_config.json"
) -> str:
    """
    Create an example configuration file.

    Parameters
    ----------
    headers : List[str]
        List of column headers to create config for.
    output_path : str, optional
        Output path for configuration file.

    Returns
    -------
    str
        Path to created configuration file.
    """
    config = create_config_template(headers, interactive=False)

    with open(output_path, "w") as f:
        json.dump(config, f, indent=2)

    logger.info(f"Created example configuration: {output_path}")
    return output_path


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate configuration dictionary.

    Parameters
    ----------
    config : Dict[str, Any]
        Configuration to validate.

    Returns
    -------
    bool
        True if configuration is valid.
    """
    required_keys = ["column_mapping"]

    for key in required_keys:
        if key not in config:
            logger.error(f"Missing required configuration key: {key}")
            return False

    # Validate column mapping
    mapping = config["column_mapping"]
    if not isinstance(mapping, dict):
        logger.error("column_mapping must be a dictionary")
        return False

    # Check if we can map required columns
    mapped_standards = set(mapping.values())
    required_standards = set(ColName.mandatory_cols)
    missing = required_standards - mapped_standards

    if missing:
        logger.error(f"Configuration missing mappings for required columns: {missing}")
        return False

    logger.info("Configuration validation passed")
    return True
