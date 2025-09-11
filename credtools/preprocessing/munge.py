"""Munge summary statistics using smunger integration.

This module provides functionality to reformat and standardize GWAS summary
statistics from various formats into a consistent format suitable for fine-mapping.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
from tqdm import tqdm

logger = logging.getLogger("Preprocessing")


def munge_sumstats(
    input_files: Union[str, List[str], Dict[str, str]],
    output_dir: str,
    config_file: Optional[str] = None,
    force_overwrite: bool = False,
    **kwargs,
) -> Dict[str, str]:
    """
    Munge summary statistics files using smunger integration.

    Parameters
    ----------
    input_files : Union[str, List[str], Dict[str, str]]
        Input summary statistics file(s). Can be:
        - Single file path (str)
        - List of file paths (List[str])
        - Dictionary mapping ancestry/cohort names to file paths (Dict[str, str])
    output_dir : str
        Output directory for munged files.
    config_file : Optional[str], optional
        Path to configuration file specifying column mappings, by default None.
    force_overwrite : bool, optional
        Whether to overwrite existing output files, by default False.
    **kwargs
        Additional arguments passed to smunger functions.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping input identifiers to output file paths.

    Raises
    ------
    ImportError
        If smunger package is not available.
    FileNotFoundError
        If input files do not exist.
    ValueError
        If input format is invalid.

    Examples
    --------
    >>> # Single file
    >>> result = munge_sumstats("gwas_eur.txt", "output/")
    >>>
    >>> # Multiple files with ancestry labels
    >>> files = {"EUR": "gwas_eur.txt", "ASN": "gwas_asn.txt"}
    >>> result = munge_sumstats(files, "output/")
    """
    # Use internal munging module (adapted from smunger)
    from .munging import load_and_munge, munge, read_config

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Normalize input files to dictionary format
    if isinstance(input_files, str):
        # Single file - use filename as key
        file_key = Path(input_files).stem
        input_dict = {file_key: input_files}
    elif isinstance(input_files, list):
        # List of files - use filenames as keys
        input_dict = {Path(f).stem: f for f in input_files}
    elif isinstance(input_files, dict):
        # Already a dictionary
        input_dict = input_files
    else:
        raise ValueError("input_files must be a string, list of strings, or dictionary")

    # Validate input files exist
    for identifier, file_path in input_dict.items():
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Input file not found: {file_path}")

    # Load configuration if provided
    config = {}
    if config_file and os.path.exists(config_file):
        logger.info(f"Loading configuration from {config_file}")
        config = read_config(config_file)

    # Process each file
    output_files = {}
    for identifier, input_file in tqdm(input_dict.items(), desc="Munging files"):
        logger.info(f"Processing {identifier}: {input_file}")

        # Define output file path
        output_file = os.path.join(output_dir, f"{identifier}.munged.txt.gz")
        output_files[identifier] = output_file

        # Skip if output exists and not overwriting
        if os.path.exists(output_file) and not force_overwrite:
            logger.info(f"Output file exists, skipping: {output_file}")
            continue

        try:
            # Apply file-specific config if available
            file_config = config.get(identifier, {})

            # Munge the file using internal munging module
            logger.info(f"Munging {input_file} -> {output_file}")
            munged_data = load_and_munge(input_file, config=file_config, **kwargs)

            # Save munged data
            logger.info(f"Saving munged data to {output_file}")
            munged_data.to_csv(output_file, sep="\t", index=False, compression="gzip")

        except Exception as e:
            logger.error(f"Failed to process {input_file}: {str(e)}")
            # Remove failed output file if it was created
            if os.path.exists(output_file):
                os.remove(output_file)
            raise

    logger.info(f"Successfully munged {len(output_files)} files")
    return output_files


def create_munge_config(
    sample_files: Dict[str, str], output_config: str, interactive: bool = True
) -> None:
    """
    Create configuration file for munging by examining sample files.

    Parameters
    ----------
    sample_files : Dict[str, str]
        Dictionary mapping identifiers to sample file paths.
    output_config : str
        Output path for the configuration file.
    interactive : bool, optional
        Whether to use interactive mode for column mapping, by default True.

    Notes
    -----
    This function helps users create configuration files by examining
    the headers of input files and providing suggested column mappings.
    """
    # Use internal munging module (adapted from smunger)
    from .munging import create_config_template, inspect_headers

    config = {}

    for identifier, file_path in sample_files.items():
        logger.info(f"Examining headers for {identifier}: {file_path}")

        # Inspect file headers
        headers = inspect_headers(file_path)
        logger.info(f"Detected headers: {headers}")

        if interactive:
            # Interactive column mapping
            logger.info(f"Please specify column mappings for {identifier}")
            file_config = create_config_template(headers, interactive=True)
        else:
            # Auto-detect column mappings
            file_config = create_config_template(headers, interactive=False)

        config[identifier] = file_config

    # Save configuration
    import json

    with open(output_config, "w") as f:
        json.dump(config, f, indent=2)

    logger.info(f"Configuration saved to {output_config}")


def validate_munged_files(
    munged_files: Dict[str, str], required_columns: Optional[List[str]] = None
) -> Dict[str, Dict]:
    """
    Validate munged summary statistics files.

    Parameters
    ----------
    munged_files : Dict[str, str]
        Dictionary mapping identifiers to munged file paths.
    required_columns : Optional[List[str]], optional
        List of required column names to check, by default None.

    Returns
    -------
    Dict[str, Dict]
        Dictionary with validation results for each file.
    """
    if required_columns is None:
        required_columns = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]

    validation_results = {}

    for identifier, file_path in munged_files.items():
        logger.info(f"Validating {identifier}: {file_path}")

        result = {
            "file_exists": os.path.exists(file_path),
            "n_variants": 0,
            "columns": [],
            "missing_columns": [],
            "validation_passed": False,
        }

        if result["file_exists"]:
            try:
                # Read first few lines to check structure
                df_sample = pd.read_csv(
                    file_path, sep="\t", nrows=5, compression="gzip"
                )
                result["columns"] = df_sample.columns.tolist()
                result["missing_columns"] = [
                    col for col in required_columns if col not in result["columns"]
                ]

                # Count total variants
                df_full = pd.read_csv(file_path, sep="\t", compression="gzip")
                result["n_variants"] = len(df_full)

                # Check if validation passed
                result["validation_passed"] = len(result["missing_columns"]) == 0

                logger.info(
                    f"{identifier}: {result['n_variants']} variants, "
                    f"columns: {result['columns']}"
                )

                if result["missing_columns"]:
                    logger.warning(
                        f"{identifier}: Missing required columns: {result['missing_columns']}"
                    )

            except Exception as e:
                logger.error(f"Error validating {file_path}: {str(e)}")
                result["error"] = str(e)

        validation_results[identifier] = result

    return validation_results
