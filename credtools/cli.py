"""Console script for credtools."""

import json
import logging
import logging.handlers
import os
from enum import Enum
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeRemainingColumn,
)

from credtools import __version__
from credtools.credtools import fine_map, pipeline
from credtools.locus import check_loci_info, load_locus_set
from credtools.meta import meta_loci
from credtools.qc import loci_qc

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
app = typer.Typer(
    context_settings=CONTEXT_SETTINGS,
    add_completion=False,
    pretty_exceptions_enable=False,
)


class MetaMethod(str, Enum):
    """The method to perform meta-analysis."""

    meta_all = "meta_all"
    meta_by_population = "meta_by_population"
    no_meta = "no_meta"


class Tool(str, Enum):
    """The tool to perform fine-mapping."""

    abf = "abf"
    abf_cojo = "abf_cojo"
    finemap = "finemap"
    rsparsepro = "rsparsepro"
    susie = "susie"
    multisusie = "multisusie"
    susiex = "susiex"


class CombineCred(str, Enum):
    """Method to combine credible sets from multiple analyses."""

    union = "union"
    intersection = "intersection"
    cluster = "cluster"


class CombinePIP(str, Enum):
    """Method to combine posterior inclusion probabilities."""

    max = "max"
    min = "min"
    mean = "mean"
    meta = "meta"


def setup_file_logging(log_file: Optional[str], verbose: bool = False) -> None:
    """Set up file and console logging configuration."""
    if log_file is None:
        return

    # Create formatter
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Set up file handler
    try:
        file_handler = logging.FileHandler(log_file, mode="a")
        file_handler.setFormatter(formatter)

        # Set logging level
        if verbose:
            file_handler.setLevel(logging.DEBUG)
        else:
            file_handler.setLevel(logging.INFO)

        # Add file handler to root logger and specific loggers
        root_logger = logging.getLogger()
        root_logger.addHandler(file_handler)

        # Add to specific credtools loggers
        for name in [
            "CREDTOOLS",
            "FINEMAP",
            "RSparsePro",
            "COJO",
            "SuSiE",
            "MULTISUSIE",
            "SUSIEX",
            "ABF",
            "ABF_COJO",
            "Locus",
            "LDMatrix",
            "QC",
            "Sumstats",
            "Utils",
        ]:
            logger = logging.getLogger(name)
            logger.addHandler(file_handler)

        logging.info(f"Logging to file: {log_file}")

    except (OSError, IOError) as e:
        console = Console()
        console.print(f"[red]Warning: Could not set up log file {log_file}: {e}[/red]")


@app.callback(invoke_without_command=True, no_args_is_help=True)
def main(
    version: bool = typer.Option(False, "--version", "-V", help="Show version."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show verbose info."),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """CREDTOOLS: Credible Set Tools for fine-mapping analysis."""
    console = Console()
    console.rule("[bold blue]CREDTOOLS[/bold blue]")
    console.print(f"Version: {__version__}", justify="center")
    console.print("Author: Jianhua Wang", justify="center")
    console.print("Email: jianhua.mert@gmail.com", justify="center")
    if version:
        typer.echo(f"CREDTOOLS version: {__version__}")
        raise typer.Exit()
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.info("Verbose mode is on.")
    else:
        for name in [
            "CREDTOOLS",
            "FINEMAP",
            "RSparsePro",
            "COJO",
            "SuSiE",
            "MULTISUSIE",
            "SUSIEX",
            "ABF",
            "ABF_COJO",
            "Locus",
            "LDMatrix",
            "QC",
            "Sumstats",
            "Utils",
        ]:
            logging.getLogger(name).setLevel(logging.INFO)
        # logging.getLogger().setLevel(logging.INFO)

    # Set up file logging if requested
    setup_file_logging(log_file, verbose)


def parse_population_config_file(
    config_file_path: str,
) -> tuple[Dict[str, str], Dict[str, str], pd.DataFrame]:
    """
    Parse population configuration file with columns: popu, cohort, sample_size, path, ld_ref.

    Parameters
    ----------
    config_file_path : str
        Path to the configuration file.

    Returns
    -------
    tuple[Dict[str, str], Dict[str, str], pd.DataFrame]
        Tuple containing:
        - Dictionary mapping population identifiers to sumstats file paths
        - Dictionary mapping population identifiers to LD reference file paths
        - Original DataFrame for later updating
    """
    if not os.path.exists(config_file_path):
        raise FileNotFoundError(f"Configuration file not found: {config_file_path}")

    try:
        # Read the configuration file
        config_df = pd.read_csv(config_file_path, sep="\t")

        # Check required columns
        required_cols = ["popu", "cohort", "sample_size", "path", "ld_ref"]
        missing_cols = [col for col in required_cols if col not in config_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in config file: {missing_cols}")

        # Create identifier -> path mappings
        sumstats_dict = {}
        ld_ref_dict = {}
        for _, row in config_df.iterrows():
            # Create identifier from population and cohort
            identifier = f"{row['popu']}_{row['cohort']}"
            sumstats_dict[identifier] = row["path"]
            ld_ref_dict[identifier] = row["ld_ref"]

            # Check if files exist
            if not os.path.exists(row["path"]):
                raise FileNotFoundError(
                    f"Summary statistics file not found: {row['path']}"
                )

            # For LD reference, check for common PLINK file extensions
            ld_base = row["ld_ref"]
            if not any(
                os.path.exists(f"{ld_base}.{ext}") for ext in ["bed", "bim", "fam"]
            ):
                raise FileNotFoundError(
                    f"LD reference files not found: {ld_base}.[bed/bim/fam]"
                )

        return sumstats_dict, ld_ref_dict, config_df

    except Exception as e:
        raise ValueError(f"Error parsing configuration file: {e}")


def parse_population_config_file_munge_only(
    config_file_path: str,
) -> tuple[Dict[str, str], pd.DataFrame]:
    """
    Parse population configuration file for munge command (backward compatibility).

    Only requires: popu, cohort, sample_size, path.

    Parameters
    ----------
    config_file_path : str
        Path to the configuration file.

    Returns
    -------
    tuple[Dict[str, str], pd.DataFrame]
        Tuple containing:
        - Dictionary mapping population identifiers to file paths
        - Original DataFrame for later updating
    """
    if not os.path.exists(config_file_path):
        raise FileNotFoundError(f"Configuration file not found: {config_file_path}")

    try:
        # Read the configuration file
        config_df = pd.read_csv(config_file_path, sep="\t")

        # Check required columns
        required_cols = ["popu", "cohort", "sample_size", "path"]
        missing_cols = [col for col in required_cols if col not in config_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in config file: {missing_cols}")

        # Create identifier -> path mapping
        input_dict = {}
        for _, row in config_df.iterrows():
            # Create identifier from population and cohort
            identifier = f"{row['popu']}_{row['cohort']}"
            input_dict[identifier] = row["path"]

            # Check if file exists
            if not os.path.exists(row["path"]):
                raise FileNotFoundError(
                    f"Summary statistics file not found: {row['path']}"
                )

        return input_dict, config_df

    except Exception as e:
        raise ValueError(f"Error parsing configuration file: {e}")


def create_updated_sumstat_info(
    original_config_df: pd.DataFrame, munged_files: Dict[str, str], output_path: str
) -> str:
    """
    Create an updated sumstat info file with new paths pointing to munged files.

    Parameters
    ----------
    original_config_df : pd.DataFrame
        Original configuration DataFrame
    munged_files : Dict[str, str]
        Dictionary mapping identifiers to munged file paths
    output_path : str
        Path for the output file

    Returns
    -------
    str
        Path to the created file
    """
    # Create a copy of the original config
    updated_config = original_config_df.copy()

    # Update the path column with munged file paths
    for idx, row in updated_config.iterrows():
        identifier = f"{row['popu']}_{row['cohort']}"
        if identifier in munged_files:
            updated_config.at[idx, "path"] = munged_files[identifier]

    # Save the updated configuration
    updated_config.to_csv(output_path, sep="\t", index=False)

    return output_path


def create_updated_chunk_info(
    original_config_df: pd.DataFrame, chunk_info_df: pd.DataFrame, output_path: str
) -> str:
    """
    Create an updated sumstat info file with paths pointing to chunked files.

    Parameters
    ----------
    original_config_df : pd.DataFrame
        Original configuration DataFrame
    chunk_info_df : pd.DataFrame
        DataFrame from chunk_sumstats with chunked file information
    output_path : str
        Path for the output file

    Returns
    -------
    str
        Path to the created file
    """
    # Group chunk_info_df by ancestry to get the base directory for each ancestry
    chunk_files_by_ancestry = {}

    for _, row in chunk_info_df.iterrows():
        ancestry = row["ancestry"]
        if ancestry not in chunk_files_by_ancestry:
            # Get the directory containing chunked files for this ancestry
            chunk_dir = os.path.dirname(row["sumstats_file"])
            chunk_files_by_ancestry[ancestry] = chunk_dir

    # Create a copy of the original config
    updated_config = original_config_df.copy()

    # Update the path column with chunk directory paths
    for idx, row in updated_config.iterrows():
        identifier = f"{row['popu']}_{row['cohort']}"

        # Try to match the identifier with ancestry in chunk_files
        if identifier in chunk_files_by_ancestry:
            updated_config.at[idx, "path"] = chunk_files_by_ancestry[identifier]

    # Save the updated configuration
    updated_config.to_csv(output_path, sep="\t", index=False)

    return output_path


@app.command(
    name="munge",
    help="Reformat and standardize GWAS summary statistics.",
)
def run_munge(
    input_config: str = typer.Argument(
        ...,
        help="Input configuration file with columns: popu, cohort, sample_size, path.",
    ),
    output_dir: str = typer.Argument(..., help="Output directory for munged files."),
    config_file: Optional[str] = typer.Option(
        None, "--config", "-c", help="Configuration file for column mappings."
    ),
    force_overwrite: bool = typer.Option(
        False, "--force", "-f", help="Overwrite existing output files."
    ),
    interactive_config: bool = typer.Option(
        False, "--interactive", "-i", help="Create configuration interactively."
    ),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Reformat and standardize GWAS summary statistics from population configuration file."""
    setup_file_logging(log_file)

    try:
        from credtools.preprocessing import munge_sumstats
        from credtools.preprocessing.munge import (
            create_munge_config,
            validate_munged_files,
        )
    except ImportError as e:
        console = Console()
        console.print("[red]Error: Preprocessing dependencies not found.[/red]")
        console.print("Please ensure smunger is installed: pip install smunger")
        raise typer.Exit(1) from e

    console = Console()
    console.print("[cyan]Munging summary statistics...[/cyan]")

    original_config_df: Optional[pd.DataFrame] = None

    def _parse_direct_inputs(input_spec: str) -> Dict[str, str]:
        normalized = input_spec.replace("\n", ",")
        paths = [part.strip() for part in normalized.split(",") if part.strip()]
        if not paths:
            raise ValueError("No input files provided.")
        input_mapping: Dict[str, str] = {}
        for idx, candidate in enumerate(paths, start=1):
            expanded = os.path.expanduser(candidate)
            if not os.path.exists(expanded):
                raise FileNotFoundError(f"Input file not found: {expanded}")
            key = Path(expanded).stem
            if key in input_mapping:
                key = f"{key}_{idx}"
            input_mapping[key] = expanded
        return input_mapping

    input_dict: Dict[str, str]
    config_candidate = os.path.expanduser(input_config)
    config_loaded = False

    if (
        "," not in input_config
        and "\n" not in input_config
        and os.path.exists(config_candidate)
    ):
        try:
            input_dict, original_config_df = parse_population_config_file_munge_only(
                config_candidate
            )
        except Exception:
            config_loaded = False
        else:
            config_loaded = True
            console.print(
                f"[green]Loaded {len(input_dict)} population files from config[/green]"
            )

    if not config_loaded:
        try:
            input_dict = _parse_direct_inputs(input_config)
        except (FileNotFoundError, ValueError) as direct_error:
            console.print(f"[red]Error parsing input files: {direct_error}[/red]")
            raise typer.Exit(1) from direct_error
        else:
            console.print(
                f"[green]Loaded {len(input_dict)} input file(s) from direct paths[/green]"
            )

    # Create interactive config if requested
    if interactive_config:
        config_output = config_file or os.path.join(output_dir, "munge_config.json")
        console.print(f"[yellow]Creating configuration file: {config_output}[/yellow]")
        create_munge_config(input_dict, config_output, interactive=True)
        config_file = config_output

    # Perform munging
    try:
        result = munge_sumstats(
            input_files=input_dict,
            output_dir=output_dir,
            config_file=config_file,
            force_overwrite=force_overwrite,
        )

        # Validate results with updated required columns
        required_columns = [
            "CHR",
            "BP",
            "SNPID",
            "EA",
            "NEA",
            "EAF",
            "BETA",
            "SE",
            "P",
            "N",
            "RSID",
        ]
        validation = validate_munged_files(result, required_columns=required_columns)

        # Create updated sumstat info file when a configuration was provided
        if original_config_df is not None:
            updated_info_path = os.path.join(output_dir, "sumstat_info_updated.txt")
            try:
                created_info_file = create_updated_sumstat_info(
                    original_config_df, result, updated_info_path
                )
                console.print(
                    f"[green]Created updated sumstat info file: {created_info_file}[/green]"
                )
            except Exception as exc:
                console.print(
                    f"[yellow]Warning: Could not create updated info file: {exc}[/yellow]"
                )

        # Print summary
        console.print(f"[green]Successfully munged {len(result)} files[/green]")
        for identifier, file_path in result.items():
            val = validation[identifier]
            status = "✓" if val["validation_passed"] else "✗"
            console.print(
                f"  {status} {identifier}: {val['n_variants']} variants -> {file_path}"
            )

    except Exception as e:
        console.print(f"[red]Error during munging: {e}[/red]")
        raise typer.Exit(1)


def _load_custom_chunks(custom_chunks_file: str) -> pd.DataFrame:
    """
    Load custom chunk definitions from file.

    Parameters
    ----------
    custom_chunks_file : str
        Path to custom chunks file with chr, start, end columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with loci coordinates.
    """
    if not os.path.exists(custom_chunks_file):
        raise FileNotFoundError(f"Custom chunks file not found: {custom_chunks_file}")

    try:
        chunks_df = pd.read_csv(custom_chunks_file, sep="\t")

        # Check required columns
        required_cols = ["chr", "start", "end"]
        missing_cols = [col for col in required_cols if col not in chunks_df.columns]
        if missing_cols:
            raise ValueError(
                f"Missing required columns in custom chunks file: {missing_cols}"
            )

        # Create locus_id
        chunks_df["locus_id"] = [
            f"chr{row['chr']}_{row['start']}_{row['end']}"
            for _, row in chunks_df.iterrows()
        ]

        # Add placeholder columns to match identify_independent_loci output
        chunks_df["lead_snp"] = None
        chunks_df["lead_bp"] = (chunks_df["start"] + chunks_df["end"]) // 2
        chunks_df["lead_p"] = None
        chunks_df["ancestry"] = "custom"
        chunks_df["n_variants"] = 0

        return chunks_df[
            [
                "chr",
                "start",
                "end",
                "locus_id",
                "lead_snp",
                "lead_bp",
                "lead_p",
                "ancestry",
                "n_variants",
            ]
        ]

    except Exception as e:
        raise ValueError(f"Error parsing custom chunks file: {e}")


def _prepare_ld_matrices(
    chunk_info_df: pd.DataFrame,
    ld_ref_dict: Dict[str, str],
    output_dir: str,
    threads: int = 1,
    ld_format: str = "plink",
    keep_intermediate: bool = False,
) -> pd.DataFrame:
    """
    Prepare LD matrices for chunked files.

    Parameters
    ----------
    chunk_info_df : pd.DataFrame
        DataFrame from chunk_sumstats with chunked file information.
    ld_ref_dict : Dict[str, str]
        Dictionary mapping ancestry names to LD reference file prefixes.
    output_dir : str
        Output directory for prepared files.
    threads : int, optional
        Number of threads, by default 1.
    ld_format : str, optional
        LD format, by default "plink".
    keep_intermediate : bool, optional
        Keep intermediate files, by default False.

    Returns
    -------
    pd.DataFrame
        DataFrame with prepared file information.
    """
    try:
        from credtools.preprocessing.prepare import prepare_finemap_inputs
    except ImportError as e:
        raise ImportError("Could not import prepare function") from e

    os.makedirs(output_dir, exist_ok=True)

    # Convert chunk_info_df to the format expected by prepare_finemap_inputs
    # We need to map ancestry names to genotype file prefixes
    genotype_files = {}
    for ancestry in chunk_info_df["ancestry"].unique():
        # Find matching LD reference for this ancestry
        matching_key = None
        for key in ld_ref_dict.keys():
            if ancestry in key:
                matching_key = key
                break

        if matching_key:
            genotype_files[ancestry] = ld_ref_dict[matching_key]
        else:
            raise ValueError(f"No LD reference found for ancestry: {ancestry}")

    # Rename columns to match prepare function expectations
    prep_chunk_df = chunk_info_df.copy()
    prep_chunk_df = prep_chunk_df.rename(columns={"ancestry": "popu"})

    # Add required columns
    prep_chunk_df["cohort"] = prep_chunk_df["popu"]
    prep_chunk_df["sample_size"] = 50000  # Placeholder

    # Add prefix column based on sumstats_file
    prep_chunk_df["prefix"] = prep_chunk_df["sumstats_file"].apply(
        lambda x: str(Path(x).with_suffix("")).replace(".sumstats", "")
    )

    # Call the prepare function
    prepared_df = prepare_finemap_inputs(
        chunk_info_df=prep_chunk_df,
        genotype_files=genotype_files,
        output_dir=output_dir,
        threads=threads,
        ld_format=ld_format,
        keep_intermediate=keep_intermediate,
    )

    return prepared_df


def _update_chunk_info_with_prepared(
    chunk_info_df: pd.DataFrame, prepared_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Update chunk info DataFrame with prepared file information.

    Parameters
    ----------
    chunk_info_df : pd.DataFrame
        Original chunk info DataFrame with sumstats_file column.
    prepared_df : pd.DataFrame
        Prepared files DataFrame with prefix column.

    Returns
    -------
    pd.DataFrame
        Updated chunk info DataFrame with prepared file prefixes.
    """
    updated_df = chunk_info_df.copy()

    # Create a mapping from locus_id + ancestry to prepared prefix
    prepared_mapping = {}
    for _, row in prepared_df.iterrows():
        key = (row["locus_id"], row["popu"])
        prepared_mapping[key] = row["prefix"]

    # Update sumstats_file with prepared prefix (add extensions for credtools compatibility)
    for idx, row in updated_df.iterrows():
        key = (row["locus_id"], row["ancestry"])
        if key in prepared_mapping:
            # Update to point to prepared files
            updated_df.at[idx, "sumstats_file"] = prepared_mapping[key] + ".sumstats.gz"

    return updated_df


@app.command(
    name="chunk",
    help="Identify independent loci, chunk summary statistics, and extract LD matrices.",
)
def run_chunk(
    input_config: str = typer.Argument(
        ...,
        help="Input configuration file with columns: popu, cohort, sample_size, path, ld_ref.",
    ),
    output_dir: str = typer.Argument(..., help="Output directory for chunked files."),
    distance_threshold: int = typer.Option(
        500000, "--distance", "-d", help="Distance threshold for independence (bp)."
    ),
    pvalue_threshold: float = typer.Option(
        5e-8, "--pvalue", "-p", help="P-value threshold for significance."
    ),
    merge_overlapping: bool = typer.Option(
        True,
        "--merge-overlapping",
        "-m",
        help="Merge overlapping loci across ancestries.",
    ),
    use_most_sig_if_no_sig: bool = typer.Option(
        True,
        "--use-most-sig",
        "-u",
        help="Use most significant SNP if no significant SNPs.",
    ),
    min_variants_per_locus: int = typer.Option(
        10, "--min-variants", "-v", help="Minimum variants per locus."
    ),
    custom_chunks: Optional[str] = typer.Option(
        None,
        "--custom-chunks",
        "-cc",
        help="Custom chunk file with chr, start, end columns.",
    ),
    ld_format: str = typer.Option(
        "plink", "--ld-format", "-f", help="LD computation format (plink/vcf)."
    ),
    keep_intermediate: bool = typer.Option(
        False, "--keep-intermediate", "-k", help="Keep intermediate files."
    ),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Identify independent loci, chunk summary statistics, and extract LD matrices from GWAS info configuration file."""
    setup_file_logging(log_file)

    try:
        from credtools.preprocessing import chunk_sumstats, identify_independent_loci
        from credtools.preprocessing.chunk import create_loci_list_for_credtools
    except ImportError as e:
        console = Console()
        console.print("[red]Error: Preprocessing module not found.[/red]")
        raise typer.Exit(1) from e

    console = Console()

    original_config_df: Optional[pd.DataFrame] = None
    sumstats_dict: Dict[str, str]
    ld_ref_dict: Dict[str, str] = {}

    def _parse_direct_chunk_inputs(input_spec: str) -> Dict[str, str]:
        normalized = input_spec.replace("\n", ",")
        paths = [part.strip() for part in normalized.split(",") if part.strip()]
        if not paths:
            raise ValueError("No input files provided.")
        input_mapping: Dict[str, str] = {}
        for idx, candidate in enumerate(paths, start=1):
            expanded = os.path.expanduser(candidate)
            if not os.path.exists(expanded):
                raise FileNotFoundError(f"Input file not found: {expanded}")
            key_name = Path(expanded).name
            if key_name.endswith(".txt.gz"):
                key_name = key_name[: -len(".txt.gz")]
            elif key_name.endswith(".gz"):
                key_name = key_name[: -len(".gz")]
            if key_name.endswith(".sumstats"):
                key_name = key_name[: -len(".sumstats")]
            if ".munged" in key_name:
                key_name = key_name.split(".munged")[0]
            if key_name in input_mapping:
                key_name = f"{key_name}_{idx}"
            input_mapping[key_name] = expanded
        return input_mapping

    config_candidate = os.path.expanduser(input_config)
    config_loaded = False
    if (
        "," not in input_config
        and "\n" not in input_config
        and os.path.exists(config_candidate)
    ):
        try:
            sumstats_dict, ld_ref_dict, original_config_df = (
                parse_population_config_file(config_candidate)
            )
        except Exception:
            config_loaded = False
        else:
            config_loaded = True
            console.print(
                f"[green]Loaded {len(sumstats_dict)} population files from config[/green]"
            )

    if not config_loaded:
        try:
            sumstats_dict = _parse_direct_chunk_inputs(input_config)
        except (FileNotFoundError, ValueError) as direct_error:
            console.print(f"[red]Error parsing input files: {direct_error}[/red]")
            raise typer.Exit(1) from direct_error
        else:
            console.print(
                f"[green]Loaded {len(sumstats_dict)} input file(s) from direct paths[/green]"
            )

    try:
        # Load or identify loci
        if custom_chunks:
            console.print(f"[cyan]Loading custom chunks from {custom_chunks}...[/cyan]")
            loci_df = _load_custom_chunks(custom_chunks)
        else:
            console.print("[cyan]Identifying independent loci...[/cyan]")
            loci_df = identify_independent_loci(
                sumstats_files=sumstats_dict,
                output_dir=output_dir,
                distance_threshold=distance_threshold,
                pvalue_threshold=pvalue_threshold,
                merge_overlapping=merge_overlapping,
                use_most_sig_if_no_sig=use_most_sig_if_no_sig,
                min_variants_per_locus=min_variants_per_locus,
            )

        if len(loci_df) == 0:
            console.print("[yellow]No loci identified[/yellow]")
            return

        # Chunk summary statistics and optionally extract LD matrices
        console.print(f"[cyan]Chunking {len(loci_df)} loci...[/cyan]")
        chunk_info_df = chunk_sumstats(
            loci_df=loci_df,
            sumstats_files=sumstats_dict,
            output_dir=os.path.join(output_dir, "chunks"),
            threads=threads,
        )

        if ld_ref_dict:
            console.print("[cyan]Extracting LD matrices...[/cyan]")
            try:
                prepared_df = _prepare_ld_matrices(
                    chunk_info_df=chunk_info_df,
                    ld_ref_dict=ld_ref_dict,
                    output_dir=os.path.join(output_dir, "prepared"),
                    threads=threads,
                    ld_format=ld_format,
                    keep_intermediate=keep_intermediate,
                )
            except Exception as e:
                console.print(
                    f"[yellow]Warning: LD extraction had issues: {e}[/yellow]"
                )
                console.print("[cyan]Continuing with chunk files only...[/cyan]")
                prepared_df = chunk_info_df
        else:
            console.print(
                "[yellow]No LD reference information provided; skipping LD matrix extraction[/yellow]"
            )
            prepared_df = chunk_info_df

        # Create credtools-compatible loci list from prepared files
        loci_list_file = os.path.join(output_dir, "loci_list.txt")
        # Update chunk_info_df with prepared file prefixes (if LD extraction succeeded)
        if "prefix" in prepared_df.columns:
            updated_chunk_df = _update_chunk_info_with_prepared(
                chunk_info_df, prepared_df
            )
        else:
            updated_chunk_df = chunk_info_df  # Use original chunk files
        credtools_df = create_loci_list_for_credtools(
            chunk_info_df=updated_chunk_df, output_file=loci_list_file
        )

        # Create updated sumstat info file when configuration is available
        if original_config_df is not None:
            updated_info_path = os.path.join(output_dir, "sumstat_info_updated.txt")
            try:
                created_info_file = create_updated_chunk_info(
                    original_config_df, chunk_info_df, updated_info_path
                )
                console.print(
                    f"[green]Created updated sumstat info file: {created_info_file}[/green]"
                )
            except Exception as exc:
                console.print(
                    f"[yellow]Warning: Could not create updated info file: {exc}[/yellow]"
                )
        else:
            console.print(
                "[yellow]Skipping creation of updated sumstat info file (no configuration provided)[/yellow]"
            )

        # Print summary
        console.print(f"[green]Successfully processed {len(loci_df)} loci[/green]")
        console.print(f"[green]Generated {len(chunk_info_df)} chunked files[/green]")
        if "prefix" in prepared_df.columns:
            console.print(
                f"[green]Generated {len(prepared_df)} prepared files with LD matrices[/green]"
            )
        else:
            console.print(
                "[yellow]LD matrix extraction failed, using chunked files only[/yellow]"
            )
        console.print(f"[green]Credtools loci list: {loci_list_file}[/green]")

    except Exception as e:
        console.print(f"[red]Error during chunking: {e}[/red]")
        raise typer.Exit(1)


@app.command(
    name="meta",
    help="Meta-analysis of summary statistics and LD matrices.",
)
def run_meta(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
    meta_method: MetaMethod = typer.Option(
        MetaMethod.meta_all, "--meta-method", "-m", help="Meta-analysis method."
    ),
    calculate_lambda_s: bool = typer.Option(
        False,
        "--calculate-lambda-s",
        "-cls",
        help="Calculate lambda_s parameter using estimate_s_rss function.",
    ),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Meta-analysis of summary statistics and LD matrices."""
    setup_file_logging(log_file)
    meta_loci(inputs, outdir, threads, meta_method, calculate_lambda_s)


@app.command(
    name="qc",
    help="Quality control of summary statistics and LD matrices.",
)
def run_qc(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Quality control of summary statistics and LD matrices."""
    setup_file_logging(log_file)
    run_summary = loci_qc(inputs, outdir, threads)

    console = Console()
    if run_summary["failed_loci"] > 0:
        console.print(
            f"[yellow]QC completed with {run_summary['failed_loci']} failed loci[/yellow]"
        )
    else:
        console.print(
            f"[green]QC completed successfully for all {run_summary['successful_loci']} loci[/green]"
        )
    console.print(f"[dim]QC run summary saved to {run_summary['log_path']}[/dim]")


@app.command(
    name="finemap",
    help="Perform fine-mapping analysis on genetic loci.",
)
def run_fine_map(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    tool: Tool = typer.Option(
        Tool.susie,
        "--tool",
        "-t",
        help="Fine-mapping tool. Single-input tools (abf, susie, etc.) process each locus individually. "
        "Multi-input tools (susiex, multisusie) process all loci together. "
        "When using single-input tools with multiple loci, results are automatically combined.",
    ),
    max_causal: int = typer.Option(
        5, "--max-causal", "-c", help="Maximum number of causal SNPs."
    ),
    adaptive_max_causal: bool = typer.Option(
        False,
        "--adaptive-max-causal",
        "-amc",
        help="Enable adaptive max_causal parameter tuning.",
    ),
    set_L_by_cojo: bool = typer.Option(
        False, "--set-L-by-cojo", "-sl", help="Set L by COJO."
    ),
    p_cutoff: float = typer.Option(
        5e-8, "--p-cutoff", "-pc", help="P-value cutoff for COJO."
    ),
    collinear_cutoff: float = typer.Option(
        0.9, "--collinear-cutoff", "-cc", help="Collinearity cutoff for COJO."
    ),
    window_size: int = typer.Option(
        10000000, "--window-size", "-ws", help="Window size for COJO."
    ),
    maf_cutoff: float = typer.Option(
        0.01, "--maf-cutoff", "-mc", help="MAF cutoff for COJO."
    ),
    diff_freq_cutoff: float = typer.Option(
        0.2,
        "--diff-freq-cutoff",
        "-dfc",
        help="Difference in frequency cutoff for COJO.",
    ),
    coverage: float = typer.Option(
        0.95, "--coverage", "-cv", help="Coverage of the credible set."
    ),
    timeout_minutes: float = typer.Option(
        30.0,
        "--timeout-minutes",
        "-tm",
        help="Maximum runtime per locus in minutes when running FINEMAP.",
    ),
    combine_cred: CombineCred = typer.Option(
        CombineCred.union,
        "--combine-cred",
        "-cc",
        help="Method to combine credible sets when using single-input tools with multiple loci.",
    ),
    combine_pip: CombinePIP = typer.Option(
        CombinePIP.max,
        "--combine-pip",
        "-cp",
        help="Method to combine PIPs when using single-input tools with multiple loci.",
    ),
    jaccard_threshold: float = typer.Option(
        0.1,
        "--jaccard-threshold",
        "-j",
        help="Jaccard threshold for combining credible sets.",
    ),
    # susie parameters
    max_iter: int = typer.Option(
        100, "--max-iter", "-i", help="Maximum number of iterations."
    ),
    estimate_residual_variance: bool = typer.Option(
        False, "--estimate-residual-variance", "-er", help="Estimate residual variance."
    ),
    min_abs_corr: float = typer.Option(
        0.5, "--min-abs-corr", "-mc", help="Minimum absolute correlation."
    ),
    convergence_tol: float = typer.Option(
        1e-3, "--convergence-tol", "-ct", help="Convergence tolerance."
    ),
    calculate_lambda_s: bool = typer.Option(
        False,
        "--calculate-lambda-s",
        "-cls",
        help="Calculate lambda_s parameter using estimate_s_rss function.",
    ),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Perform fine-mapping analysis on genetic loci.

    The appropriate analysis strategy is automatically determined based on:
    - Tool type: Single-input tools (abf, susie, finemap, etc.) vs Multi-input tools (susiex, multisusie)
    - Data structure: Single locus vs multiple loci

    When using single-input tools with multiple loci, results are automatically combined.
    """
    setup_file_logging(log_file)
    import logging
    import sys
    from datetime import datetime

    from credtools.utils import create_float_format_dict

    logger = logging.getLogger("CREDTOOLS")

    loci_info = pd.read_csv(inputs, sep="\t")
    loci_info = check_loci_info(loci_info)  # Validate input data

    # Initialize run summary
    run_summary = {
        "start_time": datetime.now().isoformat(),
        "total_loci": 0,
        "successful_loci": 0,
        "failed_loci": 0,
        "errors": [],
        "tool": tool.value,
        "parameters": {
            "max_causal": max_causal,
            "adaptive_max_causal": adaptive_max_causal,
            "set_L_by_cojo": set_L_by_cojo,
            "p_cutoff": p_cutoff,
            "coverage": coverage,
            "timeout_minutes": timeout_minutes,
            "combine_cred": combine_cred.value,
            "combine_pip": combine_pip.value,
        },
    }

    # Collect all credible sets for summary
    all_credible_sets = []

    # Create progress bar
    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
    )

    # Get total number of loci
    locus_groups = list(loci_info.groupby("locus_id"))
    total_loci = len(locus_groups)
    run_summary["total_loci"] = total_loci

    with progress:
        task = progress.add_task("[cyan]Fine-mapping loci...", total=total_loci)

        for locus_id, locus_info in locus_groups:
            try:
                locus_set = load_locus_set(
                    locus_info, calculate_lambda_s=calculate_lambda_s
                )
                creds = fine_map(
                    locus_set,
                    tool=tool,
                    max_causal=max_causal,
                    adaptive_max_causal=adaptive_max_causal,
                    set_L_by_cojo=set_L_by_cojo,
                    p_cutoff=p_cutoff,
                    collinear_cutoff=collinear_cutoff,
                    window_size=window_size,
                    maf_cutoff=maf_cutoff,
                    diff_freq_cutoff=diff_freq_cutoff,
                    coverage=coverage,
                    timeout_minutes=timeout_minutes,
                    combine_cred=combine_cred,
                    combine_pip=combine_pip,
                    jaccard_threshold=jaccard_threshold,
                    # susie parameters
                    max_iter=max_iter,
                    estimate_residual_variance=estimate_residual_variance,
                    min_abs_corr=min_abs_corr,
                    convergence_tol=convergence_tol,
                )

                # Create enhanced PIPs DataFrame
                enhanced_pips = creds.create_enhanced_pips_df(locus_set)

                # Get appropriate float formats for columns
                format_dict = create_float_format_dict(enhanced_pips)

                # Format numeric columns
                for col, fmt in format_dict.items():
                    if col in enhanced_pips.columns:
                        if fmt == "%.3e":
                            enhanced_pips[col] = enhanced_pips[col].apply(
                                lambda x: f"{x:.3e}" if pd.notna(x) else ""
                            )
                        elif fmt == "%.4f":
                            enhanced_pips[col] = enhanced_pips[col].apply(
                                lambda x: f"{x:.4f}" if pd.notna(x) else ""
                            )

                # Save enhanced PIPs with compression
                out_dir = f"{outdir}/{locus_id}"
                os.makedirs(out_dir, exist_ok=True)
                output_file = f"{out_dir}/pips.txt.gz"
                enhanced_pips.to_csv(
                    output_file,
                    sep="\t",
                    index=False,
                    compression="gzip",
                )

                # Collect credible sets for summary
                cs_summary = enhanced_pips[enhanced_pips["CRED"] != 0].copy()
                if len(cs_summary) > 0:
                    cs_summary["locus_id"] = locus_id
                    all_credible_sets.append(cs_summary)

                run_summary["successful_loci"] += 1

            except Exception as e:
                error_msg = f"Locus {locus_id} failed: {str(e)}"
                logger.error(error_msg)
                print(f"ERROR: {error_msg}", file=sys.stderr)
                run_summary["failed_loci"] += 1
                run_summary["errors"].append(error_msg)
                # Continue to next locus instead of failing

            progress.advance(task)

    # Save combined credible sets summary
    if all_credible_sets:
        combined_cs = pd.concat(all_credible_sets, ignore_index=True)
        combined_cs.to_csv(
            f"{outdir}/credible_sets_summary.txt.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    # Save parameters (simplified version)
    if run_summary["successful_loci"] > 0:
        parameters_dict = {
            "tool": tool.value,
            "n_loci_processed": run_summary["successful_loci"],
            "coverage": coverage,
            "parameters": run_summary["parameters"],
        }
        with open(f"{outdir}/parameters.json", "w") as f:
            json.dump(parameters_dict, f, indent=4)

    # Generate run summary
    run_summary["end_time"] = datetime.now().isoformat()
    with open(f"{outdir}/run_summary.log", "w") as f:
        f.write("=== CREDTOOLS FINE-MAPPING RUN SUMMARY ===\n")
        f.write(f"Start Time: {run_summary['start_time']}\n")
        f.write(f"End Time: {run_summary['end_time']}\n")
        f.write(f"Total Loci: {run_summary['total_loci']}\n")
        f.write(f"Successful: {run_summary['successful_loci']}\n")
        f.write(f"Failed: {run_summary['failed_loci']}\n")
        f.write("\n")

        if run_summary["errors"]:
            f.write("Error Details:\n")
            for error in run_summary["errors"]:
                f.write(f"  - {error}\n")
            f.write("\n")

        f.write("Parameters Used:\n")
        for key, value in run_summary["parameters"].items():
            f.write(f"  {key}: {value}\n")

    # Print summary to console
    console = Console()
    if run_summary["failed_loci"] > 0:
        console.print(
            f"[yellow]Completed with {run_summary['failed_loci']} failed loci[/yellow]"
        )
    else:
        console.print(
            f"[green]Successfully processed all {run_summary['successful_loci']} loci[/green]"
        )


@app.command(
    name="pipeline",
    help="Run whole fine-mapping pipeline on a list of loci.",
)
def run_pipeline(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    meta_method: MetaMethod = typer.Option(
        MetaMethod.meta_all, "--meta-method", "-m", help="Meta-analysis method."
    ),
    skip_qc: bool = typer.Option(
        False, "--skip-qc", "-q", help="Skip quality control."
    ),
    tool: Tool = typer.Option(
        Tool.susie,
        "--tool",
        "-t",
        help="Fine-mapping tool. Single-input tools process each locus individually, "
        "multi-input tools process all loci together.",
    ),
    max_causal: int = typer.Option(
        5, "--max-causal", "-c", help="Maximum number of causal SNPs."
    ),
    adaptive_max_causal: bool = typer.Option(
        False,
        "--adaptive-max-causal",
        "-amc",
        help="Enable adaptive max_causal parameter tuning.",
    ),
    set_L_by_cojo: bool = typer.Option(
        True, "--set-L-by-cojo", "-sl", help="Set L by COJO."
    ),
    coverage: float = typer.Option(
        0.95, "--coverage", "-cv", help="Coverage of the credible set."
    ),
    combine_cred: CombineCred = typer.Option(
        CombineCred.union,
        "--combine-cred",
        "-cc",
        help="Method to combine credible sets when using single-input tools with multiple loci.",
    ),
    combine_pip: CombinePIP = typer.Option(
        CombinePIP.max,
        "--combine-pip",
        "-cp",
        help="Method to combine PIPs when using single-input tools with multiple loci.",
    ),
    jaccard_threshold: float = typer.Option(
        0.1,
        "--jaccard-threshold",
        "-j",
        help="Jaccard threshold for combining credible sets.",
    ),
    # ABF parameters
    var_prior: float = typer.Option(
        0.2,
        "--var-prior",
        "-vp",
        help="Variance prior, by default 0.2, usually set to 0.15 for quantitative traits and 0.2 for binary traits.",
        rich_help_panel="ABF",
    ),
    # FINEMAP parameters
    n_iter: int = typer.Option(
        100000,
        "--n-iter",
        "-ni",
        help="Number of iterations.",
        rich_help_panel="FINEMAP",
    ),
    n_threads: int = typer.Option(
        1, "--n-threads", "-nt", help="Number of threads.", rich_help_panel="FINEMAP"
    ),
    # susie parameters
    max_iter: int = typer.Option(
        100,
        "--max-iter",
        "-i",
        help="Maximum number of iterations.",
        rich_help_panel="SuSie",
    ),
    estimate_residual_variance: bool = typer.Option(
        False,
        "--estimate-residual-variance",
        "-er",
        help="Estimate residual variance.",
        rich_help_panel="SuSie",
    ),
    min_abs_corr: float = typer.Option(
        0.5,
        "--min-abs-corr",
        "-mc",
        help="Minimum absolute correlation.",
        rich_help_panel="SuSie",
    ),
    convergence_tol: float = typer.Option(
        1e-3,
        "--convergence-tol",
        "-ct",
        help="Convergence tolerance.",
        rich_help_panel="SuSie",
    ),
    # RSparsePro parameters
    eps: float = typer.Option(
        1e-5, "--eps", "-e", help="Convergence criterion.", rich_help_panel="RSparsePro"
    ),
    ubound: int = typer.Option(
        100000,
        "--ubound",
        "-ub",
        help="Upper bound for convergence.",
        rich_help_panel="RSparsePro",
    ),
    cthres: float = typer.Option(
        0.7,
        "--cthres",
        "-ct",
        help="Threshold for coverage.",
        rich_help_panel="RSparsePro",
    ),
    eincre: float = typer.Option(
        1.5,
        "--eincre",
        "-ei",
        help="Adjustment for error parameter.",
        rich_help_panel="RSparsePro",
    ),
    minldthres: float = typer.Option(
        0.7,
        "--minldthres",
        "-mlt",
        help="Threshold for minimum LD within effect groups.",
        rich_help_panel="RSparsePro",
    ),
    maxldthres: float = typer.Option(
        0.2,
        "--maxldthres",
        "-mlt",
        help="Threshold for maximum LD across effect groups.",
        rich_help_panel="RSparsePro",
    ),
    varemax: float = typer.Option(
        100.0,
        "--varemax",
        "-vm",
        help="Maximum error parameter.",
        rich_help_panel="RSparsePro",
    ),
    varemin: float = typer.Option(
        1e-3,
        "--varemin",
        "-vm",
        help="Minimum error parameter.",
        rich_help_panel="RSparsePro",
    ),
    # SuSiEx parameters
    # pval_thresh: float = typer.Option(1e-5, "--pval-thresh", "-pt", help="P-value threshold for SuSiEx.", rich_help_panel="SuSiEx"),
    # maf_thresh: float = typer.Option(0.005, "--maf-thresh", "-mt", help="MAF threshold for SuSiEx.", rich_help_panel="SuSiEx"),
    mult_step: bool = typer.Option(
        False,
        "--mult-step",
        "-ms",
        help="Whether to use multiple steps in SuSiEx.",
        rich_help_panel="SuSiEx",
    ),
    keep_ambig: bool = typer.Option(
        True,
        "--keep-ambig",
        "-ka",
        help="Whether to keep ambiguous SNPs in SuSiEx.",
        rich_help_panel="SuSiEx",
    ),
    # n_threads: int = typer.Option(1, "--n-threads", "-nt", help="Number of threads.", rich_help_panel="SuSiEx"),
    min_purity: float = typer.Option(
        0.5,
        "--min-purity",
        "-mp",
        help="Minimum purity for SuSiEx.",
        rich_help_panel="SuSiEx",
    ),
    # max_iter: int = typer.Option(100, "--max-iter", "-i", help="Maximum number of iterations.", rich_help_panel="SuSiEx"),
    tol: float = typer.Option(
        1e-3, "--tol", "-t", help="Convergence tolerance.", rich_help_panel="SuSiEx"
    ),
    # MULTISUSIE parameters
    rho: float = typer.Option(
        0.75,
        "--rho",
        "-r",
        help="The prior correlation between causal variants.",
        rich_help_panel="MULTISUSIE",
    ),
    scaled_prior_variance: float = typer.Option(
        0.2,
        "--scaled-prior-variance",
        "-spv",
        help="The scaled prior variance.",
        rich_help_panel="MULTISUSIE",
    ),
    standardize: bool = typer.Option(
        False,
        "--standardize",
        "-s",
        help="Whether to standardize the summary statistics.",
        rich_help_panel="MULTISUSIE",
    ),
    pop_spec_standardization: bool = typer.Option(
        True,
        "--pop-spec-standardization",
        "-pss",
        help="Whether to use population-specific standardization.",
        rich_help_panel="MULTISUSIE",
    ),
    # estimate_residual_variance: bool = typer.Option(True, "--estimate-residual-variance", "-er", help="Estimate residual variance.", rich_help_panel="MULTISUSIE"),
    estimate_prior_variance: bool = typer.Option(
        True,
        "--estimate-prior-variance",
        "-epv",
        help="Estimate prior variance.",
        rich_help_panel="MULTISUSIE",
    ),
    estimate_prior_method: str = typer.Option(
        "early_EM",
        "--estimate-prior-method",
        "-epm",
        help="Method to estimate prior variance.",
        rich_help_panel="MULTISUSIE",
    ),
    pop_spec_effect_priors: bool = typer.Option(
        True,
        "--pop-spec-effect-priors",
        "-pesp",
        help="Whether to use population-specific effect priors.",
        rich_help_panel="MULTISUSIE",
    ),
    iter_before_zeroing_effects: int = typer.Option(
        5,
        "--iter-before-zeroing-effects",
        "-ibe",
        help="Number of iterations before zeroing out effects.",
        rich_help_panel="MULTISUSIE",
    ),
    prior_tol: float = typer.Option(
        1e-9,
        "--prior-tol",
        "-pt",
        help="Tolerance for prior variance.",
        rich_help_panel="MULTISUSIE",
    ),
    # min_abs_corr: float = typer.Option(0, "--min-abs-corr", "-mc", help="Minimum absolute correlation.", rich_help_panel="MULTISUSIE"),
    # max_iter: int = typer.Option(100, "--max-iter", "-i", help="Maximum number of iterations.", rich_help_panel="MULTISUSIE"),
    # tol: float = typer.Option(1e-3, "--tol", "-t", help="Convergence tolerance.", rich_help_panel="MULTISUSIE"),
    calculate_lambda_s: bool = typer.Option(
        False,
        "--calculate-lambda-s",
        "-cls",
        help="Calculate lambda_s parameter using estimate_s_rss function.",
    ),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Run whole fine-mapping pipeline on a list of loci."""
    setup_file_logging(log_file)
    import logging
    import sys
    from datetime import datetime

    logger = logging.getLogger("CREDTOOLS")

    loci_info = pd.read_csv(inputs, sep="\t")
    loci_info = check_loci_info(loci_info)  # Validate input data

    # Initialize overall run summary
    overall_summary = {
        "start_time": datetime.now().isoformat(),
        "total_loci": len(loci_info.groupby("locus_id")),
        "successful_loci": 0,
        "failed_loci": 0,
        "errors": [],
    }

    console = Console()

    for locus_id, locus_info in loci_info.groupby("locus_id"):
        out_dir = f"{outdir}/{locus_id}"
        os.makedirs(out_dir, exist_ok=True)

        try:
            console.print(f"[cyan]Processing locus {locus_id}...[/cyan]")
            pipeline(
                locus_info,
                outdir=out_dir,
                meta_method=meta_method,
                skip_qc=skip_qc,
                tool=tool,
                max_causal=max_causal,
                adaptive_max_causal=adaptive_max_causal,
                set_L_by_cojo=set_L_by_cojo,
                coverage=coverage,
                combine_cred=combine_cred,
                combine_pip=combine_pip,
                jaccard_threshold=jaccard_threshold,
                # susie parameters
                max_iter=max_iter,
                estimate_residual_variance=estimate_residual_variance,
                min_abs_corr=min_abs_corr,
                convergence_tol=convergence_tol,
                # ABF parameters
                var_prior=var_prior,
                # FINEMAP parameters
                n_iter=n_iter,
                n_threads=n_threads,
                # RSparsePro parameters
                eps=eps,
                ubound=ubound,
                cthres=cthres,
                eincre=eincre,
                minldthres=minldthres,
                maxldthres=maxldthres,
                varemax=varemax,
                varemin=varemin,
                # SuSiEx parameters
                mult_step=mult_step,
                keep_ambig=keep_ambig,
                min_purity=min_purity,
                tol=tol,
                # MULTISUSIE parameters
                rho=rho,
                scaled_prior_variance=scaled_prior_variance,
                standardize=standardize,
                pop_spec_standardization=pop_spec_standardization,
                estimate_prior_variance=estimate_prior_variance,
                estimate_prior_method=estimate_prior_method,
                pop_spec_effect_priors=pop_spec_effect_priors,
                iter_before_zeroing_effects=iter_before_zeroing_effects,
                prior_tol=prior_tol,
                calculate_lambda_s=calculate_lambda_s,
            )
            overall_summary["successful_loci"] += 1
            console.print(f"[green]✓ Locus {locus_id} completed successfully[/green]")

        except Exception as e:
            error_msg = f"Locus {locus_id} failed: {str(e)}"
            logger.error(error_msg)
            console.print(f"[red]✗ {error_msg}[/red]")
            overall_summary["failed_loci"] += 1
            overall_summary["errors"].append(error_msg)
            # Continue to next locus

    # Generate overall summary
    overall_summary["end_time"] = datetime.now().isoformat()
    summary_file = f"{outdir}/overall_run_summary.log"
    with open(summary_file, "w") as f:
        f.write("=== CREDTOOLS PIPELINE OVERALL SUMMARY ===\n")
        f.write(f"Start Time: {overall_summary['start_time']}\n")
        f.write(f"End Time: {overall_summary['end_time']}\n")
        f.write(f"Total Loci: {overall_summary['total_loci']}\n")
        f.write(f"Successful: {overall_summary['successful_loci']}\n")
        f.write(f"Failed: {overall_summary['failed_loci']}\n")
        if overall_summary["errors"]:
            f.write("\nError Details:\n")
            for error in overall_summary["errors"]:
                f.write(f"  - {error}\n")

    # Print final summary
    if overall_summary["failed_loci"] > 0:
        console.print(
            f"[yellow]Pipeline completed with {overall_summary['failed_loci']} failed loci[/yellow]"
        )
    else:
        console.print(
            f"[green]Pipeline completed successfully for all {overall_summary['successful_loci']} loci[/green]"
        )


@app.command()
def plot(
    input_path: str = typer.Argument(..., help="Path to QC file or locus directory"),
    plot_type: Optional[str] = typer.Option(
        None,
        "--type",
        "-t",
        help="Plot type: summary, locus_qc, locusplot, lambda_s, maf_corr, lambda_s_outliers, dentist_s_outliers, locus_pvalues, zscore_qq, ld_decay, ld_4th_moment, snp_missingness. Auto-detected if not specified.",
    ),
    output: Optional[str] = typer.Option(
        None, "--output", "-o", help="Output file path (PNG, PDF, SVG)"
    ),
    figsize_width: float = typer.Option(16, "--width", help="Figure width in inches"),
    figsize_height: float = typer.Option(
        12, "--height", help="Figure height in inches"
    ),
    dpi: int = typer.Option(300, "--dpi", help="DPI for output file"),
    include_upset: bool = typer.Option(
        True,
        "--include-upset/--no-upset",
        help="Include SNP missingness upset plot for locus plots",
    ),
):
    """Create QC plots from credtools results."""
    try:
        import matplotlib.pyplot as plt

        import credtools.plot as plot_mod
    except ImportError as e:
        console = Console()
        console.print("[red]Error: Plotting dependencies not found.[/red]")
        console.print("Please install plotting dependencies: uv add seaborn upsetplot")
        raise typer.Exit(1) from e

    console = Console()

    def _get_plotter(name: str):
        if not hasattr(plot_mod, name):
            raise AttributeError(
                f"Plot function '{name}' is not available. Update credtools or install optional plotting extras."
            )
        return getattr(plot_mod, name)

    input_path_obj = Path(input_path)

    # Auto-detect plot type if not specified
    if plot_type is None:
        if input_path_obj.is_dir():
            # Check if it's a locus directory with fine-mapping results or QC files
            if (input_path_obj / "pips.txt.gz").exists():
                plot_type = "locusplot"
            else:
                plot_type = "locus_qc"
        elif input_path_obj.name.endswith(("qc.txt.gz", "qc.txt")):
            plot_type = "summary"
        elif input_path_obj.name.endswith("compare_maf.txt.gz"):
            plot_type = "maf_corr"
        else:
            console.print(
                "[red]Cannot auto-detect plot type. Please specify --type[/red]"
            )
            console.print(
                "Available types: summary, locus_qc, locusplot, lambda_s, maf_corr, lambda_s_outliers, dentist_s_outliers, locus_pvalues, zscore_qq, ld_decay, ld_4th_moment, snp_missingness"
            )
            raise typer.Exit(1)

    console.print(f"[cyan]Creating {plot_type} plot(s)...[/cyan]")

    try:
        # Multi-panel plots
        if plot_type == "summary":
            fig = _get_plotter("plot_summary_qc")(
                qc_file=input_path,
                output_file=output,
                figsize=(figsize_width, figsize_height),
                dpi=dpi,
            )
            if output:
                console.print(f"[green]Summary QC plots saved to: {output}[/green]")
            else:
                console.print("[yellow]Displaying plots...[/yellow]")
                plt.show()

        elif plot_type == "locus_qc":
            fig = _get_plotter("plot_locus_qc")(
                locus_dir=input_path,
                output_file=output,
                figsize=(figsize_width, figsize_height),
                dpi=dpi,
                include_upset=include_upset,
            )
            if output:
                console.print(f"[green]Locus QC plots saved to: {output}[/green]")
            else:
                console.print("[yellow]Displaying plots...[/yellow]")
                plt.show()

        elif plot_type == "locusplot":
            fig = _get_plotter("plot_locusplot")(
                locus_dir=input_path,
                output_file=output,
                figsize=(figsize_width, figsize_height),
                dpi=dpi,
            )
            if output:
                console.print(f"[green]Locus plot saved to: {output}[/green]")
            else:
                console.print("[yellow]Displaying plot...[/yellow]")
                plt.show()

        # Individual plots
        else:
            # Adjust figure size for individual plots
            if figsize_width == 16 and figsize_height == 12:  # Default multi-panel size
                figsize_width = 10
                figsize_height = 6

            fig, ax = plt.subplots(figsize=(figsize_width, figsize_height))

            if plot_type == "lambda_s":
                _get_plotter("plot_lambda_s_boxplot")(input_path, ax=ax)
            elif plot_type == "maf_corr":
                _get_plotter("plot_maf_corr_barplot")(input_path, ax=ax)
            elif plot_type in ["lambda_s_outliers", "dentist_s_outliers"]:
                outlier_type = plot_type.replace("_outliers", "")
                _get_plotter("plot_outliers_barplot")(
                    input_path, outlier_type=outlier_type, ax=ax
                )
            elif plot_type == "locus_pvalues":
                _get_plotter("plot_locus_pvalues")(input_path, ax=ax)
            elif plot_type == "zscore_qq":
                _get_plotter("plot_zscore_qq")(input_path, ax=ax)
            elif plot_type == "ld_decay":
                _get_plotter("plot_ld_decay")(input_path, ax=ax)
            elif plot_type == "ld_4th_moment":
                _get_plotter("plot_ld_4th_moment")(input_path, ax=ax)
            elif plot_type == "snp_missingness":
                _get_plotter("plot_snp_missingness_upset")(input_path, ax=ax)
            else:
                console.print(f"[red]Unknown plot type: {plot_type}[/red]")
                console.print(
                    "Available types: summary, locus_qc, locusplot, lambda_s, maf_corr, lambda_s_outliers, dentist_s_outliers, locus_pvalues, zscore_qq, ld_decay, ld_4th_moment, snp_missingness"
                )
                raise typer.Exit(1)

            plt.title(f"{plot_type.replace('_', ' ').title()} Plot")
            plt.tight_layout()

            if output:
                plt.savefig(output, dpi=dpi, bbox_inches="tight")
                console.print(f"[green]Plot saved to: {output}[/green]")
            else:
                console.print("[yellow]Displaying plot...[/yellow]")
                plt.show()

    except Exception as e:
        console.print(f"[red]Error creating {plot_type} plot: {e}[/red]")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
