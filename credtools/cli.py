"""Console script for credtools."""

import json
import logging
import logging.handlers
import os
from enum import Enum
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

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
from credtools.locus import load_locus_set
from credtools.meta import meta_loci
from credtools.qc import loci_qc

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)


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


@app.command(
    name="munge",
    help="Reformat and standardize GWAS summary statistics.",
)
def run_munge(
    input_files: str = typer.Argument(
        ...,
        help="Input sumstats file(s). Can be a single file, comma-separated list, or config file.",
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
    """Reformat and standardize GWAS summary statistics using smunger integration."""
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
    console.print(f"[cyan]Munging summary statistics...[/cyan]")

    # Parse input files
    if "," in input_files:
        # Comma-separated list
        file_list = [f.strip() for f in input_files.split(",")]
        input_dict = {Path(f).stem: f for f in file_list}
    elif input_files.endswith((".json", ".yaml", ".yml")):
        # Configuration file with file mappings
        import json

        with open(input_files) as f:
            input_dict = json.load(f)
    else:
        # Single file
        input_dict = input_files

    # Create interactive config if requested
    if interactive_config and isinstance(input_dict, dict):
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

        # Validate results
        validation = validate_munged_files(result)

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


@app.command(
    name="chunk",
    help="Identify independent loci and chunk summary statistics.",
)
def run_chunk(
    input_files: str = typer.Argument(
        ...,
        help="Munged sumstats file(s). Can be single file, comma-separated list, or config file.",
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
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Identify independent loci and chunk summary statistics for fine-mapping."""
    setup_file_logging(log_file)

    try:
        from credtools.preprocessing import chunk_sumstats, identify_independent_loci
        from credtools.preprocessing.chunk import create_loci_list_for_credtools
    except ImportError as e:
        console = Console()
        console.print("[red]Error: Preprocessing module not found.[/red]")
        raise typer.Exit(1) from e

    from pathlib import Path

    console = Console()
    console.print(f"[cyan]Identifying independent loci...[/cyan]")

    # Parse input files
    if "," in input_files:
        file_list = [f.strip() for f in input_files.split(",")]
        input_dict = {Path(f).stem.replace(".munged", ""): f for f in file_list}
    elif input_files.endswith((".json", ".yaml", ".yml")):
        import json

        with open(input_files) as f:
            input_dict = json.load(f)
    else:
        input_dict = input_files

    try:
        # Identify loci
        loci_df = identify_independent_loci(
            sumstats_files=input_dict,
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

        # Chunk summary statistics
        console.print(f"[cyan]Chunking {len(loci_df)} loci...[/cyan]")
        chunk_info_df = chunk_sumstats(
            loci_df=loci_df,
            sumstats_files=input_dict,
            output_dir=os.path.join(output_dir, "chunks"),
            threads=threads,
        )

        # Create credtools-compatible loci list
        loci_list_file = os.path.join(output_dir, "loci_list.txt")
        credtools_df = create_loci_list_for_credtools(
            chunk_info_df=chunk_info_df, output_file=loci_list_file
        )

        # Print summary
        console.print(f"[green]Successfully processed {len(loci_df)} loci[/green]")
        console.print(f"[green]Generated {len(chunk_info_df)} chunked files[/green]")
        console.print(f"[green]Credtools loci list: {loci_list_file}[/green]")

    except Exception as e:
        console.print(f"[red]Error during chunking: {e}[/red]")
        raise typer.Exit(1)


@app.command(
    name="prepare",
    help="Prepare LD matrices and final fine-mapping inputs.",
)
def run_prepare(
    chunk_info: str = typer.Argument(..., help="Chunk info file from 'chunk' command."),
    genotype_config: str = typer.Argument(
        ...,
        help="Genotype configuration file (JSON) mapping ancestries to file prefixes.",
    ),
    output_dir: str = typer.Argument(..., help="Output directory for prepared files."),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
    ld_format: str = typer.Option(
        "plink", "--ld-format", "-f", help="LD computation format (plink/vcf)."
    ),
    keep_intermediate: bool = typer.Option(
        False, "--keep-intermediate", "-k", help="Keep intermediate files."
    ),
    log_file: Optional[str] = typer.Option(
        None, "--log-file", "-l", help="Log output to specified file."
    ),
):
    """Prepare LD matrices and final fine-mapping input files."""
    setup_file_logging(log_file)
    try:
        from credtools.preprocessing import prepare_finemap_inputs
        from credtools.preprocessing.chunk import create_loci_list_for_credtools
    except ImportError as e:
        console = Console()
        console.print("[red]Error: Preprocessing module not found.[/red]")
        raise typer.Exit(1) from e

    import json

    console = Console()
    console.print(f"[cyan]Preparing fine-mapping inputs...[/cyan]")

    # Load chunk info
    if not os.path.exists(chunk_info):
        console.print(f"[red]Chunk info file not found: {chunk_info}[/red]")
        raise typer.Exit(1)

    chunk_info_df = pd.read_csv(chunk_info, sep="\t")
    console.print(f"Loaded {len(chunk_info_df)} chunks")

    # Load genotype configuration
    if not os.path.exists(genotype_config):
        console.print(f"[red]Genotype config file not found: {genotype_config}[/red]")
        raise typer.Exit(1)

    with open(genotype_config) as f:
        genotype_files = json.load(f)

    console.print(f"Genotype files for {len(genotype_files)} ancestries")

    try:
        # Prepare files
        prepared_df = prepare_finemap_inputs(
            chunk_info_df=chunk_info_df,
            genotype_files=genotype_files,
            output_dir=output_dir,
            threads=threads,
            ld_format=ld_format,
            keep_intermediate=keep_intermediate,
        )

        # Create final loci list for credtools
        final_loci_file = os.path.join(output_dir, "final_loci_list.txt")
        final_df = create_loci_list_for_credtools(
            chunk_info_df=prepared_df, output_file=final_loci_file
        )

        # Print summary
        console.print(f"[green]Successfully prepared {len(prepared_df)} files[/green]")
        console.print(f"[green]Final loci list: {final_loci_file}[/green]")

        # Print ancestry summary
        ancestry_summary = prepared_df.groupby("ancestry").size()
        for ancestry, count in ancestry_summary.items():
            console.print(f"  {ancestry}: {count} loci")

    except Exception as e:
        console.print(f"[red]Error during preparation: {e}[/red]")
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
    loci_qc(inputs, outdir, threads)


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
            "When using single-input tools with multiple loci, results are automatically combined."
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
    combine_cred: CombineCred = typer.Option(
        CombineCred.union, 
        "--combine-cred", 
        "-cc", 
        help="Method to combine credible sets when using single-input tools with multiple loci."
    ),
    combine_pip: CombinePIP = typer.Option(
        CombinePIP.max, 
        "--combine-pip", 
        "-cp", 
        help="Method to combine PIPs when using single-input tools with multiple loci."
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
    loci_info = pd.read_csv(inputs, sep="\t")
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

    with progress:
        task = progress.add_task("[cyan]Fine-mapping loci...", total=total_loci)

        for locus_id, locus_info in locus_groups:
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
                combine_cred=combine_cred,
                combine_pip=combine_pip,
                jaccard_threshold=jaccard_threshold,
                # susie parameters
                max_iter=max_iter,
                estimate_residual_variance=estimate_residual_variance,
                min_abs_corr=min_abs_corr,
                convergence_tol=convergence_tol,
            )
            out_dir = f"{outdir}/{locus_id}"
            os.makedirs(out_dir, exist_ok=True)
            creds.pips.to_csv(f"{out_dir}/pips.txt", sep="\t", header=False, index=True)
            with open(f"{out_dir}/creds.json", "w") as f:
                json.dump(creds.to_dict(), f, indent=4)

            progress.advance(task)


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
             "multi-input tools process all loci together."
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
        help="Method to combine credible sets when using single-input tools with multiple loci."
    ),
    combine_pip: CombinePIP = typer.Option(
        CombinePIP.max, 
        "--combine-pip", 
        "-cp", 
        help="Method to combine PIPs when using single-input tools with multiple loci."
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
    loci_info = pd.read_csv(inputs, sep="\t")
    for locus_id, locus_info in loci_info.groupby("locus_id"):
        out_dir = f"{outdir}/{locus_id}"
        os.makedirs(out_dir, exist_ok=True)
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


if __name__ == "__main__":
    app()
