# Changelog

## [0.0.38] (2025-09-21)

### Added
- **FINEMAP timeout control**: `credtools finemap` now accepts `--timeout-minutes` (default 30) to cap per-locus runs, with timeout failures reported in `run_summary.log`.

### Changed
- **External tool execution**: FINEMAP wrapper and tool manager enforce the timeout and append `[timeout]` markers to logs when the limit is exceeded, ensuring stalled runs surface immediately.

## [0.0.37] (2025-09-20)

### Fixed
- **LD matrix ingestion**: Replace NaN entries with zeros when loading from lower-triangle text or `.npz` files, and add regression tests to guard behavior.
- **QC pipeline robustness**: Locus-level exceptions no longer abort `credtools qc`; failures are logged per locus, success/failure counts are summarized to `qc_run_summary.log`, and CLI feedback reflects mixed outcomes.
- **Error messaging**: Intersecting sumstats/LD now reports the offending locus when common variants are missing to simplify debugging.

## [0.0.36] (2025-09-18)

### Changed
- **CLI preprocessing flows**: `credtools munge` and `credtools chunk` now accept
  comma- or newline-separated input paths in addition to TSV config files,
  with friendlier validation and logging when LD references are omitted.
- **Preprocessing outputs**: Munge results retain the sample-size `N` column in
  the standard export schema and downstream validators require it.
- **Optional plotting dependency**: Centralized the `upsetplot` import so the
  plotting module gatekeeps the dependency once and falls back gracefully when
  the package is absent.


## [0.0.35] (2025-09-17)

### Added
- **Comprehensive Plotting Module**: New visualization capabilities for QC results
  - Added `credtools.plot` module with publication-quality plotting functions
  - **Summary QC Plots** (2x2 layout from `qc.txt.gz`):
    - Lambda-s distribution boxplot by cohort
    - MAF correlation barplot between summary statistics and LD reference
    - Lambda-s outliers count barplot by cohort
    - Dentist-s outliers count barplot by cohort
  - **Locus-specific Plots** (2x2 layout from individual QC files):
    - Locus p-value plot with credible set annotations (from `expected_z.txt.gz`)
    - Observed vs expected z-scores QQ plot (from `expected_z.txt.gz`)
    - LD decay curve plots by cohort (from `ld_decay.txt.gz`)
    - LD 4th moment boxplots by cohort (from `ld_4th_moment.txt.gz`)
    - SNP missingness upset plot showing overlap patterns (from `snp_missingness.txt.gz`)
  - **CLI Command**:
    - `credtools plot`: Unified plotting command with auto-detection of plot type
    - Supports all plot types: summary (2x2), locus (2x2), and individual plots
    - Smart auto-detection based on input path (directory for locus, qc.txt.gz for summary)
  - **Features**:
    - Population-aware color schemes for consistent visualization
    - Support for PNG, PDF, SVG output formats
    - Customizable figure sizes and DPI settings
    - Graceful handling of missing files and optional dependencies
    - Professional styling with seaborn integration

### Dependencies
- Added `seaborn>=0.11.0` for enhanced statistical plotting
- Added `upsetplot>=0.6.0` for intersection visualization (optional for SNP missingness plots)

### Enhanced
- QC workflow now generates comprehensive visualizations alongside numerical results
- All plotting functions integrate seamlessly with existing QC output file formats
- Error handling with informative messages for missing dependencies or data files

## [0.0.34] (2025-09-14)

### Fixed
- **CI Test Suite**: Fixed failing CI tests and errors
  - Fixed `TypeError` in test fixtures by adding required `locus_start` and `locus_end` parameters to `Locus` constructor calls
  - Removed unused `toml` import that was causing `ModuleNotFoundError` in integration tests
  - Fixed linting issues: converted f-strings without placeholders to regular strings
  - Added missing docstrings for mock classes in QC module
  - All 70 tests now pass successfully

### Improved
- Enhanced test coverage and reliability of CI pipeline
- Improved code quality with resolved linting and documentation issues

## [0.0.33] (2025-09-12)

### Added
- **QC Summary Statistics**: Enhanced quality control with comprehensive summary reports
  - Added locus-level QC summary files (`qc.txt.gz`) in each locus directory
  - Added global QC summary file (`qc.txt.gz`) in output root directory
  - Summary includes key QC metrics: SNP counts, significance thresholds, MAF correlations, lambda-s values, and outlier counts
  - Configurable thresholds for flip detection, lambda-s outliers, and Dentist-S outliers
  - New `locus_qc_summary()` function to generate summary statistics from detailed QC results
  - Enhanced QC metrics:
    - `n_1e-5`, `n_5e-8`: Count of SNPs below significance thresholds
    - `maf_corr`: Correlation between summary statistics and LD reference MAF
    - `n_flip`: Count of potential allele flips (logLR > 2 AND |z| > 2)
    - `n_lambda_s_outlier`: Count of lambda-s outliers (|z_std_diff| > 3)
    - `n_dentist_s_outlier`: Count of Dentist-S outliers (-log10p ≥ 4 AND r² ≥ 0.6)

### Changed
- Modified `locus_qc()` function to accept threshold parameters for outlier detection
- Updated `qc_locus_cli()` to return both locus ID and summary statistics
- Enhanced `loci_qc()` function to aggregate and save global QC summary across all loci
- Improved QC workflow output structure with hierarchical summary files

## [0.0.32] (2025-09-12)

### Changed
- **BREAKING**: Removed `strategy` parameter from fine-mapping interface
  - Fine-mapping strategy is now automatically determined based on tool type and data structure
  - Multi-input tools (susiex, multisusie) automatically process all loci together
  - Single-input tools automatically combine results when multiple loci are provided
  - Added deprecation warning for backward compatibility
- Enhanced CLI with enum validation for combination methods
  - Added `CombineCred` enum for credible set combination methods (union, intersection, cluster)
  - Added `CombinePIP` enum for PIP combination methods (max, min, mean, meta)
  - Improved input validation and auto-completion support

### Removed
- Web visualization feature moved to v2 (will be available in future release)
  - Removed `credtools web` command documentation
  - Removed web-related installation instructions
  - Removed web tutorial files and examples
  - Updated all workflow examples to reference output files instead

### Improved
- Simplified user interface with automatic strategy selection
- Better CLI help with enum option display
- Updated documentation to reflect streamlined workflow

## [0.0.31] (2025-09-11)

### Fixed
- CI error

## [0.0.30] (2025-09-11)

### Added
- ABF+COJO
- adaptive causal

## [0.0.28] (2025-06-13)

### Added
- add api docs

## [0.0.27] (2025-06-12)

### Added
- add set_L_by_cojo to cli:pipeline

## [0.0.26] (2025-06-02)

### Added
- add web visualization

## [0.0.25] (2025-06-02)

### Added
- add tutorial

## [0.0.23] (2025-02-01)

### Fixed
- fix finemap cred bug

## [0.0.21] (2025-01-20)

### Fixed
- fix no install error for carma

## [0.0.20] (2025-01-20)

### Fixed
- fix zero maf in finemap

## [0.0.19] (2025-01-20)

### Added
- qc support for multiprocessing

## [0.0.18] (2025-01-19)

### Fixed
- fix the bug of no credible set

## [0.0.17] (2025-01-18)

### Added
- support for multiprocessing
- add progress bar

## [0.0.16] (2025-01-18)

### Added
- support for sumstats.gz and ldmap.gz


## [0.0.15] (2024-12-17)

### Added
- cli args

## [0.0.14] (2024-12-16)

### Added
- cli

## [0.0.13] (2024-12-16)

### Added
- pipeline

## [0.0.12] (2024-12-15)

### Added
- ensemble fine-mapping

## [0.0.11] (2024-12-15)

### Added
- multisusie

## [0.0.10] (2024-12-15)

### Added
- susiex
- Rsparseld
- CARMA

## [0.0.9] (2024-10-21)

### Added
- abf
- susie
- finemap

## [0.0.8] (2024-10-10)

### Added
- load ld matrix and ld map
- munge sumstat
- example data

## [0.0.7] (2024-10-09)

### Added
- test for ldmatrix

## [0.0.6] (2024-10-09)

### Added
- functions for load LD
- test for ColName


## [0.0.5] (2024-10-08)

* First release on PyPI.
