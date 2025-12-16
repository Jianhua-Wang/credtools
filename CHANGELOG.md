# Changelog

## [0.3.1] (2025-12-16)

### Fixed
- **QC threshold parameter consistency**: Unified default threshold values between CLI and internal functions
  - Changed CLI defaults: `logLR_threshold` (1.0→2.0), `z_threshold` (1.0→2.0), `r_threshold` (0.1→0.8)
  - Fixed missing threshold parameters in `remove_outliers_and_rerun_qc()` causing inconsistent outlier detection
  - Cleaned data now correctly shows `n_lambda_s_outlier = 0` when using same thresholds for detection and re-QC

## [0.3.0] (2025-12-15)

### Fixed
- **QC outlier removal output structure**: Reorganized cleaned data output to improve clarity and downstream usability
  - **BREAKING**: All cleaned data now organized under `{out_dir}/cleaned/` directory instead of `{out_dir}/{locus_id}/cleaned/`
  - **Fixed missing files for loci without outliers**: Previously, when no outliers were detected, cleaned data files were not saved, causing `cleaned_loci_info.txt.gz` paths to be invalid. Now all loci have cleaned data saved (identical to original if no outliers detected)
  - **Added missing cleaned QC summaries**:
    - Per-locus cleaned QC summary: `{out_dir}/cleaned/{locus_id}/qc.txt.gz`
    - Global cleaned QC summary: `{out_dir}/cleaned/qc_cleaned.txt.gz`

### Changed
- **New output directory structure for cleaned data**:
  ```
  {out_dir}/
  ├── qc.txt.gz                              # Original QC summary
  ├── {locus_id}/                            # Original QC results
  │   ├── qc.txt.gz
  │   ├── expected_z.txt.gz
  │   └── ...
  └── cleaned/                               # NEW: All cleaned results
      ├── qc_cleaned.txt.gz                  # Global cleaned QC summary
      ├── outlier_removal_summary.txt.gz     # Global outlier removal summary
      ├── cleaned_loci_info.txt.gz           # Cleaned loci info for downstream
      └── {locus_id}/                        # Per-locus cleaned data
          ├── {prefix}.sumstats.gz
          ├── {prefix}.ld.npz
          ├── {prefix}.ldmap.gz
          ├── qc.txt.gz                      # Per-locus cleaned QC
          ├── expected_z.txt.gz
          └── outlier_removal_summary.txt.gz
  ```

### Migration Guide

#### For users of `credtools qc --remove-outlier`

**Old paths** (no longer valid):
```
{out_dir}/qc.txt.gz
{out_dir}/outlier_removal_summary.txt.gz
{out_dir}/cleaned_loci_info.txt.gz
{out_dir}/{locus_id}/cleaned/{prefix}.sumstats.gz
{out_dir}/{locus_id}/cleaned/outlier_removal_summary.txt.gz
```

**New paths**:
```
{out_dir}/qc.txt.gz                                    # Original QC (unchanged)
{out_dir}/cleaned/qc_cleaned.txt.gz                    # NEW
{out_dir}/cleaned/outlier_removal_summary.txt.gz       # MOVED
{out_dir}/cleaned/cleaned_loci_info.txt.gz             # MOVED
{out_dir}/cleaned/{locus_id}/{prefix}.sumstats.gz      # MOVED
{out_dir}/cleaned/{locus_id}/qc.txt.gz                 # NEW
{out_dir}/cleaned/{locus_id}/outlier_removal_summary.txt.gz  # MOVED
```

**For downstream fine-mapping**:
- Use `{out_dir}/cleaned/cleaned_loci_info.txt.gz` as input (path updated but still valid)
- All paths in this file now correctly point to `{out_dir}/cleaned/{locus_id}/`

### Benefits
- **Clearer organization**: Original and cleaned data completely separated
- **Reliable paths**: All cleaned data files guaranteed to exist, even when no outliers detected
- **Better QC tracking**: Can compare `qc.txt.gz` vs `cleaned/qc_cleaned.txt.gz` to assess cleaning impact
- **Consistent structure**: All cleaned outputs follow same directory pattern

## [0.2.0] (2025-11-25)

### BREAKING CHANGES
- **Renamed `min_abs_corr` parameter to `purity`** across all CLI commands and Python API
  - CLI: `--min-abs-corr` flag replaced with `--purity` (short form: `-p`)
  - Python API: `min_abs_corr=` parameter renamed to `purity=` in all function signatures
  - **Migration required**: Update all scripts and workflows to use the new parameter name
- **Changed default purity threshold from 0.5/0.1 to 0.0** (no filtering by default)
  - Previous behavior: SuSiE family tools had `min_abs_corr=0.1` (function) or 0.5 (CLI)
  - New behavior: `purity=0.0` means no filtering unless explicitly set by user
  - Other tools (ABF, FINEMAP, RSparsePro) previously had no purity filtering
  - **Migration impact**: Users must explicitly set `--purity` to filter credible sets

### Added
- **Unified purity calculation and filtering for all fine-mapping tools**
  - All 7 tools now calculate and store purity values in CredibleSet objects
  - Purity definition: minimum absolute LD correlation between all variant pairs in a credible set
  - Tools now include: ABF, ABF+COJO, FINEMAP, RSparsePro (previously unsupported)
  - Multi-ancestry purity: Element-wise maximum LD across populations, then minimum
- **Centralized purity filtering logic**
  - New `filter_credset_by_purity()` function in `credibleset.py`
  - Filtering applied in `fine_map()` and `combine_creds()` functions
  - Consistent behavior across all tools and analysis strategies
- **Enhanced `combine_creds()` function**
  - Added `min_purity` parameter for filtering combined credible sets
  - Purity filtering applied after merging credible sets from multiple tools/loci

### Changed
- **More intuitive parameter naming**: "purity" more clearly describes the metric than "min_abs_corr"
- **Opt-in filtering model**: Users choose when to apply purity filtering (default: no filtering)
- **Transparent purity values**: All tools now expose purity calculations for user inspection

### Migration Guide

#### For CLI Users
```bash
# Old command (no longer works)
credtools finemap --min-abs-corr 0.5 [other options]

# New command
credtools finemap --purity 0.5 [other options]

# Or use short form
credtools finemap -p 0.5 [other options]
```

#### For Python API Users
```python
# Old code (no longer works)
from credtools import fine_map
result = fine_map(locus_set, tool="susie", min_abs_corr=0.5)

# New code
from credtools import fine_map
result = fine_map(locus_set, tool="susie", purity=0.5)
```

#### Recommended Actions
1. **Update all scripts**: Replace `--min-abs-corr` with `--purity` and `min_abs_corr=` with `purity=`
2. **Review filtering needs**: Since default changed to 0.0 (no filtering), explicitly set `--purity` if needed
3. **For old SuSiE behavior**: Use `--purity 0.1` to match previous SuSiE function default
4. **Check purity values**: All tools now store purity - inspect `credset.purity` for quality assessment

### Benefits
- **Consistency**: Same purity filtering available for all 7 fine-mapping tools
- **Flexibility**: Users control filtering threshold per analysis
- **Transparency**: Purity values accessible for all credible sets
- **Better quality control**: Unified metric for assessing credible set reliability

## [0.1.0] (2025-11-25)

### Added
- **Adaptive max_causal for multi-input tools**: Extended adaptive max_causal parameter tuning to multi-input fine-mapping tools (multisusie and susiex)
  - Added `_adaptive_fine_map_multi()` function for LocusSet-level adaptive logic
  - Same algorithm as single-input tools: increase max_causal by 5 (up to 20) when saturated, decrease by 1 (down to 1) on failure
  - Multi-input tools now support `--adaptive-max-causal` flag for automatic parameter optimization
  - Unified adaptive strategy across all supported fine-mapping tools

### Changed
- **SuSiEx mult_step default**: Changed default value from `True` to `False` to avoid conflicts with credtools' adaptive logic
  - Users can still enable SuSiEx's internal multi-step refinement by explicitly setting `--mult-step`
  - Updated documentation to clarify interaction between `--mult-step` and `--adaptive-max-causal`

### Enhanced
- **Documentation improvements**:
  - Updated `adaptive_max_causal` parameter documentation to include multi-input tool support
  - Enhanced CLI help text for `--mult-step` option with usage guidance
  - Added detailed docstrings explaining adaptive logic for both single-locus and LocusSet processing

### Technical Details
- Single-input tools (finemap, susie, rsparsepro): Adaptive logic applies per-locus
- Multi-input tools (multisusie, susiex): Adaptive logic applies to entire LocusSet
- Backward compatible: `adaptive_max_causal` defaults to `False`, existing workflows unchanged

## [0.0.41] (2025-11-22)

### Fixed
- **Multi-cohort CS purity calculation**: Fix ValueError when calculating credible set purity across cohorts with non-overlapping SNP sets. The function now takes the union of all CS SNPs across cohorts and expands each cohort's LD matrix to the union size (following MultiSuSiE's approach), with missing SNPs represented as uncorrelated (LD=0, diagonal=1). This resolves broadcasting errors like "operands could not be broadcast together with shapes (96,96) (95,95)".

## [0.0.40] (2025-10-30)
### Added
- add cred and PIP of each cohort in multi cohort fine-mapping
- output zero credible sets when n_sig 5e-8 = 0

### Fixed
- fix multisusie doesn't filter purity issue

## [0.0.39] (2025-09-22)

### Added
- **Parallel FINEMAP execution**: `credtools finemap` accepts `--processes/-np` to fan out loci across worker processes, keeping sequential behaviour available with the default of 1.

### Changed
- **Per-locus task coordination**: Fine-mapping runs now stream through a worker-safe task wrapper, preserving rich progress reporting, capturing tracebacks, and writing locus outputs from whichever worker completes first.
- **CLI artefacts**: Run summaries and combined credible set exports remain deterministic under parallel execution, with parameters.json recording the requested worker count.

### Tests
- Added regression coverage for the CLI multiprocessing path by mocking the worker pool to assert identical outputs for 1 vs N processes.

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
