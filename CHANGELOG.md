# Changelog

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
