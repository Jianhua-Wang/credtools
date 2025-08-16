# credtools


[![pypi](https://img.shields.io/pypi/v/credtools.svg)](https://pypi.org/project/credtools/)
[![python](https://img.shields.io/pypi/pyversions/credtools.svg)](https://pypi.org/project/credtools/)
[![Build Status](https://github.com/Jianhua-Wang/credtools/actions/workflows/dev.yml/badge.svg)](https://github.com/Jianhua-Wang/credtools/actions/workflows/dev.yml)
[![codecov](https://codecov.io/gh/Jianhua-Wang/credtools/branch/main/graphs/badge.svg)](https://codecov.io/github/Jianhua-Wang/credtools)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



Multi-ancestry fine-mapping pipeline with interactive web visualization.


* Documentation: <https://Jianhua-Wang.github.io/credtools>
* GitHub: <https://github.com/Jianhua-Wang/credtools>
* PyPI: <https://pypi.org/project/credtools/>
* Free software: MIT


## Features

- **Whole-genome preprocessing**: Start from raw GWAS summary statistics and genotype data
  - Standardize and munge summary statistics from various formats
  - Prepare LD matrices and fine-mapping inputs automatically
- **Multi-ancestry fine-mapping**: Support for multiple fine-mapping tools (SuSiE, FINEMAP, etc.)
- **Meta-analysis capabilities**: Combine results across populations and cohorts
- **Quality control**: Built-in QC metrics and visualizations
- **Interactive web interface**: Explore results through a modern web dashboard
- **Command-line interface**: Easy-to-use CLI for all operations

## Installation

### Basic Installation
```bash
pip install credtools
```

### With Web Visualization
To use the interactive web interface, install with web dependencies:
```bash
pip install credtools[web]
```

## Quick Start

### Command Line Usage

```bash
# Complete workflow: from whole-genome data to fine-mapping results
# Step 1: Standardize summary statistics
credtools munge raw_gwas_eur.txt,raw_gwas_asn.txt output/munged/

# Step 2: Identify independent loci and chunk data
credtools chunk output/munged/*.munged.txt.gz output/chunks/

# Step 3: Prepare LD matrices and final inputs
credtools prepare output/chunks/chunk_info.txt genotype_config.json output/prepared/

# Step 4: Run fine-mapping pipeline
credtools pipeline output/prepared/final_loci_list.txt output/results/

# Step 5: Launch web visualization interface
credtools web output/results --port 8080
```

### Advanced Usage

```bash
# Launch web visualization interface
credtools web /path/to/results --port 8080

# View specific loci files
credtools web /path/to/data \
  --allmeta-loci data/allmeta_loci.txt \
  --popumeta-loci data/popumeta_loci.txt \
  --nometa-loci data/nometa_loci.txt
```

### Web Interface

The web interface provides:
- **Home page**: Overview of all loci with interactive filtering
- **Locus pages**: Detailed views with LocusZoom-style plots
- **Quality control**: Comprehensive QC metrics and visualizations
- **Multi-tool comparison**: Compare results across different fine-mapping methods

Access the web interface at `http://localhost:8080` after running `credtools web`.

## Preprocessing Workflow

credtools now supports starting from whole-genome summary statistics and genotype data, eliminating the need for manual preprocessing:

### Step 1: Munge Summary Statistics (`credtools munge`)
- **Purpose**: Standardize and clean GWAS summary statistics from various formats
- **Features**: 
  - Automatic header detection and mapping
  - Data validation and quality control
  - Support for multiple file formats
- **Input**: Raw GWAS files with various column headers
- **Output**: Standardized `.munged.txt.gz` files

### Step 2: Chunk Loci (`credtools chunk`)  
- **Purpose**: Identify independent loci and create regional chunks for fine-mapping
- **Features**:
  - Distance-based independent SNP identification
  - Cross-ancestry loci coordination
  - Configurable significance thresholds
- **Input**: Munged summary statistics files
- **Output**: Locus-specific chunked files and metadata

### Step 3: Prepare Inputs (`credtools prepare`)
- **Purpose**: Generate LD matrices and final fine-mapping input files
- **Features**:
  - LD matrix computation from genotype data
  - Variant intersection and quality control  
  - Multi-threaded processing
- **Input**: Chunked files + genotype data configuration
- **Output**: credtools-ready input files (`.sumstats.gz`, `.ld.npz`, `.ldmap.gz`)

### Multi-Ancestry Support
- **Consistent loci definition**: Union approach across ancestries
- **Flexible input formats**: Support for various GWAS summary statistics formats
- **Coordinated processing**: Ensure compatibility across populations

## Documentation

For detailed documentation, see <https://Jianhua-Wang.github.io/credtools>

## Web Visualization

The web module (`credtools.web`) provides interactive visualization of fine-mapping results. See [credtools/web/README.md](credtools/web/README.md) for detailed usage instructions.
