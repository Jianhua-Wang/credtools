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
# Run complete fine-mapping pipeline
credtools pipeline input_loci.txt output_dir/

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

## Documentation

For detailed documentation, see <https://Jianhua-Wang.github.io/credtools>

## Web Visualization

The web module (`credtools.web`) provides interactive visualization of fine-mapping results. See [credtools/web/README.md](credtools/web/README.md) for detailed usage instructions.
