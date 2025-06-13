# CREDTOOLS Web Visualization

This module provides a web-based visualization interface for CREDTOOLS (Credible Set Tools) results.

## Installation

To use the web visualization features, you need to install additional dependencies:

```bash
pip install credtools[web]
# or
pip install dash dash-bootstrap-components dash-mantine-components plotly
```

## Usage

### Starting the Web Interface

```bash
# Basic usage - uses current directory as data directory
credtools web

# Specify data directory and custom options
credtools web /path/to/data --webdata-dir custom_webdata --port 8080

# Force regeneration of web data with specific loci files
credtools web /path/to/data \
  --allmeta-loci data/real/meta/all/all_meta_loci_sig.txt \
  --popumeta-loci data/real/meta/ancestry/loci_info_sig.txt \
  --nometa-loci data/real/all_loci_list_sig.txt \
  --force-regenerate
```

### Command Options

- `data_dir`: Base directory containing fine-mapping data (default: current directory)
- `--webdata-dir, -w`: Directory for processed web data (default: "webdata")
- `--allmeta-loci, -a`: Path to allmeta loci info file
- `--popumeta-loci, -p`: Path to popumeta loci info file  
- `--nometa-loci, -n`: Path to nometa loci info file
- `--force-regenerate, -f`: Force regeneration of web data
- `--threads, -t`: Number of threads for data processing (default: 10)
- `--port`: Port to run web server on (default: 8080)
- `--host`: Host to bind web server to (default: "0.0.0.0")
- `--debug`: Run in debug mode

## Features

### Home Page
- Overview of all loci with filtering by meta-analysis method and fine-mapping tool
- Interactive plots showing credible set statistics
- Sortable table with locus information
- Clickable links to detailed locus views

### Locus Page
- Detailed visualization for individual loci
- LocusZoom-style plots showing association signals and LD structure
- Quality control metrics display
- Multiple fine-mapping tool comparison

## Data Structure

The web interface expects processed data in the following structure:

```
webdata/
├── all_loci_info.txt          # Summary information for all loci
├── allmeta/                   # All-ancestry meta-analysis results
│   └── [locus_id]/
│       ├── [pop]_[cohort].res.gz
│       ├── qc.txt.gz
│       └── ...
├── popumeta/                  # Population-specific meta-analysis results
│   └── [locus_id]/
│       └── ...
└── nometa/                    # No meta-analysis results
    └── [locus_id]/
        └── ...
```

## Development

The web module consists of:

- `export.py`: Data processing and export functionality
- `app.py`: Main Dash application
- `pages/`: Individual page components
  - `home.py`: Home page with overview
  - `locus.py`: Detailed locus view

## Troubleshooting

### Missing Dependencies
If you see import errors, ensure web dependencies are installed:
```bash
pip install credtools[web]
```

### No Data Found
Ensure you have run fine-mapping analysis and have the expected data structure. The web command will try to find default loci files in standard locations, or you can specify them explicitly.

### Port Already in Use
Use a different port with `--port` option:
```bash
credtools web --port 8081
``` 