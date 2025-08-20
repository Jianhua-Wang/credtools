# Web Visualization

The `credtools web` command launches an interactive web dashboard for exploring and visualizing fine-mapping results. This provides a user-friendly interface for examining posterior inclusion probabilities, credible sets, and comparing results across loci and ancestries.

## Overview

The web interface transforms static fine-mapping results into an interactive exploration platform. The web visualization:

- **Provides interactive plots** of posterior inclusion probabilities and association signals
- **Enables locus browsing** with dynamic filtering and selection
- **Compares results** across different methods, ancestries, and meta-analysis strategies
- **Displays credible sets** with customizable confidence levels
- **Offers data export** capabilities for figures and tables
- **Supports real-time analysis** with parameter adjustment
- **Handles large datasets** efficiently with optimized data loading

## When to Use

Use `credtools web` when you have:

- Completed fine-mapping analyses that need exploration
- Results from multiple methods, ancestries, or strategies to compare
- Need to create publication-quality figures
- Want to share results with collaborators interactively
- Large numbers of loci requiring systematic review
- Need to identify high-priority variants for follow-up

## Basic Usage

### Launch with Default Settings

```bash
credtools web finemap_results/
```

### Custom Port and Host

```bash
credtools web finemap_results/ --port 8080 --host 0.0.0.0
```

### Pre-process Data for Web

```bash
credtools web data_directory/ \
  --allmeta-loci meta/all_meta_loci.txt \
  --popumeta-loci meta/pop_meta_loci.txt \
  --nometa-loci nometa/loci_list.txt \
  --force-regenerate
```

## Command Options

### Data and Processing Options

| Option | Description | Default |
|--------|-------------|---------|
| `--webdata-dir` / `-w` | Directory for processed web data | webdata |
| `--allmeta-loci` / `-a` | Path to allmeta loci info file | None |
| `--popumeta-loci` / `-p` | Path to popumeta loci info file | None |
| `--nometa-loci` / `-n` | Path to nometa loci info file | None |
| `--force-regenerate` / `-f` | Force regeneration of web data | False |
| `--threads` / `-t` | Threads for data processing | 10 |

### Server Options

| Option | Description | Default |
|--------|-------------|---------|
| `--port` | Port to run web server | 8080 |
| `--host` | Host to bind server | 0.0.0.0 |
| `--debug` | Run in debug mode | False |

## Data Preparation

The web interface requires processed data files optimized for interactive visualization:

### Automatic Data Detection

If no loci files are specified, credtools web automatically searches for:
- `data/real/meta/all/all_meta_loci_sig.txt` (allmeta results)
- `data/real/meta/ancestry/loci_info_sig.txt` (popumeta results)  
- `data/real/all_loci_list_sig.txt` (nometa results)

### Manual Data Specification

For custom data organization, specify loci files explicitly:

```bash
credtools web my_results/ \
  --allmeta-loci my_results/allmeta/loci_list.txt \
  --popumeta-loci my_results/popumeta/loci_list.txt \
  --nometa-loci my_results/nometa/loci_list.txt
```

### Force Data Regeneration

When underlying results change, regenerate web data:

```bash
credtools web results/ --force-regenerate --threads 16
```

## Web Interface Features

### Main Dashboard

**Locus Browser**
- Searchable table of all analyzed loci
- Sortable by various metrics (max PIP, credible set size, etc.)
- Filterable by chromosome, ancestry, significance

**Summary Statistics**
- Overall analysis metrics and success rates
- Distribution of credible set sizes
- Cross-ancestry comparison summaries

### Locus Detail View

**Manhattan-style Plots**
- Posterior inclusion probabilities vs. genomic position
- Association p-values vs. genomic position  
- Linkage disequilibrium heatmaps

**Credible Set Visualization**
- Interactive credible set boundaries
- Variant annotation and details
- Downloadable high-resolution plots

**Multi-Ancestry Comparison**
- Side-by-side PIP comparisons
- Ancestry-specific vs. meta-analysis results
- Effect size and frequency comparisons

### Interactive Features

**Dynamic Filtering**
- Filter by PIP thresholds
- Select specific chromosomes or regions
- Focus on high-confidence variants

**Customizable Views**
- Adjust plot parameters in real-time
- Toggle data layers on/off
- Export custom figure configurations

**Data Export**
- Download plots as PNG/PDF
- Export filtered data tables
- Generate summary reports

## Expected Output

### Web Data Directory Structure

```
webdata/
├── all_loci_info.txt          # Master loci summary
├── locus_data/                # Individual locus files
│   ├── locus_1.json
│   ├── locus_2.json
│   └── ...
├── summary_stats.json         # Global summary statistics
├── plots/                     # Pre-generated plot data
│   ├── manhattan_data.json
│   └── summary_plots.json
└── metadata.json              # Processing metadata
```

### Web Server Output

When launched, the web server provides:
- **Local URL**: http://localhost:8080 (or specified port)
- **Network URL**: http://[host-ip]:[port] for remote access
- **Real-time logs**: Server activity and user interactions
- **Performance metrics**: Page load times and data processing

## Examples

### Example 1: Basic Web Launch

```bash
# Launch web interface on default port
credtools web finemap_results/
# Opens browser to http://localhost:8080
```

### Example 2: Multi-User Server Setup

```bash
# Launch server accessible from network
credtools web results/ \
  --host 0.0.0.0 \
  --port 8080 \
  --threads 20

# Accessible at http://server-ip:8080
```

### Example 3: Custom Data Organization

```bash
# Web interface for custom analysis structure
credtools web /path/to/analysis/ \
  --allmeta-loci /path/to/analysis/meta_all/loci.txt \
  --popumeta-loci /path/to/analysis/meta_pop/loci.txt \
  --webdata-dir /path/to/webdata \
  --force-regenerate
```

### Example 4: High-Performance Setup

```bash
# Optimized for large datasets
credtools web large_study_results/ \
  --threads 32 \
  --force-regenerate \
  --port 9000 \
  --host 0.0.0.0

# Pre-processes all data with maximum parallelization
```

### Example 5: Development and Debugging

```bash
# Debug mode with detailed logging
credtools web test_results/ \
  --debug \
  --port 5000 \
  --force-regenerate

# Useful for development and troubleshooting
```

## Navigation Guide

### Getting Started

1. **Launch the web interface** using the command above
2. **Open your browser** to the displayed URL
3. **Browse the locus table** to see all available results
4. **Click on a locus** to view detailed visualizations
5. **Use filters** to focus on specific subsets of data

### Key Interface Elements

**Locus Table**
- Sort by clicking column headers
- Filter using the search box or dropdown menus
- Click row to open detailed locus view

**Locus Plots**
- Manhattan plot shows PIPs across the locus
- LD heatmap displays correlation structure
- Credible sets highlighted with colored regions

**Control Panels**
- Adjust plot parameters with sliders and dropdowns
- Toggle ancestry-specific views
- Export plots and data

### Advanced Features

**Cross-Locus Analysis**
- Compare credible set sizes across loci
- Identify loci with ancestry-specific signals
- Generate summary statistics and reports

**Custom Visualizations**
- Combine multiple plot types
- Overlay additional annotation data
- Create publication-ready figures

## Troubleshooting

### Common Issues

**Web data not found**: Run with `--force-regenerate` to create web-optimized data files

**Server won't start**: Check if port is already in use, try different port with `--port`

**Slow loading**: Reduce number of threads or increase server memory allocation

**Browser compatibility**: Use modern browsers (Chrome, Firefox, Safari, Edge)

**Empty visualizations**: Verify fine-mapping results are in expected format and location

### Performance Optimization

**Large Datasets:**
- Use `--force-regenerate` with high thread count for initial setup
- Consider subsetting loci for faster exploration
- Use dedicated server for multi-user access

**Network Access:**
- Configure firewall rules for specified port
- Use `--host 0.0.0.0` for network accessibility
- Consider reverse proxy for production deployment

### Data Quality Issues

**Missing loci**: Check that loci files point to actual fine-mapping results
**Incomplete results**: Some loci may fail processing - check logs for details
**Inconsistent formats**: Ensure all fine-mapping outputs follow credtools format
**Large memory usage**: Monitor system resources during data processing

## Security Considerations

### Network Deployment

**Production Use:**
- Use reverse proxy (nginx, Apache) for secure deployment
- Enable HTTPS for sensitive data
- Implement access controls and authentication
- Monitor server logs for suspicious activity

**Data Privacy:**
- Ensure compliance with institutional data policies
- Consider data anonymization for shared access
- Use secure networks for sensitive genetic data
- Implement proper user access controls

## Tips for Success

1. **Pre-process data offline**: Use `--force-regenerate` before sharing with others
2. **Optimize for your use case**: Adjust thread count based on available resources
3. **Use bookmarks**: Save URLs for specific loci or views of interest
4. **Export early and often**: Download figures and tables as you explore
5. **Share responsibly**: Consider data sensitivity when providing network access
6. **Monitor performance**: Watch server logs and system resources during use
7. **Keep data updated**: Regenerate web data when underlying analyses change

## Integration with Workflow

The web interface is typically the final step for results exploration:

```bash
# Complete workflow ending with web visualization
credtools munge ancestry_files.json munged/
credtools chunk munged/ chunked/
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/
credtools pipeline prepared/final_loci_list.txt finemap_results/
credtools web finemap_results/ --port 8080

# Share URL with collaborators for interactive exploration
```