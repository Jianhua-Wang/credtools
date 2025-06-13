# Web Visualization

CREDTOOLS provides an interactive web interface for exploring fine-mapping results. This tutorial covers how to install, configure, and use the web visualization features.

## Prerequisites

1. **Completed fine-mapping analysis** with CREDTOOLS
2. **Web dependencies** installed
3. **Results directory** with fine-mapping output

## Installation

The web visualization requires additional dependencies that are not included in the base CREDTOOLS installation:

```bash
# Install CREDTOOLS with web dependencies
pip install credtools[web]
```

## Basic Usage

Get help on web visualization options:

```bash
credtools web --help
```

### Launch from Results Directory

```bash
# Navigate to your CREDTOOLS results directory
cd /path/to/your/credtools/results

# Start the web interface
credtools web
```

### Launch with Custom Port

```bash
credtools web /path/to/credtools/results --port 8080
```

## Data Structure

Your CREDTOOLS results should follow this structure:

```
credtools_results/
├── allmeta_loci.txt
├── popumeta_loci.txt
├── susie/
│   ├── susie_allmeta/
│   └── susie_popumeta/
└── finemap/
    ├── finemap_allmeta/
    └── finemap_popumeta/
```

## Advanced Usage

### Launch from Specific Directory

```bash
cd /path/to/credtools/results
credtools web
```

### Force Data Regeneration

```bash
credtools web /path/to/data \
    --allmeta-loci data/allmeta_loci.txt \
    --popumeta-loci data/popumeta_loci.txt \
    --force-regenerate
```

### Debug Mode

```bash
credtools web /path/to/data --debug --port 8081
```

## Python API

When you run `credtools web`, it automatically:

1. Processes your fine-mapping results
2. Generates web visualization data
3. Launches a web server

You can also use the Python API directly:

```python
from credtools.web.export import export_for_web
from credtools.web.app import run_app

# Process data for web visualization
export_for_web(
    data_base_dir="/path/to/credtools/results",
    allmeta_loci="allmeta_loci.txt",
    popumeta_loci="popumeta_loci.txt"
)

# Launch web server
run_app(port=8080)
```

## Customization

### Custom Web Data Directory

```bash
credtools web /path/to/data --webdata-dir /custom/path/webdata
```

### Multiple Loci Files

```bash
credtools web /path/to/data \
    --allmeta-loci data/allmeta_loci.txt \
    --popumeta-loci data/popumeta_loci.txt \
    --susie-loci data/susie_loci.txt
```

### Performance Tuning

```bash
# Use more threads for data processing
credtools web /path/to/data --threads 30

# Use fast storage for web data
credtools web /path/to/data --webdata-dir /fast/ssd/webdata
```

## Troubleshooting

### Missing Dependencies

If you see import errors:

```bash
pip install credtools[web]
```

### Custom Loci Files

If your loci files are in a different location:

```bash
credtools web /path/to/data --allmeta-loci /path/to/loci.txt
```

### Port Conflicts

Use a different port:

```bash
credtools web /path/to/data --port 8081
```

### Slow Processing

Reduce thread count:

```bash
credtools web /path/to/data --threads 5
```

### Debug Mode

Enable debug logging:

```bash
credtools web /path/to/data --debug
```

## Logging

Check CREDTOOLS logs for processing issues:

```bash
tail -f credtools.log
```

## Example Workflow

1. Run fine-mapping pipeline:

```bash
credtools pipeline input_loci.txt results/
```

2. Launch web interface:

```bash
credtools web results/ --port 8080
```

## Multiple Datasets

For multiple datasets:

```bash
for i in {1..3}; do
    dataset="dataset${i}"
    credtools web results/${dataset} --webdata-dir webdata/${dataset} --port 808${i}
done
``` 