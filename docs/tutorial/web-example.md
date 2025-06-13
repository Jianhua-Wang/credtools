# Web Visualization Example

This example demonstrates a complete CREDTOOLS workflow including web visualization.

## Sample Dataset

For this example, we'll use a hypothetical fine-mapping analysis with three loci across multiple populations.

### Directory Structure

```
example_analysis/
├── input/
│   └── loci_info.txt
├── results/
│   ├── data/
│   │   └── real/
│   │       ├── credset/
│   │       └── qc/
│   └── webdata/  # Generated automatically
└── loci_files/
    ├── allmeta_loci.txt
    ├── popumeta_loci.txt
    └── nometa_loci.txt
```

## Step-by-Step Workflow

### 1. Run Fine-mapping Pipeline

First, run the complete fine-mapping pipeline:

```bash
# Navigate to your analysis directory
cd example_analysis/

# Run the complete pipeline
credtools pipeline input/loci_info.txt results/ \
  --tool susie \
  --strategy multi_input \
  --meta-method meta_all \
  --threads 10
```

This generates:
- Fine-mapping results in `results/data/real/credset/`
- Quality control metrics in `results/data/real/qc/`
- Meta-analysis results if applicable

### 2. Prepare Loci Information Files

Create loci information files for web visualization:

**allmeta_loci.txt:**
```
locus_id	chr	start	end	prefix	popu	cohort	sample_size
locus1	1	1000000	2000000	/path/to/locus1_data	EUR	cohort1	50000
locus1	1	1000000	2000000	/path/to/locus1_data	ASN	cohort2	30000
locus2	2	5000000	6000000	/path/to/locus2_data	EUR	cohort1	50000
locus3	3	8000000	9000000	/path/to/locus3_data	AFR	cohort3	25000
```

### 3. Launch Web Interface

#### Basic Launch

```bash
# Simple launch from results directory
cd results/
credtools web
```

#### Custom Configuration

```bash
# Launch with specific settings
credtools web results/ \
  --allmeta-loci loci_files/allmeta_loci.txt \
  --popumeta-loci loci_files/popumeta_loci.txt \
  --nometa-loci loci_files/nometa_loci.txt \
  --port 8080 \
  --threads 15
```

### 4. Explore Results

Once the web interface starts, open your browser to `http://localhost:8080`.

#### Home Page Features

1. **Filter by Meta-analysis Method:**
   - Select "allmeta" to see all-ancestry meta-analysis results
   - Choose "popumeta" for population-specific results
   - Pick "nometa" for individual cohort results

2. **Filter by Fine-mapping Tool:**
   - Compare results across SuSiE, FINEMAP, etc.
   - Each tool may show different credible sets

3. **View Summary Statistics:**
   - Number of credible sets per locus
   - Total credible set sizes
   - SNPs with high posterior probabilities

#### Locus-Specific Views

Click on any locus ID to see detailed results:

1. **Association Plots:**
   - Manhattan plot with LD coloring
   - Fine-mapping posterior probabilities

2. **Quality Control:**
   - Lambda inflation values
   - DENTIST-S statistics
   - MAF correlation metrics

3. **Credible Sets:**
   - Highlighted credible variants
   - Posterior inclusion probabilities
   - Cross-tool comparisons

## Advanced Usage

### Programmatic Access

You can also process data and launch the web interface programmatically:

```python
from credtools.web.export import export_for_web
from credtools.web.app import run_app

# Process data for web visualization
export_for_web(
    data_base_dir="results/",
    webdata_dir="webdata/",
    allmeta_loci_file="loci_files/allmeta_loci.txt",
    popumeta_loci_file="loci_files/popumeta_loci.txt",
    nometa_loci_file="loci_files/nometa_loci.txt",
    threads=10
)

# Launch web application
run_app(
    webdata_dir="webdata/",
    port=8080,
    debug=True
)
```

### Batch Processing

For multiple datasets:

```bash
#!/bin/bash
# Process multiple result directories

datasets=("dataset1" "dataset2" "dataset3")
port=8080

for dataset in "${datasets[@]}"; do
    echo "Processing $dataset..."
    
    # Run fine-mapping if needed
    credtools pipeline input/${dataset}_loci.txt results/${dataset}/
    
    # Launch web interface on different ports
    credtools web results/${dataset}/ \
      --port $((port++)) \
      --webdata-dir webdata/${dataset} &
done

echo "All web interfaces launched. Check ports 8080-8082"
```

### Custom Styling

The web interface uses Bootstrap themes. You can customize the appearance by modifying the Dash app configuration in your scripts.

## Troubleshooting

### Common Issues

1. **No data appears:** Check that loci files have the correct format and paths
2. **Slow loading:** Reduce the number of threads or use SSD storage
3. **Port conflicts:** Use a different port with `--port` option

### Getting Help

```bash
# Get command help
credtools web --help

# Check version
credtools --version

# Enable debug mode
credtools web results/ --debug
```

## Next Steps

- Try filtering by different meta-analysis methods
- Explore individual loci in detail
- Export plots for presentations
- Integrate web visualization into your analysis pipeline
- Share results with collaborators via the web interface

For more information:
- [Full Web Tutorial](web-visualization.md)
- [Advanced Usage](advanced.md)
- [API Documentation](../api.md) 