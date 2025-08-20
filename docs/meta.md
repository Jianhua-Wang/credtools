# Meta-Analysis

The `credtools meta` command performs meta-analysis of summary statistics and LD matrices across multiple ancestries or studies. This step combines evidence from different populations to improve fine-mapping resolution and power.

## Overview

Meta-analysis in credtools integrates genetic evidence across populations while accounting for different LD structures and effect sizes. The meta-analysis process:

- **Combines summary statistics** using inverse variance weighting or other methods
- **Merges LD matrices** appropriately across ancestries 
- **Handles population-specific effects** and heterogeneity
- **Creates unified datasets** for downstream fine-mapping
- **Supports flexible strategies** for different study designs
- **Preserves ancestry-specific information** when needed

## When to Use

Use `credtools meta` when you have:

- Fine-mapping data from multiple ancestries or populations
- Studies with overlapping genetic signals that could benefit from combined analysis
- Need to increase statistical power through meta-analysis
- Want to identify ancestry-specific vs shared causal variants
- Prepared locus files from multiple populations ready for integration

## Basic Usage

### Standard Meta-Analysis

```bash
credtools meta prepared_loci_list.txt meta_output/ --meta-method meta_all
```

### Population-Specific Meta-Analysis

```bash
credtools meta prepared_loci_list.txt meta_output/ --meta-method meta_by_population
```

### Skip Meta-Analysis

```bash
credtools meta prepared_loci_list.txt meta_output/ --meta-method no_meta
```

## Command Options

| Option | Description | Default |
|--------|-------------|---------|
| `--meta-method` / `-m` | Meta-analysis method | meta_all |
| `--threads` / `-t` | Number of parallel threads | 1 |

## Meta-Analysis Methods

### meta_all
Combines all ancestries into a single meta-analyzed dataset for each locus. This approach:
- Maximizes statistical power by pooling all available data
- Assumes similar effect sizes across populations
- Creates one unified summary statistic per locus
- Best for traits with consistent effects across ancestries

### meta_by_population  
Performs meta-analysis within ancestry groups while keeping populations separate. This approach:
- Maintains population-specific effect estimates
- Allows for heterogeneity between ancestries
- Creates separate datasets for each major ancestry group
- Suitable when effect sizes differ significantly between populations

### no_meta
Processes each ancestry/study separately without meta-analysis. This approach:
- Preserves all population-specific information
- Enables ancestry-specific fine-mapping
- Useful for identifying population-specific causal variants
- Required when populations have very different LD structures

## Expected Output

The meta-analysis process creates:

- **Meta-analyzed summary statistics**: Combined effect estimates for each method
- **Merged LD matrices**: Appropriately weighted correlation matrices
- **Loci information files**: Updated metadata for downstream analysis
- **Processing logs**: Detailed information about meta-analysis steps

### Output Structure

```
meta_output/
├── meta_all/           # Results from meta_all method
│   ├── locus_1/
│   │   ├── sumstats.txt
│   │   ├── ldmatrix.txt
│   │   └── locus_info.json
│   └── locus_2/
├── meta_by_population/ # Results from meta_by_population method
│   ├── EUR/
│   ├── ASN/
│   └── AFR/
└── no_meta/           # Results from no_meta method
    ├── study1/
    ├── study2/
    └── study3/
```

## Examples

### Example 1: Standard Multi-Ancestry Meta-Analysis

```bash
# Perform meta-analysis across all ancestries
credtools meta prepared/final_loci_list.txt meta_results/ \
  --meta-method meta_all \
  --threads 4

# Output: Combined datasets maximizing power across populations
```

### Example 2: Population-Specific Analysis

```bash
# Keep ancestry groups separate but meta-analyze within groups
credtools meta prepared/final_loci_list.txt meta_results/ \
  --meta-method meta_by_population \
  --threads 8

# Useful when you have multiple studies per ancestry
```

### Example 3: No Meta-Analysis (Individual Studies)

```bash
# Process each study/ancestry independently
credtools meta prepared/final_loci_list.txt meta_results/ \
  --meta-method no_meta \
  --threads 2

# Preserves all study-specific information
```

### Example 4: High-Performance Processing

```bash
# Use maximum threads for large datasets
credtools meta prepared/final_loci_list.txt meta_results/ \
  --meta-method meta_all \
  --threads 16

# Recommended for genome-wide studies with many loci
```

## Integration with Workflow

Meta-analysis fits into the credtools workflow after preparation:

```bash
# 1. Prepare LD matrices and inputs
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/

# 2. Perform meta-analysis
credtools meta prepared/final_loci_list.txt meta/ --meta-method meta_all

# 3. Run quality control (optional)
credtools qc meta/meta_all/loci_list.txt qc_results/

# 4. Perform fine-mapping
credtools finemap meta/meta_all/loci_list.txt finemap_results/
```

## Choosing the Right Method

### Use meta_all when:
- Effect sizes are consistent across populations
- You want maximum statistical power
- Trait biology suggests similar mechanisms across ancestries
- Sample sizes are relatively balanced

### Use meta_by_population when:
- You have multiple studies per ancestry
- Some heterogeneity in effect sizes is expected
- You want to preserve some population structure
- Studies within ancestry are more homogeneous than across ancestries

### Use no_meta when:
- Effect sizes differ substantially between populations
- You want to identify population-specific signals
- Studies have very different designs or phenotype definitions
- You plan to compare results across populations

## Troubleshooting

### Common Issues

**Memory issues with large datasets**: Reduce the number of threads or process subsets of loci separately.

**Inconsistent ancestry labels**: Ensure ancestry identifiers match between summary statistics and LD matrices.

**Missing files**: Check that all required input files from the prepare step are present and accessible.

**Convergence issues**: Some loci may fail meta-analysis due to data quality issues - these will be logged and skipped.

### Quality Checks

The meta command automatically performs several quality checks:
- Validates input file formats and required columns
- Checks for matching variants between summary statistics and LD matrices
- Identifies and handles allele mismatches
- Reports summary statistics for each locus processed

### Performance Optimization

1. **Use appropriate thread counts**: Generally use 1 thread per CPU core, but reduce for memory-intensive analyses
2. **Process by chromosome**: For very large studies, consider splitting by chromosome first
3. **Monitor memory usage**: Large LD matrices can consume significant RAM
4. **Use fast storage**: Place output directory on fast SSD storage when possible

## Tips for Success

1. **Validate inputs first**: Always run the prepare step successfully before meta-analysis
2. **Check ancestry matching**: Ensure genetic ancestry matches between summary stats and reference panels
3. **Document your choices**: Record which meta-analysis method you used and why
4. **Save intermediate results**: Keep meta-analysis outputs for reproducibility
5. **Plan for heterogeneity**: Consider population-specific analysis if you observe large effect differences
6. **Use consistent identifiers**: Maintain consistent locus and ancestry naming throughout your workflow