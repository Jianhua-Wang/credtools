# Quick Start Guide

This guide will get you running your first CREDTOOLS analysis in just a few minutes using the `credtools pipeline` command - the easiest way to perform end-to-end multi-ancestry fine-mapping.

## What is `credtools pipeline`?

The `credtools pipeline` command runs the complete CREDTOOLS workflow in a single command:

1. Quality control
2. Meta-analysis
3. Fine-mapping
4. Results aggregation

## Input Data Format

CREDTOOLS requires a tab-separated file describing your loci and studies. Here's the required format:

| Column | Description | Example |
|--------|-------------|---------|
| `chr` | Chromosome | `8` |
| `start` | Start position (bp) | `41242482` |
| `end` | End position (bp) | `42492482` |
| `popu` | Population/ancestry | `EUR`, `AFR`, `SAS`, `HIS` |
| `sample_size` | Sample size | `337465` |
| `cohort` | Cohort/study name | `UKBB`, `MVP` |
| `prefix` | File path prefix | `/path/to/data/EUR.UKBB.chr8_41242482_42492482` |
| `locus_id` | Locus identifier | `chr8_41242482_42492482` |

!!! info "File Structure"
    For each `prefix`, CREDTOOLS expects these files:
    
    - `{prefix}.sumstats` - Summary statistics
    - `{prefix}.ld` or `{prefix}.ld.npz` - LD matrix  
    - `{prefix}.ldmap` - LD matrix variant map

## Example Input File

```bash title="my_loci.txt"
chr	start	end	popu	sample_size	cohort	prefix	locus_id
8	41242482	42492482	AFR	89499	MVP	data/AFR.MVP.chr8_41242482_42492482	chr8_41242482_42492482
8	41242482	42492482	EUR	337465	MVP	data/EUR.MVP.chr8_41242482_42492482	chr8_41242482_42492482
8	41242482	42492482	EUR	442817	UKBB	data/EUR.UKBB.chr8_41242482_42492482	chr8_41242482_42492482
```

## Basic Usage

### Simple Cross-Ancestry Analysis

```bash
credtools pipeline my_loci.txt output_dir \
    --tool susie \    --threads 4
```

This command:

- Combines all studies across ancestries (`meta_all`)
- Uses multi-input strategy with MultiSuSiE
- Outputs results to `output_dir/`

### Population-Specific Analysis  

```bash
credtools pipeline my_loci.txt output_dir \
    --tool susie \    --threads 4 \
    --max-causal 5 \
    --credible-level 0.95
```

This command:

- Meta-analyzes within each ancestry separately (`meta_by_population`)
- Runs SuSiE on each population, then combines results (`post_hoc_combine`)

## Understanding the Output

After running `credtools pipeline`, you'll find these files in your output directory:

### Meta-Analysis Results
```
output_dir/
├── {locus_id}.{popu}.{cohort}.sumstat    # Meta-analyzed summary stats
├── {locus_id}.{popu}.{cohort}.ld.npz     # Meta-analyzed LD matrix
└── {locus_id}.{popu}.{cohort}.ldmap      # LD variant mapping
```

### Quality Control Reports
```
output_dir/
├── s_estimate.txt        # Inconsistency parameter estimates
├── kriging_rss.txt       # Allele switch detection
├── maf_comparison.txt    # MAF consistency across studies
├── cochran_q.txt         # Heterogeneity testing
└── ld_structure.txt      # LD matrix eigenanalysis
```

### Fine-Mapping Results
```
output_dir/
├── pips.txt             # Posterior inclusion probabilities
└── creds.json           # Credible sets information
```

## Interpreting Results

### Posterior Inclusion Probabilities (PIPs)

The `pips.txt` file contains PIPs for each variant:

```bash title="pips.txt"
8-41234567-A-G	0.0234
8-41235678-C-T	0.8765
8-41236789-G-A	0.0456
```

- Values range from 0 to 1
- Higher values indicate stronger evidence for causality
- Typically, variants with PIP > 0.1 are considered noteworthy

### Credible Sets

The `creds.json` file contains credible sets - groups of variants that collectively have high probability of containing the causal variant:

```json title="creds.json"
{
  "credible_sets": {
    "cs1": {
      "variants": ["8-41235678-C-T", "8-41235680-A-G"],
      "coverage": 0.95,
      "total_pip": 0.96
    }
  }
}
```

## Common Options

### Meta-Analysis Methods

```bash
# Combine all studies regardless of ancestry
--meta-method meta_all

# Combine studies within each ancestry separately  
--meta-method meta_by_population

# Keep all studies separate (no meta-analysis)
--meta-method no_meta
```

### Fine-Mapping Tools

```bash
# General purpose, robust
--tool susie

# Multi-ancestry designed tools
--tool multisusie
--tool susiex

# Bayesian model averaging
--tool finemap

# Simple Bayes factors  
--tool abf
```

### Quality Control

```bash
# Skip QC (faster but not recommended)
--skip-qc

# Include QC (default, recommended)
# No flag needed - QC runs by default
```

## Troubleshooting

!!! warning "Common Issues"
    
    **File not found errors**
    : Check that your file paths in the input table are correct
    : Ensure summary statistics and LD files exist for each prefix
    
    **Memory errors**
    : Large LD matrices can consume significant memory
    : Consider analyzing loci one at a time for very large regions
    
    **Tool-specific errors**
    : Some tools have specific requirements (see tool documentation)
    : Try SuSiE first as it's the most robust default option

!!! tip "Performance Tips"
    
    - Start with smaller regions to test your setup
    - Use `--tool susie` for initial exploration (fastest, most reliable)
    - Save QC results to identify problematic studies before fine-mapping

## Next Steps

Once you've run your first analysis:

- **[Single-Input Fine-Mapping](single-input.md)** - Learn about analyzing individual studies
- **[Multi-Input Fine-Mapping](multi-input.md)** - Deep dive into multi-ancestry analysis  
- **[Advanced Topics](advanced.md)** - Customize parameters and understand tool options

## Example with Real Data

Using the included example data:

```bash
# Navigate to example data directory
cd exampledata/

# Run pipeline on example locus
credtools pipeline test_loci.txt results/ \
    --tool susie \    --threads 4
```

This will analyze the multi-ancestry example data and produce results in the `results/` directory. 