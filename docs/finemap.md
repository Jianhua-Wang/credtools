# Fine-Mapping

The `credtools finemap` command performs statistical fine-mapping analysis to identify causal genetic variants within loci. This is the core analysis step that implements various fine-mapping algorithms and strategies to prioritize likely causal variants.

## Overview

Fine-mapping in credtools uses sophisticated statistical methods to narrow down large association signals to smaller sets of likely causal variants. The fine-mapping process:

- **Implements multiple algorithms** including SuSiE, ABF, FINEMAP, RSparsePro, and others
- **Supports three strategies** for handling multi-ancestry data
- **Calculates posterior inclusion probabilities (PIPs)** for each variant
- **Generates credible sets** of variants likely to contain causal variants
- **Handles multi-signal loci** with conditional analysis approaches
- **Provides comprehensive output** including detailed statistics and visualizations

## When to Use

Use `credtools finemap` when you have:

- Quality-controlled summary statistics and LD matrices ready for analysis
- Defined loci from genome-wide association studies
- Need to identify the most likely causal variants within association signals
- Want to compare results across different fine-mapping methods
- Multi-ancestry data requiring sophisticated modeling approaches

## Supported Tools

### SuSiE (Default)
**Sum of Single Effects model**
- Handles multiple causal variants naturally
- Provides interpretable credible sets
- Fast and robust for most applications
- Good for loci with unclear causal architecture

### ABF
**Approximate Bayes Factor**
- Simple and fast single-variant method
- Good baseline for comparison
- Assumes single causal variant per locus
- Suitable for well-defined single signals

### ABF+COJO
**ABF with Conditional Analysis**
- Combines COJO stepwise selection with ABF
- Handles multi-signal loci effectively
- Uses conditional analysis for each independent signal
- Good compromise between simplicity and sophistication

### FINEMAP
**Shotgun Stochastic Search**
- Considers multiple causal variant configurations
- Provides model-averaged results
- Computationally intensive but thorough
- Excellent for complex multi-signal loci

### RSparsePro
**Robust Sparse Regression**
- Handles model uncertainty explicitly
- Good for noisy or heterogeneous data
- Robust to model misspecification
- Suitable for challenging datasets

### multiSuSiE
**Multi-ancestry SuSiE**
- Designed specifically for multi-population data
- Models population-specific and shared effects
- Leverages information across ancestries
- Best for multi-ancestry fine-mapping

### SuSiEx
**Cross-population SuSiE**
- Advanced multi-population method
- Accounts for different LD structures
- Models population-specific effect sizes
- State-of-the-art for multi-ancestry studies

## Fine-Mapping Strategies

### single_input
Analyzes each ancestry separately, then combines results post-hoc
- Maintains ancestry-specific effects
- Simple and interpretable
- Good when populations differ substantially
- Fastest approach for multi-ancestry data

### multi_input  
Jointly analyzes all ancestries simultaneously
- Uses methods designed for multi-population data
- Shares information across populations
- More powerful when effects are similar
- Requires compatible fine-mapping tools

### post_hoc_combine
Runs single-ancestry analysis then combines evidence
- Flexible combination strategies
- Can weight ancestries differently
- Good for heterogeneous effect sizes
- Allows for sensitivity analyses

## Basic Usage

### Single Ancestry with SuSiE

```bash
credtools finemap prepared/loci_list.txt results/ --tool susie
```

### Multi-Ancestry with multiSuSiE

```bash
credtools finemap prepared/loci_list.txt results/ \
  --strategy multi_input --tool multisusie
```

### Complex Multi-Signal Analysis

```bash
credtools finemap qc/passed_loci.txt results/ \
  --tool susie --max-causal 5 --coverage 0.95
```

## Command Options

### General Options

| Option | Description | Default |
|--------|-------------|---------|
| `--strategy` / `-s` | Fine-mapping strategy | single_input |
| `--tool` / `-t` | Fine-mapping tool | susie |
| `--max-causal` / `-c` | Maximum causal variants | 1 |
| `--coverage` / `-cv` | Credible set coverage | 0.95 |
| `--combine-cred` / `-cc` | Method to combine credible sets | union |
| `--combine-pip` / `-cp` | Method to combine PIPs | max |
| `--jaccard-threshold` / `-j` | Jaccard threshold for combination | 0.1 |

### COJO Options (for ABF+COJO)

| Option | Description | Default |
|--------|-------------|---------|
| `--set-L-by-cojo` / `-sl` | Set max causal by COJO results | True |
| `--p-cutoff` / `-pc` | P-value cutoff for COJO | 5e-8 |
| `--collinear-cutoff` / `-cc` | Collinearity threshold | 0.9 |
| `--window-size` / `-ws` | COJO window size | 10000000 |
| `--maf-cutoff` / `-mc` | MAF cutoff | 0.01 |
| `--diff-freq-cutoff` / `-dfc` | Frequency difference cutoff | 0.2 |

### Tool-Specific Options

#### SuSiE Parameters
| Option | Description | Default |
|--------|-------------|---------|
| `--max-iter` / `-i` | Maximum iterations | 100 |
| `--estimate-residual-variance` / `-er` | Estimate residual variance | False |
| `--min-abs-corr` / `-mc` | Minimum absolute correlation | 0.5 |
| `--convergence-tol` / `-ct` | Convergence tolerance | 1e-3 |

#### ABF Parameters  
| Option | Description | Default |
|--------|-------------|---------|
| `--var-prior` / `-vp` | Variance prior | 0.2 |

#### FINEMAP Parameters
| Option | Description | Default |
|--------|-------------|---------|
| `--n-iter` / `-ni` | Number of iterations | 100000 |
| `--n-threads` / `-nt` | Number of threads | 1 |

#### RSparsePro Parameters
| Option | Description | Default |
|--------|-------------|---------|
| `--eps` / `-e` | Convergence criterion | 1e-5 |
| `--ubound` / `-ub` | Upper bound for convergence | 100000 |
| `--cthres` / `-ct` | Coverage threshold | 0.7 |

## Expected Output

Fine-mapping generates comprehensive results for each locus:

### Standard Output Files

```
results/
├── locus_1/
│   ├── pips.txt              # Posterior inclusion probabilities
│   ├── creds.json            # Credible sets and metadata
│   ├── tool_specific/        # Tool-specific outputs
│   │   ├── susie_results.rds
│   │   └── summary.txt
│   └── plots/                # Visualization files
│       ├── pip_plot.png
│       └── credset_plot.png
└── locus_2/
    └── ...
```

### Key Output Files

**pips.txt**: Tab-separated file with variants and their posterior inclusion probabilities
- Columns: SNPID, CHR, BP, EA, NEA, BETA, SE, P, PIP
- Sorted by decreasing PIP values
- Includes all variants analyzed

**creds.json**: JSON file containing:
- Credible sets (groups of variants likely to contain causal variant)
- Coverage probabilities for each set
- Method metadata and parameters
- Quality metrics and convergence information

## Examples

### Example 1: Standard Single-Ancestry Analysis

```bash
# Run SuSiE on single ancestry with default parameters
credtools finemap prepared/EUR_loci.txt results_eur/ \
  --tool susie --max-causal 3

# Output: PIPs and credible sets for European ancestry
```

### Example 2: Multi-Ancestry Joint Analysis

```bash
# Use multiSuSiE for joint analysis across ancestries
credtools finemap prepared/multi_loci.txt results_multi/ \
  --strategy multi_input --tool multisusie \
  --max-causal 5 --coverage 0.99

# Leverages information across all populations simultaneously
```

### Example 3: Comparative Analysis with Multiple Tools

```bash
# Compare results across different methods
for tool in susie abf finemap; do
  credtools finemap qc/loci.txt results_${tool}/ \
    --tool $tool --max-causal 2
done

# Allows method comparison and sensitivity analysis
```

### Example 4: Complex Multi-Signal Locus

```bash
# Use sophisticated method for complex locus
credtools finemap complex_locus.txt results_complex/ \
  --tool susie --max-causal 10 \
  --estimate-residual-variance \
  --coverage 0.95

# Handles loci with many independent signals
```

### Example 5: High-Confidence Fine-Mapping

```bash
# Use strict parameters for high-confidence results
credtools finemap validated_loci.txt results_strict/ \
  --tool multisusie --strategy multi_input \
  --max-causal 3 --coverage 0.99 \
  --jaccard-threshold 0.2

# More stringent credible sets and combination
```

## Integration with Workflow

Fine-mapping is typically the final analysis step:

```bash
# Complete workflow ending with fine-mapping
credtools munge ancestry_files.json munged/
credtools chunk munged/ chunked/
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/
credtools meta prepared/final_loci_list.txt meta/
credtools qc meta/meta_all/loci_list.txt qc/
credtools finemap qc/passed_loci.txt finemap_results/ \
  --tool susie --max-causal 3

# Launch web interface to explore results
credtools web finemap_results/
```

## Choosing the Right Tool

### Use SuSiE when:
- General-purpose fine-mapping with good performance
- Loci may have multiple causal variants
- Want interpretable credible sets
- Computational efficiency is important

### Use ABF when:
- Simple baseline analysis
- Confident loci have single causal variants
- Very fast results needed
- Comparing with literature using ABF

### Use ABF+COJO when:
- Want to handle multi-signal loci with simple method
- Need conditional analysis approach
- Balancing simplicity with multi-signal capability

### Use FINEMAP when:
- Complex loci with many potential causal variants
- Want comprehensive model averaging
- Computational resources are abundant
- Need most sophisticated single-ancestry method

### Use multiSuSiE when:
- Multi-ancestry data with shared effects
- Want to leverage cross-population information
- Effects are reasonably similar across populations

### Use SuSiEx when:
- Multi-ancestry data with population-specific effects
- Different LD structures across populations
- Want state-of-the-art multi-population method

## Interpreting Results

### Posterior Inclusion Probabilities (PIPs)
- **PIP > 0.95**: Very high confidence causal variant
- **PIP 0.5-0.95**: Moderate confidence, worth following up
- **PIP 0.1-0.5**: Low to moderate evidence
- **PIP < 0.1**: Little evidence for causality

### Credible Sets
- **95% credible set**: Contains causal variant with 95% probability
- **Multiple sets per locus**: Indicates multiple independent signals
- **Large sets (>100 variants)**: Suggests limited fine-mapping resolution
- **Small sets (<10 variants)**: Good fine-mapping resolution

### Quality Metrics
- **Convergence**: Algorithm reached stable solution
- **Model fit**: How well model explains the data
- **Coverage**: Credible set size vs. probability coverage
- **Consistency**: Results similar across methods/ancestries

## Troubleshooting

### Common Issues

**Convergence failures**: Increase max iterations or adjust convergence tolerance
**Memory errors**: Reduce max causal variants or use single-ancestry strategy
**No credible sets**: Lower coverage threshold or increase max causal variants
**Very large credible sets**: May indicate poor LD coverage or complex architecture

### Performance Optimization

1. **Choose appropriate max causal**: Start with 3-5, increase only if needed
2. **Use parallel processing**: Some tools support multithreading
3. **Optimize memory usage**: Consider single-ancestry strategy for very large studies
4. **Monitor convergence**: Check logs for convergence issues

## Tips for Success

1. **Start with SuSiE**: Good default choice for most applications
2. **Validate with multiple methods**: Compare results across tools when possible
3. **Consider ancestry differences**: Use appropriate multi-ancestry strategy
4. **Check convergence**: Ensure algorithms converged properly
5. **Interpret probabilistically**: Focus on PIPs rather than hard cutoffs
6. **Follow up high-PIP variants**: Validate top variants with additional analyses