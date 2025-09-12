# Single-Input Fine-Mapping

Single-input fine-mapping is the traditional approach where you analyze one cohort or ancestry at a time. This strategy is ideal when you have well-powered individual studies or want to understand ancestry-specific genetic architecture.

## When to Use Single-Input Strategy

!!! tip "Single-Input Use Cases"
    
    **Well-powered individual studies**
    : When each study has sufficient power for fine-mapping alone
    
    **Ancestry-specific analysis**
    : When you want to understand population-specific causal variants
    
    **Exploratory analysis**
    : When investigating individual study characteristics before meta-analysis
    
    **Heterogeneous effects**
    : When effect sizes vary significantly across populations

## Step-by-Step Workflow

### 1. Prepare Single-Study Input

For single-input analysis, your loci file should contain only one row per locus:

```bash title="single_study_loci.txt"
chr	start	end	popu	sample_size	cohort	prefix	locus_id
8	41242482	42492482	EUR	442817	UKBB	data/EUR.UKBB.chr8_41242482_42492482	chr8_41242482_42492482
```

### 2. Quality Control

Quality control for single studies focuses on internal consistency:

```bash
credtools qc single_study_loci.txt qc_output/
```

#### QC Metrics for Single Studies

**Inconsistency Parameter (s)**
: Measures consistency between z-scores and LD matrix
: Values > 0.2 suggest potential issues

**Kriging RSS**
: Detects potential allele switches or data quality issues
: Flags variants with unexpected z-scores given LD structure

**LD Structure Analysis**  
: Eigenvalue decomposition of LD matrix
: Identifies problematic LD patterns

**MAF Comparison**
: Compares allele frequencies between summary stats and LD reference
: Large differences suggest potential strand issues

!!! warning "Single-Study QC Limitations"
    
    Single-study QC cannot detect:
    - Cross-study heterogeneity  
    - Population-specific biases
    - Cohort-specific technical artifacts
    
    Consider multi-study QC when possible.

### 3. Fine-Mapping

Run fine-mapping on your single study:

```bash
credtools finemap single_study_loci.txt output/ \    --tool susie \
    --max-causal 5 \
    --coverage 0.95
```

## Supported Tools for Single-Input

### SuSiE (Recommended)

**Best for**: General purpose, robust across different scenarios

```bash
credtools finemap input.txt output/ \    --tool susie \
    --max-causal 10 \
    --max-iter 100 \
    --estimate-residual-variance \
    --convergence-tol 1e-3
```

**Key Parameters:**

- `--max-causal`: Maximum number of causal variants (default: 1)
- `--max-iter`: Maximum iterations (default: 100)  
- `--estimate-residual-variance`: Estimate phenotype variance (default: False)
- `--convergence-tol`: Convergence tolerance (default: 1e-3)

### FINEMAP

**Best for**: Bayesian model averaging, comprehensive uncertainty quantification

```bash
credtools finemap input.txt output/ \    --tool finemap \
    --max-causal 5 \
    --n-iter 100000 \
    --n-threads 4
```

**Key Parameters:**

- `--n-iter`: Number of MCMC iterations (default: 100000)
- `--n-threads`: Number of parallel threads (default: 1)

### ABF (Approximate Bayes Factors)

**Best for**: Fast, simple analysis with minimal assumptions

```bash
credtools finemap input.txt output/ \    --tool abf \
    --var-prior 0.2
```

**Key Parameters:**

- `--var-prior`: Prior variance (0.15 for quantitative, 0.2 for binary traits)

### CARMA

**Best for**: Model uncertainty quantification, outlier detection  

```bash
credtools finemap input.txt output/ \    --tool carma \
    --max-causal 10 \
    --effect-size-prior "Spike-slab" \
    --y-var 1.0
```

**Key Parameters:**

- `--effect-size-prior`: "Spike-slab" or "Cauchy" (default: "Spike-slab")
- `--y-var`: Phenotype variance (default: 1.0)
- `--bf-threshold`: Bayes factor threshold (default: 10.0)

### RSparsePro  

**Best for**: Sparse regression approach, computational efficiency

```bash
credtools finemap input.txt output/ \    --tool rsparsepro \
    --max-causal 5 \
    --eps 1e-5 \
    --cthres 0.7
```

**Key Parameters:**

- `--eps`: Convergence criterion (default: 1e-5)
- `--cthres`: Coverage threshold (default: 0.7)
- `--minldthres`: Minimum LD within effect groups (default: 0.7)

## Comparing Single-Input Tools

| Tool | Speed | Memory | Uncertainty | Best Use Case |
|------|-------|--------|-------------|---------------|
| **SuSiE** | Fast | Low | Good | General purpose |
| **FINEMAP** | Moderate | Moderate | Excellent | Comprehensive analysis |
| **ABF** | Very Fast | Very Low | Basic | Quick screening |
| **CARMA** | Slow | High | Excellent | Research applications |
| **RSparsePro** | Fast | Low | Good | Large regions |

## Automatic Parameter Setting

CREDTOOLS can automatically determine the maximum number of causal variants using COJO:

```bash
credtools finemap input.txt output/ \    --tool susie \
    --set-L-by-cojo \
    --p-cutoff 5e-8 \
    --collinear-cutoff 0.9
```

**COJO Parameters:**

- `--p-cutoff`: P-value threshold for conditioning (default: 5e-8)
- `--collinear-cutoff`: Collinearity threshold (default: 0.9)
- `--window-size`: Window for conditional analysis (default: 10Mb)
- `--maf-cutoff`: MAF cutoff (default: 0.01)

!!! tip "COJO Integration"
    
    COJO (Conditional and Joint analysis) estimates the number of independent signals in a region. CREDTOOLS uses this to set `--max-causal` automatically, which often works better than arbitrary values.

## Output Files

### Posterior Inclusion Probabilities
```bash title="output/pips.txt"
# SNPID	PIP
8-41234567-A-G	0.0234
8-41235678-C-T	0.8765
8-41236789-G-A	0.0456
```

### Credible Sets
```json title="output/creds.json"
{
  "credible_sets": {
    "cs1": {
      "variants": ["8-41235678-C-T", "8-41235680-A-G"],
      "coverage": 0.95,
      "total_pip": 0.96,
      "min_abs_corr": 0.34
    }
  },
  "tool": "susie",
    "parameters": {...}
}
```

## Common Issues and Solutions

!!! warning "Troubleshooting Single-Input Analysis"
    
    **No credible sets found**
    : Try increasing `--max-causal` 
    : Check if region has sufficient signal (`--p-cutoff`)
    : Verify LD matrix quality
    
    **Very large credible sets**
    : May indicate weak signal or LD issues
    : Try more stringent `--coverage` (e.g., 0.99)
    : Consider region subdivision
    
    **Tool convergence issues**
    : Increase `--max-iter` for iterative methods
    : Try different `--convergence-tol` values
    : Switch to more robust tool (SuSiE)
    
    **Memory issues with large regions**
    : Use RSparsePro for efficiency
    : Consider region subdivision
    : Reduce LD matrix precision if possible

## Example: Complete Single-Study Analysis

```bash
# 1. Quality control
credtools qc my_study.txt qc_results/

# 2. Review QC metrics (check s_estimate.txt, kriging_rss.txt)

# 3. Fine-mapping with automatic L setting
credtools finemap my_study.txt results/ \    --tool susie \
    --set-L-by-cojo \
    --coverage 0.95 \
    --max-iter 100

# 4. Review results (pips.txt, creds.json)
```

## Next Steps

- **[Multi-Input Fine-Mapping](multi-input.md)** - Learn how to analyze multiple studies together
- **[Advanced Topics](advanced.md)** - Deep dive into tool-specific parameters and optimization 