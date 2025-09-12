# Advanced Topics

This section covers advanced CREDTOOLS usage including detailed parameter optimization, custom workflows, and troubleshooting complex scenarios.

## Tool-Specific Parameter Optimization

### SuSiE Advanced Parameters

SuSiE is the most commonly used tool in CREDTOOLS. Here's how to optimize its performance:

```bash
credtools finemap input.txt output/ \    --tool susie \
    --max-causal 10 \
    --max-iter 100 \
    --estimate-residual-variance \
    --min-abs-corr 0.5 \
    --convergence-tol 1e-3 \
    --coverage 0.95
```

!!! tip "SuSiE Parameter Tuning"
    
    **`--max-causal`** (L parameter)
    : Start with COJO-estimated value (`--set-L-by-cojo`)
    : Increase if credible sets seem too restrictive
    : Rule of thumb: L ≈ number of genome-wide significant hits in region
    
    **`--max-iter`**
    : Default 100 is usually sufficient
    : Increase to 500+ for complex regions with many causal variants
    : Monitor convergence warnings
    
    **`--estimate-residual-variance`**
    : Use `True` when phenotype variance is unknown
    : Use `False` for standardized effect sizes
    : Can improve convergence in some cases
    
    **`--min-abs-corr`**
    : Minimum correlation threshold for credible sets
    : Lower values (0.1-0.3) for diverse LD patterns
    : Higher values (0.5-0.8) for strong LD regions

### FINEMAP Advanced Configuration

FINEMAP offers extensive Bayesian model configuration:

```bash
credtools finemap input.txt output/ \    --tool finemap \
    --max-causal 5 \
    --n-iter 1000000 \
    --n-threads 8 \
    --coverage 0.95
```

!!! note "FINEMAP Considerations"
    
    **Computational Requirements**
    : Memory usage scales with region size² 
    : Use more iterations (1M+) for stable results
    : Parallel threads improve speed significantly
    
    **Model Space Exploration**
    : FINEMAP explores all possible causal combinations
    : Exponential complexity limits max-causal to ~5-8
    : Consider region subdivision for larger L

### MultiSuSiE Population Parameters

Fine-tune multi-population analysis:

```bash
credtools finemap input.txt output/ \    --tool multisusie \
    --max-causal 10 \
    --rho 0.75 \
    --scaled-prior-variance 0.2 \
    --pop-spec-standardization \
    --estimate-prior-variance \
    --pop-spec-effect-priors \
    --iter-before-zeroing-effects 5 \
    --prior-tol 1e-9
```

!!! tip "MultiSuSiE Optimization"
    
    **`--rho` (correlation parameter)**
    : 0.9-0.95: Strong sharing across populations (most traits)
    : 0.7-0.8: Moderate sharing with some population specificity
    : 0.5-0.6: Weak sharing, mostly population-specific effects
    
    **Population-specific options**
    : Use `--pop-spec-standardization` when sample sizes vary >10x
    : Use `--pop-spec-effect-priors` for very different populations
    : Monitor convergence with complex population structure

### CARMA Model Selection

CARMA offers sophisticated model uncertainty quantification:

```bash
credtools finemap input.txt output/ \    --tool carma \
    --max-causal 10 \
    --effect-size-prior "Spike-slab" \
    --y-var 1.0 \
    --bf-threshold 10.0 \
    --outlier-bf-threshold 0.31 \
    --max-model-dim 200000 \
    --tau 0.04 \
    --em-dist "logistic"
```

!!! note "CARMA Model Configuration"
    
    **Prior Selection**
    : "Spike-slab" for sparse genetic architecture
    : "Cauchy" for more diffuse effects
    
    **Bayes Factor Thresholds**
    : Lower `--bf-threshold` includes more models
    : Higher values focus on strongest evidence
    : `--outlier-bf-threshold` controls outlier detection

## Custom Workflow Development

### Manual Pipeline Execution

For maximum control, run pipeline steps separately:

```bash
# Step 1: Meta-analysis
credtools meta input_loci.txt meta_output/ \
    --meta-method meta_by_population \
    --threads 4

# Step 2: Quality control  
credtools qc meta_output/updated_loci.txt qc_output/ \
    --threads 4

# Step 3: Review QC and filter problematic studies
# (Manual inspection of QC outputs)

# Step 4: Fine-mapping with optimized parameters
credtools finemap filtered_loci.txt final_output/ \    --tool susie \
    --max-causal 15 \
    --combine-cred cluster \
    --jaccard-threshold 0.2
```

### Parallel Processing for Multiple Loci

Process many loci efficiently:

```bash
# Create per-locus input files
split -l 2 all_loci.txt locus_

# Process in parallel
for locus_file in locus_*; do
    locus_id=$(tail -1 $locus_file | cut -f8)
    credtools pipeline $locus_file results/$locus_id/ \
        --tool susie --max-causal 5 &
done
wait

# Combine results
find results/ -name "pips.txt" -exec cat {} \; > combined_pips.txt
```

### Custom Combination Strategies

Implement custom result combination:

```python
import pandas as pd
import numpy as np
from pathlib import Path

def custom_pip_combination(pip_files, method="harmonic_mean"):
    """Custom PIP combination across studies."""
    all_pips = []
    
    for file in pip_files:
        pips = pd.read_csv(file, sep='\t', header=None, 
                          names=['SNP', 'PIP'], index_col=0)
        all_pips.append(pips)
    
    # Align all studies to same SNP set
    common_snps = set.intersection(*[set(p.index) for p in all_pips])
    aligned_pips = [p.loc[common_snps] for p in all_pips]
    
    if method == "harmonic_mean":
        # Harmonic mean of PIPs
        pip_matrix = pd.concat(aligned_pips, axis=1)
        combined = len(pip_matrix.columns) / (1/pip_matrix).sum(axis=1)
    elif method == "geometric_mean":
        # Geometric mean of PIPs  
        pip_matrix = pd.concat(aligned_pips, axis=1)
        combined = pip_matrix.prod(axis=1) ** (1/pip_matrix.shape[1])
    
    return combined.sort_values(ascending=False)
```

## Quality Control Deep Dive

### Interpreting QC Metrics

#### S Parameter Interpretation

The s parameter measures data quality:

```python
import pandas as pd

s_estimates = pd.read_csv('qc_output/s_estimate.txt', sep='\t')
print(f"Mean s: {s_estimates['s'].mean():.3f}")
print(f"Max s: {s_estimates['s'].max():.3f}")

# Flag problematic studies
problematic = s_estimates[s_estimates['s'] > 0.2]
print(f"Problematic studies: {len(problematic)}")
```

!!! warning "S Parameter Thresholds"
    
    - **s < 0.1**: Excellent data quality
    - **0.1 < s < 0.2**: Acceptable quality  
    - **s > 0.2**: Potential issues (investigate further)
    - **s > 0.5**: Serious problems (consider excluding)

#### Cochran's Q Analysis

Assess effect size heterogeneity:

```python
cochran_q = pd.read_csv('qc_output/cochran_q.txt', sep='\t')

# High Q indicates heterogeneity
high_het = cochran_q[cochran_q['Q_pval'] < 0.05]
print(f"SNPs with significant heterogeneity: {len(high_het)}")

# Examine I² statistic
mean_i2 = cochran_q['I2'].mean()
print(f"Mean I²: {mean_i2:.1f}%")
```

### Advanced QC Filtering

```bash
# Filter studies based on QC metrics
python filter_studies.py \
    --input original_loci.txt \
    --s-threshold 0.2 \
    --het-threshold 0.05 \
    --output filtered_loci.txt
```

```python title="filter_studies.py"
#!/usr/bin/env python3
import pandas as pd
import argparse

def filter_studies(input_file, s_threshold, het_threshold, output_file):
    # Load original loci
    loci = pd.read_csv(input_file, sep='\t')
    
    # Load QC metrics
    s_est = pd.read_csv('qc_output/s_estimate.txt', sep='\t')
    
    # Filter based on s parameter
    good_studies = s_est[s_est['s'] <= s_threshold]['study_id']
    filtered_loci = loci[loci['prefix'].isin(good_studies)]
    
    # Additional filtering logic here...
    
    filtered_loci.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered from {len(loci)} to {len(filtered_loci)} studies")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--s-threshold', type=float, default=0.2)
    parser.add_argument('--het-threshold', type=float, default=0.05)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    filter_studies(args.input, args.s_threshold, args.het_threshold, args.output)
```

## Performance Optimization

### Memory Management

For large regions or many populations:

```bash
# Monitor memory usage
credtools finemap input.txt output/ \
    --tool susie \
    --max-causal 5 \
    --verbose 2>&1 | grep -i memory

# Reduce memory footprint
export OMP_NUM_THREADS=1  # Limit thread memory
ulimit -v 32000000       # Set memory limit (32GB)
```

### Computational Scaling

```bash
# Scale across multiple nodes
for chr in {1..22}; do
    sbatch --job-name=chr${chr} run_chr.sh $chr
done
```

```bash title="run_chr.sh"
#!/bin/bash
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

chr=$1
credtools pipeline chr${chr}_loci.txt results/chr${chr}/ \
    --threads 8 \
    --tool susie \
    --max-causal 10
```

## Troubleshooting Complex Issues

### Convergence Problems

When tools fail to converge:

!!! warning "Convergence Issues"
    
    **SuSiE not converging**
    : Increase `--max-iter` to 500+
    : Try `--estimate-residual-variance`
    : Reduce `--max-causal` if very high
    : Check for numerical instability in LD matrix
    
    **FINEMAP stuck**
    : Increase `--n-iter` substantially
    : Reduce `--max-causal` 
    : Check for perfect correlation in LD matrix
    
    **MultiSuSiE unstable**
    : Reduce `--rho` if populations very different
    : Use `--pop-spec-effect-priors`
    : Check population-specific sample sizes

### LD Matrix Issues

Common LD matrix problems:

```python
import numpy as np
from scipy.linalg import LinAlgError

def diagnose_ld_matrix(ld_file):
    """Diagnose LD matrix problems."""
    ld = np.load(ld_file)['ld']
    
    # Check for numerical issues
    print(f"Matrix shape: {ld.shape}")
    print(f"Diagonal range: {np.diag(ld).min():.3f} - {np.diag(ld).max():.3f}")
    print(f"Off-diagonal range: {ld[~np.eye(ld.shape[0], dtype=bool)].min():.3f} - {ld[~np.eye(ld.shape[0], dtype=bool)].max():.3f}")
    
    # Check positive definiteness
    try:
        eigvals = np.linalg.eigvals(ld)
        min_eigval = eigvals.min()
        print(f"Minimum eigenvalue: {min_eigval:.6f}")
        
        if min_eigval < -1e-8:
            print("WARNING: Matrix not positive semidefinite")
        
        # Check condition number
        cond_num = np.linalg.cond(ld)
        print(f"Condition number: {cond_num:.2e}")
        
        if cond_num > 1e12:
            print("WARNING: Matrix is near-singular")
            
    except LinAlgError:
        print("ERROR: Cannot compute eigenvalues")

# Usage
diagnose_ld_matrix('data/EUR.UKBB.chr8_41242482_42492482.ld.npz')
```

### Missing Data Patterns

Handle systematic missingness:

```python
def analyze_missingness(loci_file):
    """Analyze variant missingness patterns across studies."""
    import pandas as pd
    from pathlib import Path
    
    loci = pd.read_csv(loci_file, sep='\t')
    
    all_variants = set()
    study_variants = {}
    
    for _, row in loci.iterrows():
        # Load variant list for each study
        ldmap_file = f"{row['prefix']}.ldmap"
        if Path(ldmap_file).exists():
            variants = pd.read_csv(ldmap_file, sep='\t')['SNPID'].tolist()
            study_variants[row['prefix']] = set(variants)
            all_variants.update(variants)
    
    # Create missingness matrix
    missingness = pd.DataFrame(index=sorted(all_variants), 
                             columns=study_variants.keys())
    
    for study, variants in study_variants.items():
        missingness[study] = missingness.index.isin(variants)
    
    # Summary statistics
    variant_coverage = missingness.sum(axis=1)
    study_coverage = missingness.sum(axis=0)
    
    print(f"Total variants: {len(all_variants)}")
    print(f"Variants in all studies: {sum(variant_coverage == len(study_variants))}")
    print(f"Study coverage range: {study_coverage.min()} - {study_coverage.max()}")
    
    return missingness

# Usage
missingness = analyze_missingness('input_loci.txt')
```

## Best Practices Summary

### Analysis Strategy Selection

!!! tip "Strategy Decision Tree"
    
    **Single well-powered study** →  `susie`
    
    **Multiple studies, similar populations** → `meta_all` +  `multisusie`
    
    **Multiple studies, different populations** → `meta_by_population` +  `susie`
    
    **Exploratory analysis** → `no_meta` +  `susie`
    
    **Research/publication** → Multiple strategies + comparison

### Parameter Selection Guidelines

1. **Start conservative**: Use default parameters initially
2. **Validate with simulations**: Test on simulated data when possible
3. **Compare strategies**: Run multiple approaches and compare
4. **Document choices**: Keep detailed records of parameter decisions
5. **Iterate based on results**: Adjust parameters based on initial findings

### Quality Control Workflow

1. **Always run QC first**: Never skip quality control
2. **Manual review**: Don't rely solely on automated flags
3. **Population-specific checks**: Different populations may have different issues
4. **Iterative filtering**: Remove problematic studies and re-run
5. **Document exclusions**: Keep records of why studies were excluded

This completes the comprehensive CREDTOOLS tutorial series. Each section builds upon the previous ones to provide a complete guide to multi-ancestry fine-mapping with CREDTOOLS. 