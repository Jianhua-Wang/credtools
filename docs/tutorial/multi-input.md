# Multi-Input Analysis

Multi-input fine-mapping analyzes multiple cohorts and ancestries simultaneously to leverage shared genetic architecture while accounting for population differences. This is where CREDTOOLS truly shines, offering sophisticated approaches to multi-ancestry genetic analysis.

## When to Use Multi-Input Strategy

!!! tip "Multi-Input Use Cases"
    
    **Multiple ancestries available**
    : Leverage power across populations while modeling LD differences
    
    **Shared causal architecture**
    : When you expect similar causal variants across populations
    
    **Increased statistical power**
    : Combine sample sizes for improved fine-mapping resolution
    
    **Cross-population validation**
    : Identify variants with consistent effects across ancestries

## Multi-Input Workflow Components

### 1. Meta-Analysis Strategies

The first step is deciding how to combine your studies:

#### Cross-Ancestry Meta-Analysis (`meta_all`)

Combines all studies regardless of ancestry:

```bash
credtools meta input_loci.txt meta_output/ \
    --meta-method meta_all \
    --threads 4
```

**When to use:**
- Strong prior belief in shared causal variants
- Large effect sizes relative to population differences
- Increased power is the primary goal

**What it does:**
- Performs inverse-variance weighted meta-analysis of summary statistics
- Sample-size weighted averaging of LD matrices
- Creates single meta-analyzed dataset per locus

#### Within-Ancestry Meta-Analysis (`meta_by_population`)

Combines studies within each ancestry separately:

```bash
credtools meta input_loci.txt meta_output/ \
    --meta-method meta_by_population \
    --threads 4
```

**When to use:**
- Population-specific effect sizes expected
- Want to preserve ancestry-specific signals
- Balanced approach between power and specificity

**What it does:**
- Meta-analyzes EUR studies together, AFR studies together, etc.
- Maintains population-specific LD structure
- Creates separate datasets for each ancestry

#### No Meta-Analysis (`no_meta`)

Keeps all studies separate:

```bash
credtools meta input_loci.txt meta_output/ \
    --meta-method no_meta \
    --threads 4
```

**When to use:**
- Highly heterogeneous effect sizes
- Study-specific technical factors
- Maximum preservation of individual study characteristics

### 2. Quality Control Across Studies

Multi-study QC provides insights impossible with single studies:

```bash
credtools qc input_loci.txt qc_output/ --threads 4
```

#### Cross-Study QC Metrics

**Cochran's Q Test**
: Tests for heterogeneity in effect sizes across studies
: High Q values suggest population or study-specific effects

**SNP Missingness Patterns**
: Identifies variants missing in specific populations
: Helps understand coverage differences across ancestries

**LD Structure Comparison**
: Compares LD patterns across populations
: Identifies regions with dramatically different LD structure

**Cross-Ancestry MAF Correlations**
: Compares allele frequencies across populations
: Detects potential population stratification or technical issues

!!! info "QC Output Files"
    ```
    qc_output/
    ├── cochran_q.txt           # Effect size heterogeneity
    ├── snp_missingness.txt     # Coverage patterns
    ├── ld_4th_moment.txt       # LD structure comparison
    ├── s_estimate.txt          # Consistency parameters
    └── maf_comparison.txt      # Frequency comparisons
    ```

### 3. Multi-Input Fine-Mapping Approaches

CREDTOOLS automatically selects the appropriate approach for multi-input fine-mapping:

#### A. Multi-Input Tools

Use tools specifically designed for multi-population analysis:

```bash
credtools finemap input_loci.txt output/ \
    --tool multisusie \
    --max-causal 5
```

**Supported Tools:**

- **MultiSuSiE**: Multi-population extension of SuSiE
- **SuSiEx**: Cross-ancestry fine-mapping tool

#### B. Post-hoc Combination

Run single-input tools on each study, then intelligently combine results:

```bash
credtools finemap input_loci.txt output/ \
    --tool susie \
    --combine-cred union \
    --combine-pip max \
    --jaccard-threshold 0.1
```

**Combination Methods:**

!!! note "Credible Set Combination (`--combine-cred`)"
    
    **`union`** (default)
    : Union of all credible sets from individual studies
    : Most inclusive, captures all potential causal variants
    
    **`intersection`**
    : Only variants present in credible sets from all studies
    : Most conservative, highest confidence variants only
    
    **`cluster`**
    : Groups overlapping credible sets using Jaccard similarity
    : Balanced approach, creates meta-credible sets

!!! note "PIP Combination (`--combine-pip`)"
    
    **`max`** (default)
    : Maximum PIP across all studies for each variant
    : Emphasizes strongest signals
    
    **`mean`**
    : Average PIP across studies
    : Balanced view of evidence
    
    **`meta`**
    : Meta-analysis formula: 1 - ∏(1 - PIP_i)
    : Accounts for independence assumption

## Detailed Tool Descriptions

### MultiSuSiE

**Best for**: Multi-population analysis with shared and population-specific effects

```bash
credtools finemap input_loci.txt output/ \
    --tool multisusie \
    --max-causal 10 \
    --rho 0.75 \
    --scaled-prior-variance 0.2 \
    --pop-spec-standardization \
    --estimate-prior-variance
```

**Key Parameters:**

- `--rho`: Prior correlation between causal variants across populations (0.75)
- `--scaled-prior-variance`: Prior effect size variance (0.2)  
- `--pop-spec-standardization`: Use population-specific standardization
- `--estimate-prior-variance`: Estimate rather than fix prior variance
- `--pop-spec-effect-priors`: Population-specific effect size priors

!!! tip "MultiSuSiE Guidance"
    - Higher `--rho` assumes more sharing across populations
    - Use `--pop-spec-standardization` when sample sizes vary greatly
    - `--estimate-prior-variance` is usually recommended

### SuSiEx  

**Best for**: Cross-ancestry fine-mapping with explicit modeling of population differences

```bash
credtools finemap input_loci.txt output/ \
    --tool susiex \
    --max-causal 5 \
    --mult-step \
    --keep-ambig \
    --min-purity 0.5 \
    --tol 1e-3
```

**Key Parameters:**

- `--mult-step`: Use multiple refinement steps  
- `--keep-ambig`: Keep ambiguous SNPs in analysis
- `--min-purity`: Minimum purity for credible sets (0.5)
- `--tol`: Convergence tolerance (1e-3)

## Comprehensive Multi-Input Example

Here's a complete workflow for multi-ancestry fine-mapping:

### Step 1: Input Preparation

```bash title="multi_ancestry_loci.txt"
chr	start	end	popu	sample_size	cohort	prefix	locus_id
8	41242482	42492482	AFR	89499	MVP	data/AFR.MVP.chr8_41242482_42492482	chr8_41242482_42492482
8	41242482	42492482	EUR	337465	MVP	data/EUR.MVP.chr8_41242482_42492482	chr8_41242482_42492482
8	41242482	42492482	EUR	442817	UKBB	data/EUR.UKBB.chr8_41242482_42492482	chr8_41242482_42492482
8	41242482	42492482	SAS	8253	UKBB	data/SAS.UKBB.chr8_41242482_42492482	chr8_41242482_42492482
```

### Step 2: Quality Control

```bash
# Run comprehensive QC
credtools qc multi_ancestry_loci.txt qc_results/ --threads 4

# Review heterogeneity
head qc_results/cochran_q.txt
```

### Step 3A: Cross-Ancestry Analysis with MultiSuSiE

```bash
# Meta-analyze across all populations
credtools pipeline multi_ancestry_loci.txt results_cross_ancestry/ \
    --meta-method meta_all \
    --tool multisusie \
    --max-causal 10 \
    --rho 0.8 \
    --pop-spec-standardization \
    --estimate-prior-variance
```

### Step 3B: Population-Specific Analysis

```bash
# Meta-analyze within populations, then combine
credtools pipeline multi_ancestry_loci.txt results_pop_specific/ \
    --meta-method meta_by_population \
    --tool susie \
    --combine-cred cluster \
    --combine-pip meta \
    --jaccard-threshold 0.2
```

### Step 4: Compare Approaches

```python
# Compare results from different strategies
import json
import pandas as pd

# Load results
with open('results_cross_ancestry/creds.json') as f:
    cross_ancestry = json.load(f)
    
with open('results_pop_specific/creds.json') as f:
    pop_specific = json.load(f)

# Compare credible sets and PIPs
cross_pips = pd.read_csv('results_cross_ancestry/pips.txt', 
                        sep='\t', header=None, names=['SNP', 'PIP'])
pop_pips = pd.read_csv('results_pop_specific/pips.txt', 
                      sep='\t', header=None, names=['SNP', 'PIP'])
```

## Advanced Multi-Input Considerations

### Handling Missing Data

When variants are missing in some populations:

!!! warning "Missing Data Strategies"
    
    **Complete case analysis**
    : Only use variants present in all studies
    : Reduces power but ensures consistency
    
    **Imputation**
    : Impute missing variants using population-specific references
    : Requires careful validation
    
    **Weighted analysis**
    : Weight contributions by data availability
    : Built into CREDTOOLS meta-analysis

### Population Stratification

When populations have internal structure:

```bash
# Use more conservative thresholds
credtools finemap input.txt output/ \
    --tool multisusie \
    --rho 0.5 \
    --min-abs-corr 0.8
```

### Computational Considerations

Large multi-population analyses can be computationally intensive:

!!! tip "Performance Optimization"
    
    - Use `--threads` for parallelization
    - Consider analyzing loci separately for very large studies
    - Monitor memory usage with many populations
    - Use post-hoc combination for maximum flexibility

## Interpreting Multi-Input Results

### Population-Specific Effects

Look for variants with:
- High PIPs in some populations but not others
- Different credible sets across populations
- High Cochran's Q values

### Shared Causal Variants

Evidence includes:
- Consistent PIPs across populations
- Overlapping credible sets
- Low heterogeneity in meta-analysis

### Meta-Analysis Benefits

Compare single-population vs. meta-analysis:
- Increased resolution (smaller credible sets)
- Higher PIPs for true causal variants
- Discovery of additional signals

## Common Multi-Input Issues

!!! warning "Troubleshooting Multi-Input Analysis"
    
    **Excessive heterogeneity**
    : Use `meta_by_population` instead of `meta_all`
    : Check for technical differences between studies
    : Consider population-specific analyses
    
    **No shared signals**
    : May indicate true population-specific effects
    : Check LD structure differences
    : Verify consistent variant coding
    
    **Computational issues**
    : Reduce number of populations analyzed together
    : Use post-hoc combination strategy
    : Increase available memory/compute resources

## Next Steps

- **[Advanced Topics](advanced.md)** - Deep dive into parameter optimization and custom workflows
- Review specific tool documentation for detailed parameter guidance
- Consider population genetics factors in result interpretation 