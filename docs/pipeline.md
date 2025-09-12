# Complete Pipeline

The `credtools pipeline` command runs the complete fine-mapping analysis workflow in a single command, orchestrating meta-analysis, quality control, and fine-mapping steps automatically. This is the most convenient way to run end-to-end fine-mapping analysis.

## Overview

The pipeline command integrates all major credtools analysis steps into a streamlined workflow. The pipeline process:

- **Automates workflow orchestration** from prepared inputs to final results
- **Handles meta-analysis** with configurable strategies
- **Performs quality control** with optional skip functionality
- **Executes fine-mapping** with full parameter control
- **Manages intermediate files** and directory structure
- **Provides progress tracking** across all analysis steps
- **Ensures consistency** between analysis phases

## When to Use

Use `credtools pipeline` when you have:

- Prepared locus files ready for complete analysis
- Need to run the full meta-analysis → QC → fine-mapping workflow
- Want consistent parameter settings across all steps
- Prefer automated workflow management over manual step execution
- Production analyses requiring reproducible, standardized processing

## Basic Usage

### Standard Multi-Ancestry Pipeline

```bash
credtools pipeline prepared/final_loci_list.txt results/ \--tool multisusie
```

### Single-Ancestry Pipeline

```bash
credtools pipeline prepared/loci_list.txt results/ \
  --meta-method no_meta --tool susie
```

### Skip Quality Control

```bash
credtools pipeline prepared/loci_list.txt results/ \
  --skip-qc --tool finemap
```

## Command Options

### Workflow Control

| Option | Description | Default |
|--------|-------------|---------|
| `--meta-method` / `-m` | Meta-analysis method | meta_all |
| `--skip-qc` / `-q` | Skip quality control step | False |

| `--tool` / `-t` | Fine-mapping tool | susie |

### General Fine-Mapping Options

| Option | Description | Default |
|--------|-------------|---------|
| `--max-causal` / `-c` | Maximum causal variants | 1 |
| `--set-L-by-cojo` / `-sl` | Set max causal by COJO | True |
| `--coverage` / `-cv` | Credible set coverage | 0.95 |
| `--combine-cred` / `-cc` | Credible set combination method | union |
| `--combine-pip` / `-cp` | PIP combination method | max |
| `--jaccard-threshold` / `-j` | Jaccard threshold | 0.1 |

### Tool-Specific Parameters

The pipeline command accepts all tool-specific parameters from the individual fine-mapping tools. Parameters are organized into help panels for each tool:

#### ABF Parameters
- `--var-prior` / `-vp`: Variance prior (default: 0.2)

#### FINEMAP Parameters  
- `--n-iter` / `-ni`: Number of iterations (default: 100000)
- `--n-threads` / `-nt`: Number of threads (default: 1)

#### SuSiE Parameters
- `--max-iter` / `-i`: Maximum iterations (default: 100)
- `--estimate-residual-variance` / `-er`: Estimate residual variance (default: False)
- `--min-abs-corr` / `-mc`: Minimum absolute correlation (default: 0.5)
- `--convergence-tol` / `-ct`: Convergence tolerance (default: 1e-3)

#### RSparsePro Parameters
- `--eps` / `-e`: Convergence criterion (default: 1e-5)
- `--ubound` / `-ub`: Upper bound (default: 100000)
- `--cthres` / `-ct`: Coverage threshold (default: 0.7)
- Various other specialized parameters

#### SuSiEx Parameters
- `--mult-step` / `-ms`: Use multiple steps (default: False)
- `--keep-ambig` / `-ka`: Keep ambiguous SNPs (default: True)
- `--min-purity` / `-mp`: Minimum purity (default: 0.5)
- `--tol` / `-t`: Convergence tolerance (default: 1e-3)

#### multiSuSiE Parameters
- `--rho` / `-r`: Prior correlation between causal variants (default: 0.75)
- `--scaled-prior-variance` / `-spv`: Scaled prior variance (default: 0.2)
- `--standardize` / `-s`: Standardize summary statistics (default: False)
- `--pop-spec-standardization` / `-pss`: Population-specific standardization (default: True)
- `--estimate-prior-variance` / `-epv`: Estimate prior variance (default: True)
- `--estimate-prior-method` / `-epm`: Prior estimation method (default: "early_EM")
- Various other specialized parameters

## Pipeline Workflow

The pipeline executes these steps in sequence:

### 1. Meta-Analysis Phase
- Processes input loci according to specified meta-analysis method
- Creates meta-analyzed summary statistics and LD matrices
- Generates updated loci information files

### 2. Quality Control Phase (Optional)
- Validates meta-analysis outputs
- Checks data consistency and statistical validity
- Filters out loci that fail quality checks
- Creates QC reports and warnings

### 3. Fine-Mapping Phase
- Runs fine-mapping using specified tool and strategy
- Applies all tool-specific parameters
- Generates posterior inclusion probabilities and credible sets
- Creates comprehensive output files

## Expected Output

The pipeline creates a comprehensive directory structure:

```
results/
├── meta/                    # Meta-analysis results
│   ├── meta_all/           # Results by meta method
│   ├── meta_by_population/
│   └── no_meta/
├── qc/                     # Quality control results (if not skipped)
│   ├── summary_report.txt
│   ├── passed_loci_list.txt
│   └── failed_loci.txt
├── finemap/               # Fine-mapping results
│   ├── locus_1/
│   │   ├── pips.txt
│   │   ├── creds.json
│   │   └── tool_output/
│   └── locus_2/
└── pipeline_log.txt       # Complete pipeline log
```

## Examples

### Example 1: Standard Multi-Ancestry Analysis

```bash
# Complete multi-ancestry pipeline with multiSuSiE
credtools pipeline prepared/multi_ancestry_loci.txt results/ \
  --meta-method meta_all \  --tool multisusie \
  --max-causal 3 \
  --coverage 0.95

# Produces complete analysis from meta-analysis through fine-mapping
```

### Example 2: Single-Ancestry High-Performance Analysis

```bash
# Single ancestry with parallel FINEMAP
credtools pipeline prepared/eur_loci.txt results_eur/ \
  --meta-method no_meta \
  --tool finemap \
  --max-causal 5 \
  --n-threads 8 \
  --n-iter 200000

# Optimized for computational performance
```

### Example 3: Conservative Analysis with Strict QC

```bash
# Conservative analysis with thorough quality control
credtools pipeline prepared/loci.txt results_conservative/ \
  --meta-method meta_by_population \
  --tool susie \
  --max-causal 2 \
  --coverage 0.99 \
  --estimate-residual-variance

# More stringent parameters for high-confidence results
```

### Example 4: Fast Screening Analysis

```bash
# Quick analysis skipping QC for initial screening
credtools pipeline prepared/loci.txt results_quick/ \
  --skip-qc \
  --tool abf \
  --max-causal 1 \
  --meta-method meta_all

# Fastest option for preliminary analysis
```

### Example 5: Production Multi-Ancestry Pipeline

```bash
# Full-featured production analysis
credtools pipeline prepared/production_loci.txt results_production/ \
  --meta-method meta_all \  --tool susiex \
  --max-causal 5 \
  --coverage 0.95 \
  --mult-step \
  --min-purity 0.7

# State-of-the-art multi-ancestry fine-mapping
```

## Workflow Comparison

### Pipeline vs Individual Commands

**Use Pipeline When:**
- Running standard analysis workflows
- Want automated step coordination
- Need consistent parameter application
- Prefer single-command execution
- Production or batch processing

**Use Individual Commands When:**
- Need custom intermediate steps
- Want to inspect outputs between phases
- Debugging analysis issues
- Comparing different meta-analysis strategies
- Maximum control over each step

### Example Equivalent Workflows

**Pipeline Command:**
```bash
credtools pipeline loci.txt results/ --tool susie --max-causal 3
```

**Individual Commands:**
```bash
credtools meta loci.txt results/meta/ --meta-method meta_all
credtools qc results/meta/meta_all/loci_list.txt results/qc/
credtools finemap results/qc/passed_loci.txt results/finemap/ \
  --tool susie --max-causal 3
```

## Performance Considerations

### Computational Resources

**Memory Requirements:**
- Depends on largest LD matrix size and number of loci
- Peak usage during fine-mapping phase
- Consider memory-efficient tools for large studies

**CPU Utilization:**
- Meta-analysis: Limited parallelization
- QC: Embarrassingly parallel across loci
- Fine-mapping: Tool-dependent parallelization

**Storage Requirements:**
- Intermediate files can be substantial
- Consider cleanup strategies for large studies
- Plan for both temporary and permanent storage needs

### Optimization Strategies

1. **Choose appropriate tools**: Balance accuracy vs. computational cost
2. **Optimize max causal**: Start conservative, increase if needed
3. **Use tool-specific threading**: Enable when available (FINEMAP, etc.)
4. **Monitor progress**: Check logs for bottlenecks or failures
5. **Plan disk space**: Ensure adequate storage for all outputs

## Troubleshooting

### Common Pipeline Issues

**Pipeline stops at meta-analysis**: Check input file format and loci definitions
**QC failures cause empty fine-mapping**: Use `--skip-qc` or fix data quality issues  
**Fine-mapping convergence failures**: Adjust tool-specific convergence parameters
**Memory errors during pipeline**: Reduce max causal variants or use simpler tools

### Debugging Strategies

1. **Check pipeline logs**: Review complete log file for error messages
2. **Run steps individually**: Isolate which phase is causing issues
3. **Validate inputs**: Ensure input files are properly formatted
4. **Monitor resources**: Check memory and disk usage during execution
5. **Test with subset**: Try pipeline on small subset of loci first

### Recovery from Failures

**Partial pipeline completion**: Resume from last successful step
**Individual locus failures**: Pipeline continues with remaining loci
**Resource exhaustion**: Restart with more conservative parameters
**Data quality issues**: Fix inputs and restart pipeline

## Tips for Success

1. **Start simple**: Begin with default parameters, optimize later
2. **Monitor progress**: Check logs regularly for long-running analyses
3. **Plan resources**: Ensure adequate compute resources before starting
4. **Save intermediate results**: Keep meta-analysis and QC outputs
5. **Document parameters**: Record all parameter choices for reproducibility
6. **Test first**: Run pipeline on small subset before full analysis
7. **Consider alternatives**: Individual commands may be better for complex workflows