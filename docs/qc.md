# Quality Control

The `credtools qc` command performs comprehensive quality control checks on summary statistics and LD matrices to ensure data integrity and suitability for fine-mapping analysis. This step helps identify and resolve potential issues before running computationally intensive fine-mapping.

## Overview

Quality control in credtools validates the consistency and quality of your fine-mapping inputs. The QC process:

- **Validates summary statistics** for missing values, outliers, and format consistency
- **Checks LD matrix integrity** including positive definiteness and reasonable correlation values
- **Verifies variant matching** between summary statistics and LD matrices
- **Identifies allele mismatches** and strand orientation issues
- **Reports data quality metrics** for each locus and ancestry
- **Flags problematic loci** that may cause fine-mapping failures
- **Provides recommendations** for resolving common issues

## When to Use

Use `credtools qc` when you have:

- Prepared fine-mapping inputs that need validation before analysis
- Concerns about data quality or consistency across ancestries
- New datasets or reference panels that haven't been validated
- Previous fine-mapping runs that failed due to data issues
- Need to generate quality reports for publication or sharing

## Basic Usage

### Standard QC Check

```bash
credtools qc prepared/final_loci_list.txt qc_output/
```

### Multi-threaded QC

```bash
credtools qc prepared/final_loci_list.txt qc_output/ --threads 8
```

### QC After Meta-Analysis

```bash
credtools qc meta/meta_all/loci_list.txt qc_output/
```

## Command Options

| Option | Description | Default |
|--------|-------------|---------|
| `--threads` / `-t` | Number of parallel threads | 1 |

## Quality Control Checks

### Summary Statistics Validation

**Format Consistency**
- Required columns present (CHR, BP, EA, NEA, BETA, SE, P)
- Appropriate data types and value ranges
- No missing values in critical columns

**Statistical Validity**
- P-values in valid range (0, 1]
- Standard errors are positive
- Effect sizes within reasonable bounds
- Allele frequencies between 0 and 1 (if present)

**Biological Plausibility**
- Chromosome and position coordinates valid
- Alleles are valid DNA bases
- No duplicate variants within locus

### LD Matrix Validation

**Matrix Properties**
- Square symmetric matrix
- Positive definite (all eigenvalues > 0)
- Diagonal elements equal to 1
- Off-diagonal correlations in [-1, 1]

**Variant Matching**
- All summary statistic variants present in LD matrix
- Consistent variant ordering
- Matching allele orientations

**Data Quality**
- Reasonable correlation structure
- No excessive missing data
- Appropriate matrix conditioning

### Cross-Dataset Consistency

**Multi-Ancestry Checks**
- Consistent variant sets across populations
- Similar allele frequencies where expected
- Reasonable effect size correlations

**Longitudinal Consistency**
- Variant positions match reference genome
- Allele definitions consistent with standards
- No systematic biases across chromosomes

## Expected Output

The QC process creates detailed reports and summaries:

### Quality Control Reports

```
qc_output/
├── summary_report.txt          # Overall QC summary
├── locus_reports/              # Individual locus details
│   ├── locus_1_qc.txt
│   ├── locus_2_qc.txt
│   └── ...
├── failed_loci.txt             # Loci that failed QC
├── warnings.txt                # Non-critical issues found
└── qc_metrics.json             # Machine-readable results
```

### QC Summary Metrics

- **Pass/Fail counts** for each check type
- **Quality scores** for each locus
- **Recommended actions** for failed loci
- **Data completeness** statistics
- **Cross-ancestry comparisons** (if applicable)

## Understanding QC Results

### QC Status Categories

**PASS**: Locus passes all quality checks and is ready for fine-mapping
**WARN**: Minor issues detected but locus may still be usable
**FAIL**: Critical issues that will likely cause fine-mapping failure

### Common Warning Types

- Minor allele frequency differences between populations
- Slightly high correlation values (> 0.99)
- Missing optional columns (EAF, INFO scores)
- Borderline LD matrix conditioning

### Common Failure Types

- Missing required data columns
- Non-positive definite LD matrices
- Severe variant mismatches between datasets
- Invalid statistical values (negative SE, P-values > 1)

## Examples

### Example 1: Basic QC for Single Ancestry

```bash
# Run QC on prepared single-ancestry data
credtools qc prepared/final_loci_list.txt qc_results/

# Review results
cat qc_results/summary_report.txt
cat qc_results/failed_loci.txt
```

### Example 2: Multi-Ancestry QC with Parallel Processing

```bash
# QC multi-ancestry meta-analysis results
credtools qc meta/meta_all/loci_list.txt qc_meta/ --threads 12

# Check for ancestry-specific issues
grep "ancestry" qc_meta/warnings.txt
```

### Example 3: QC After Filtering

```bash
# QC after applying custom filters
credtools qc filtered/loci_list.txt qc_filtered/ --threads 4

# Compare with original QC results
diff qc_original/summary_report.txt qc_filtered/summary_report.txt
```

### Example 4: Iterative QC and Fixing

```bash
# Initial QC
credtools qc prepared/loci_list.txt qc_initial/

# Fix issues based on QC report
# ... manual data cleaning steps ...

# Re-run QC to verify fixes
credtools qc cleaned/loci_list.txt qc_final/
```

## Integration with Workflow

Quality control should be performed after data preparation or meta-analysis:

```bash
# Standard workflow with QC
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/
credtools meta prepared/final_loci_list.txt meta/ --meta-method meta_all
credtools qc meta/meta_all/loci_list.txt qc/ --threads 8
credtools finemap qc/passed_loci_list.txt finemap_results/

# Alternative: QC before meta-analysis
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/
credtools qc prepared/final_loci_list.txt qc_pre_meta/
credtools meta qc_pre_meta/passed_loci_list.txt meta/
credtools finemap meta/meta_all/loci_list.txt finemap_results/
```

## Resolving Common Issues

### Failed LD Matrix Checks

**Non-positive definite matrices**:
```bash
# Usually indicates numerical precision issues
# Try using higher-quality reference panels
# Consider increasing minimum allele frequency filters
```

**Variant mismatches**:
```bash
# Check allele flipping and strand orientation
# Verify reference genome versions match
# Update variant IDs to consistent format
```

### Summary Statistics Issues

**Missing critical columns**:
```bash
# Re-run munging with proper column mapping
# Check original GWAS file format
# Verify required statistics are available
```

**Invalid statistical values**:
```bash
# Check for data corruption during file transfers
# Verify GWAS analysis pipeline outputs
# Look for systematic issues in specific regions
```

### Performance Issues

**Slow QC processing**:
```bash
# Increase thread count up to available cores
# Process subsets of loci in parallel
# Use faster storage for temporary files
```

**Memory consumption**:
```bash
# Reduce thread count for large LD matrices
# Process chromosomes separately
# Consider using high-memory compute nodes
```

## Troubleshooting

### Common Error Messages

**"LD matrix not positive definite"**: The correlation matrix has numerical issues. Check reference panel quality and consider filtering low-frequency variants.

**"Variant mismatch between sumstats and LD"**: Summary statistics and LD matrix contain different variants. Verify data preparation steps.

**"Invalid p-values detected"**: P-values outside valid range [0,1]. Check GWAS analysis pipeline.

**"Memory allocation failed"**: Insufficient RAM for large matrices. Reduce thread count or use smaller chunks.

### Best Practices

1. **Always run QC**: Never skip quality control, especially with new datasets
2. **Review all reports**: Check both summary and detailed locus reports
3. **Fix issues systematically**: Address data quality problems before fine-mapping
4. **Document QC results**: Keep QC reports for reproducibility and troubleshooting
5. **Use appropriate resources**: Allocate sufficient compute resources for large studies

## Tips for Success

1. **Start with QC**: Run quality control early to identify issues before investing in computationally expensive steps
2. **Use parallel processing**: QC can be parallelized effectively across loci
3. **Keep detailed logs**: QC reports are valuable for debugging fine-mapping issues
4. **Iterate as needed**: Re-run QC after fixing issues to ensure problems are resolved
5. **Compare across ancestries**: Use QC to identify systematic differences between populations
6. **Plan for failures**: Some loci may fail QC - have strategies for handling incomplete datasets