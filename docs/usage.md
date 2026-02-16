# Usage Overview

credtools provides a comprehensive suite of commands for multi-ancestry fine-mapping analysis. This page provides an overview of all available commands and their typical usage patterns.

## Command Summary

credtools includes 6 main subcommands that can be used individually or as part of integrated workflows:

### Data Preparation Commands

#### [munge](munge.md) - Summary Statistics Munging
Standardizes GWAS summary statistics from various formats into credtools-compatible format.

```bash
credtools munge input_files.json output_dir/
```

**Use for:** Converting raw GWAS data, handling multi-ancestry studies, standardizing column formats.

#### [chunk](chunk.md) - Loci Identification and Chunking  
Identifies independent genetic loci and splits data into locus-specific files.

```bash
credtools chunk munged_files.json chunked_output/
```

**Use for:** Defining independent loci, creating analysis-ready chunks, handling genome-wide data.

### Analysis Commands

#### [meta](meta.md) - Meta-Analysis
Combines summary statistics and LD matrices across ancestries or studies.

```bash
credtools meta prepared_loci.txt meta_output/ --meta-method meta_all
```

**Use for:** Multi-ancestry meta-analysis, combining studies, increasing statistical power.

#### [qc](qc.md) - Quality Control
Performs comprehensive quality checks on fine-mapping inputs.

```bash
credtools qc loci_list.txt qc_output/
```

**Use for:** Validating data quality, identifying problematic loci, ensuring analysis reliability.

#### [finemap](finemap.md) - Fine-Mapping Analysis
Runs statistical fine-mapping to identify causal variants.

```bash
credtools finemap loci_list.txt finemap_output/ --tool susie --max-causal 3
```

**Use for:** Identifying causal variants, calculating posterior probabilities, generating credible sets.

### Workflow Commands

#### [pipeline](pipeline.md) - Complete Pipeline
Runs the full analysis workflow in a single command.

```bash
credtools pipeline prepared_loci.txt results/ --tool multisusie```

**Use for:** Automated end-to-end analysis, production workflows, consistent parameter application.


## Typical Workflows

### Single-Ancestry Analysis

```bash
# 1. Standardize summary statistics
credtools munge gwas_eur.txt munged/

# 2. Identify loci, chunk data, and extract LD matrices
#    (include ld_ref column in population config for automatic LD extraction)
credtools chunk population_config.txt chunked/

# 3. Run fine-mapping
credtools finemap chunked/loci_list.txt results/ --tool susie
```

### Multi-Ancestry Analysis

```bash
# 1. Munge all ancestries
credtools munge population_config.txt munged/

# 2. Identify shared loci, chunk data, and extract LD matrices
#    (population config with ld_ref column handles everything)
credtools chunk munged/sumstat_info_updated.txt chunked/ --merge-overlapping

# 3. Run complete pipeline with meta-analysis
credtools pipeline chunked/loci_list.txt results/ \
  --meta-method meta_all --tool multisusie
```

### Quality-Focused Workflow

```bash
# Standard preparation steps
credtools munge population_config.txt munged/
credtools chunk munged/sumstat_info_updated.txt chunked/

# Meta-analysis with quality control
credtools meta chunked/loci_list.txt meta/
credtools qc meta/meta_all/loci_list.txt qc/

# Fine-mapping only on QC-passed loci
credtools finemap qc/passed_loci_list.txt finemap/
```

### Comparative Analysis Workflow

```bash
# Prepare data once
credtools munge population_config.txt munged/
credtools chunk munged/sumstat_info_updated.txt chunked/

# Compare different meta-analysis strategies
for method in meta_all meta_by_population no_meta; do
  credtools pipeline chunked/loci_list.txt results_${method}/ \
    --meta-method $method --tool susie
done

# Compare different fine-mapping tools
for tool in susie abf finemap multisusie; do
  credtools finemap meta/meta_all/loci_list.txt results_${tool}/ --tool $tool
done

# Compare results in results_comparison/ directory
```

## Command Selection Guide

### Choose Commands Based on Your Needs

**Starting with raw GWAS data?**
→ Begin with `munge` to standardize formats

**Have genome-wide association data?**
→ Use `chunk` to identify independent loci

**Need LD information?**
→ Include `ld_ref` in your population config and `chunk` will extract LD matrices automatically

**Multiple ancestries or studies?**
→ Use `meta` to combine evidence appropriately

**Concerned about data quality?**
→ Use `qc` to validate inputs before analysis

**Ready for fine-mapping?**
→ Use `finemap` for detailed variant-level analysis

**Want automated workflow?**
→ Use `pipeline` for streamlined end-to-end analysis

**Need to explore results?**
→ Examine JSON and CSV output files

### Individual Commands vs Pipeline

**Use individual commands when:**
- Learning the credtools workflow
- Need custom intermediate processing
- Debugging analysis issues
- Comparing different strategies
- Maximum control over each step

**Use pipeline when:**
- Running standard analysis workflows
- Production or batch processing
- Want consistent parameter application
- Automated processing pipelines
- Time-sensitive analyses

## Common Parameter Patterns

### Computational Resources

```bash
# High-performance multi-threading
--threads 16

# Memory-efficient processing  
--threads 4

# Single-threaded for debugging
--threads 1
```

### Multi-Ancestry Strategies

```bash
# Maximum power through combination
--meta-method meta_all# Population-specific analysis
--meta-method meta_by_population# Individual ancestry analysis
--meta-method no_meta```

### Fine-Mapping Configuration

```bash
# Conservative single-signal analysis
--tool abf --max-causal 1

# Standard multi-signal analysis
--tool susie --max-causal 3

# Complex multi-signal analysis
--tool susie --max-causal 10 --estimate-residual-variance

# State-of-the-art multi-ancestry
--tool multisusie--max-causal 5
```

## Getting Help

### Command-Specific Help

```bash
# Get help for any command
credtools munge --help
credtools finemap --help
credtools pipeline --help

# See all available commands
credtools --help
```

### Detailed Documentation

Each command has comprehensive documentation:
- [munge](munge.md) - Summary statistics munging
- [chunk](chunk.md) - Loci identification and chunking  
- [meta](meta.md) - Meta-analysis
- [qc](qc.md) - Quality control
- [finemap](finemap.md) - Fine-mapping analysis
- [pipeline](pipeline.md) - Complete pipeline

### Additional Resources

- [Tutorial](tutorial.md) - Step-by-step guides
- [API Documentation](API/) - Programmatic usage
- [Installation](installation.md) - Setup instructions
- [Contributing](contributing.md) - Development guidelines

## Tips for Success

1. **Start simple**: Begin with default parameters and single ancestry
2. **Plan your workflow**: Understand which commands you'll need before starting
3. **Use appropriate resources**: Match computational requirements to available hardware
4. **Validate inputs**: Use QC command to catch issues early
5. **Save intermediate results**: Keep outputs from each step for troubleshooting
6. **Document parameters**: Record command-line options for reproducibility
7. **Explore results thoroughly**: Examine output files and summary statistics