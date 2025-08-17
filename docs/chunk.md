# Loci Identification and Chunking

The `credtools chunk` command identifies independent genetic loci from munged summary statistics and splits the data into locus-specific files ready for fine-mapping. This step is crucial for defining the genomic regions that will be analyzed independently.

## Overview

After munging your summary statistics, you need to identify independent genetic signals and create focused datasets for fine-mapping. The chunking process:

- **Identifies independent loci** based on distance and significance thresholds
- **Merges overlapping regions** across different ancestries when appropriate
- **Creates locus-specific files** containing only variants within each region
- **Generates credtools-compatible** loci lists for downstream analysis
- **Handles multi-ancestry data** with consistent loci across populations

## When to Use

Use `credtools chunk` when you have:

- Munged summary statistics ready for fine-mapping
- Genome-wide data that needs to be split into independent regions
- Multi-ancestry studies requiring consistent loci definitions
- Large datasets that benefit from parallel processing of smaller regions

## Basic Usage

### Single Ancestry

```bash
credtools chunk EUR.munged.txt.gz output_dir/
```

### Multiple Ancestries

```bash
credtools chunk "EUR.munged.txt.gz,ASN.munged.txt.gz,AFR.munged.txt.gz" output_dir/
```

### Using File Configuration

Create a JSON file mapping ancestries to munged files:

```json
{
  "EUR": "munged/EUR.munged.txt.gz",
  "ASN": "munged/ASN.munged.txt.gz",
  "AFR": "munged/AFR.munged.txt.gz"
}
```

Then run:

```bash
credtools chunk ancestry_files.json output_dir/
```

## Command Options

| Option | Description | Default |
|--------|-------------|---------|
| `--distance` / `-d` | Distance threshold for independence (bp) | 500000 |
| `--pvalue` / `-p` | P-value threshold for significance | 5e-8 |
| `--merge-overlapping` / `-m` | Merge overlapping loci across ancestries | True |
| `--use-most-sig` / `-u` | Use most significant SNP if no significant SNPs | True |
| `--min-variants` / `-v` | Minimum variants per locus | 10 |
| `--threads` / `-t` | Number of threads | 1 |

## Algorithm Details

### Loci Identification Process

1. **Significance filtering**: Identify SNPs below the p-value threshold
2. **Distance-based clustering**: Group significant SNPs within the distance threshold
3. **Lead SNP selection**: Choose the most significant SNP in each cluster
4. **Region definition**: Create windows around lead SNPs (±distance/2)
5. **Overlap resolution**: Merge overlapping regions across ancestries if requested
6. **Quality filtering**: Remove loci with too few variants

### Multi-Ancestry Handling

When processing multiple ancestries:

- Each ancestry is processed independently first
- Overlapping loci across ancestries can be merged into unified regions
- The most significant lead SNP across all ancestries is selected for merged loci
- Ancestry information is preserved in the output

## Expected Output

The chunking process creates several important files:

### Main Output Files

1. **`identified_loci.txt`** - Summary of all identified loci with coordinates and lead SNPs
2. **`chunks/`** - Directory containing locus-specific summary statistics files
3. **`chunk_info.txt`** - Metadata about all generated chunk files
4. **`loci_list.txt`** - Credtools-compatible loci list for fine-mapping

### Chunk File Structure

Each locus generates ancestry-specific files:
```
chunks/
├── EUR.chr1_12345_67890.sumstats.gz
├── ASN.chr1_12345_67890.sumstats.gz
├── AFR.chr1_12345_67890.sumstats.gz
└── ...
```

### Loci List Format

The `loci_list.txt` file contains:

| Column | Description |
|--------|-------------|
| locus_id | Unique identifier (chr_start_end) |
| chr | Chromosome number |
| start | Start position (bp) |
| end | End position (bp) |
| popu | Population/ancestry code |
| cohort | Cohort identifier |
| sample_size | Sample size (placeholder) |
| prefix | File prefix for credtools |

## Examples

### Example 1: Conservative Loci Definition

```bash
# Use stricter thresholds for fewer, more significant loci
credtools chunk munged_files.json output/ \
  --distance 1000000 \
  --pvalue 1e-8 \
  --min-variants 50
```

### Example 2: Liberal Loci Definition

```bash
# Use relaxed thresholds for more comprehensive coverage
credtools chunk munged_files.json output/ \
  --distance 250000 \
  --pvalue 1e-5 \
  --min-variants 5
```

### Example 3: Population-Specific Analysis

```bash
# Don't merge overlapping loci across ancestries
credtools chunk ancestry_files.json output/ \
  --merge-overlapping false
```

### Example 4: Parallel Processing

```bash
# Use multiple threads for faster processing
credtools chunk large_dataset.json output/ \
  --threads 8
```

## Parameter Guidelines

### Distance Threshold (`--distance`)

- **500kb (default)**: Balanced approach, suitable for most analyses
- **250kb**: More granular loci, better for dense association regions
- **1Mb**: Conservative approach, reduces computational burden
- **Consider LD patterns**: Longer distance in populations with extended LD

### P-value Threshold (`--pvalue`)

- **5e-8 (default)**: Genome-wide significance threshold
- **1e-5**: Suggestive significance, more inclusive
- **1e-8**: Very stringent, fewer but highly significant loci
- **Population-specific**: Consider ancestry-specific significance levels

### Minimum Variants (`--min-variants`)

- **10 (default)**: Ensures reasonable LD computation
- **5**: More inclusive, useful for sparse regions
- **20-50**: Conservative, better for high-quality fine-mapping

## Integration with Workflow

Chunked data feeds directly into the preparation step:

```bash
# 1. Munge summary statistics  
credtools munge ancestry_files.json munged/

# 2. Identify loci and chunk data
credtools chunk munged/ chunked/

# 3. Prepare LD matrices (uses chunk_info.txt)
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/

# 4. Run fine-mapping (uses final_loci_list.txt)
credtools finemap prepared/final_loci_list.txt results/
```

## Quality Control

### Reviewing Loci

After chunking, examine the results:

```bash
# Check number of loci identified
wc -l chunked/identified_loci.txt

# Review loci summary
head -20 chunked/identified_loci.txt

# Check chunk file counts
ls chunked/chunks/ | wc -l
```

### Ancestry Coverage

Verify consistent coverage across ancestries:

```bash
# Count chunks per ancestry
cut -f2 chunked/chunk_info.txt | sort | uniq -c
```

## Troubleshooting

### Common Issues

**No loci identified**: Lower the p-value threshold or check that munged files contain significant associations.

**Too many small loci**: Increase the distance threshold or minimum variants requirement.

**Memory issues with large datasets**: Use parallel processing with `--threads` and ensure sufficient RAM.

**Inconsistent loci across ancestries**: Check that munged files use consistent chromosome and position formats.

### Performance Optimization

**Large datasets**: Use more threads (`--threads 4-8`) for faster processing.

**Memory constraints**: Process ancestries separately by providing single files rather than multiple files.

**Storage space**: Consider the trade-off between number of loci and storage requirements for chunk files.

## Advanced Usage

### Custom Significance Levels

For population-specific analysis, you might use different p-value thresholds:

```bash
# European ancestry (well-powered)
credtools chunk EUR.munged.txt.gz eur_chunks/ --pvalue 5e-8

# Smaller ancestries (more lenient)  
credtools chunk ASN.munged.txt.gz asn_chunks/ --pvalue 1e-6
```

### Region-Specific Analysis

To focus on specific genomic regions:

```bash
# Pre-filter summary statistics to chromosome 1 
zcat EUR.munged.txt.gz | awk '$1==1' | gzip > EUR_chr1.munged.txt.gz
credtools chunk EUR_chr1.munged.txt.gz chr1_chunks/
```

## Tips for Success

1. **Start with defaults**: Use default parameters initially, then optimize based on results
2. **Consider your goals**: More loci = more comprehensive but computationally intensive
3. **Check overlap**: Review identified loci for known associations in your trait
4. **Balance precision vs coverage**: Stricter thresholds give higher confidence but may miss signals
5. **Document parameters**: Keep track of thresholds used for reproducibility
6. **Validate with literature**: Compare identified loci with known associations for your trait