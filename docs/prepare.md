# LD Matrix Preparation

The `credtools prepare` command extracts linkage disequilibrium (LD) matrices from genotype data and creates the final input files needed for fine-mapping analysis. This step bridges chunked summary statistics with reference genotype data.

## Overview

Fine-mapping requires accurate LD information to model the correlation structure between variants. The preparation process:

- **Extracts LD matrices** from reference genotype panels for each locus
- **Matches variants** between summary statistics and genotype data
- **Handles allele flipping** and strand orientation issues
- **Creates optimized files** for fast fine-mapping computation
- **Supports multiple formats** including PLINK and VCF genotype files
- **Enables parallel processing** for computational efficiency

## When to Use

Use `credtools prepare` when you have:

- Chunked summary statistics from the previous step
- Reference genotype data (PLINK .bed/.bim/.fam or VCF files)
- Matched ancestry between summary statistics and reference panels
- Need to create final inputs for credtools fine-mapping

## Basic Usage

### Standard Workflow

```bash
# Create genotype configuration file
echo '{
  "EUR": "/path/to/eur_reference",
  "ASN": "/path/to/asn_reference", 
  "AFR": "/path/to/afr_reference"
}' > genotype_config.json

# Prepare LD matrices
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/
```

### Single Ancestry

```bash
# For single ancestry, still use JSON format
echo '{"EUR": "/path/to/reference"}' > genotype_config.json
credtools prepare chunk_info.txt genotype_config.json prepared/
```

## Command Options

| Option | Description | Default |
|--------|-------------|---------|
| `--threads` / `-t` | Number of threads for parallel processing | 1 |
| `--ld-format` / `-f` | LD computation format (plink/vcf) | plink |
| `--keep-intermediate` / `-k` | Keep intermediate files | False |

## Input Requirements

### Chunk Info File

The `chunk_info.txt` file from the previous chunking step contains:
- Locus coordinates and identifiers
- Paths to chunked summary statistics files
- Ancestry and sample information

### Genotype Configuration

JSON file mapping ancestry codes to genotype file prefixes:

```json
{
  "EUR": "/data/reference/eur_1kg_phase3",
  "ASN": "/data/reference/asn_1kg_phase3",
  "AFR": "/data/reference/afr_1kg_phase3"
}
```

For PLINK format, provide the prefix (without .bed/.bim/.fam extensions).
For VCF format, provide the full path to the VCF file.

### Reference Panel Requirements

**PLINK format** (.bed/.bim/.fam):
- Binary genotype files with complete trio
- BIM file with chromosome, SNP ID, genetic distance, position, alleles
- FAM file with sample information

**VCF format** (planned support):
- Compressed VCF files (.vcf.gz) with tabix index
- Proper chromosome and position formatting
- Consistent allele encoding

## Algorithm Details

### Processing Pipeline

1. **Variant extraction**: Extract variants within each locus region from genotype data
2. **LD computation**: Calculate correlation matrix using PLINK or custom methods
3. **Data intersection**: Match variants between summary statistics and LD data
4. **Allele alignment**: Handle flipped alleles and strand orientation
5. **Quality control**: Filter variants and validate LD matrix properties
6. **File generation**: Create compressed output files for fine-mapping

### Allele Handling

The preparation step carefully handles:
- **Strand flipping**: Automatic detection and correction of strand issues
- **Allele ordering**: Consistent alphabetical ordering of allele pairs
- **Reference matching**: Alignment between summary statistics and reference panel alleles

## Expected Output

### File Structure

For each locus and ancestry combination:

```
prepared/
├── EUR.chr1_12345_67890.sumstats.gz    # Intersected summary statistics
├── EUR.chr1_12345_67890.ld.npz         # Compressed LD matrix
├── EUR.chr1_12345_67890.ldmap.gz       # LD variant mapping
├── ASN.chr1_12345_67890.sumstats.gz
├── ASN.chr1_12345_67890.ld.npz
├── ASN.chr1_12345_67890.ldmap.gz
├── prepared_files.txt                   # Summary of all prepared files
└── final_loci_list.txt                  # Updated loci list for fine-mapping
```

### File Formats

**Sumstats files**: Tab-separated, gzipped summary statistics with only variants present in LD matrix.

**LD matrix files**: NumPy compressed arrays (.npz) containing correlation matrices optimized for memory and speed.

**LD map files**: Tab-separated mapping files linking matrix positions to genomic coordinates and allele information.

**Final loci list**: Updated credtools-compatible format ready for fine-mapping.

## Examples

### Example 1: Standard Multi-Ancestry Preparation

```bash
# Set up genotype configuration
echo '{
  "EUR": "/reference/1000G_phase3/EUR",
  "ASN": "/reference/1000G_phase3/ASN",
  "AFR": "/reference/1000G_phase3/AFR"
}' > genotypes.json

# Prepare with parallel processing
credtools prepare chunk_info.txt genotypes.json prepared/ --threads 4
```

### Example 2: Single Large Reference Panel

```bash
# Use same reference for multiple ancestries
echo '{
  "EUR": "/reference/1000G_phase3/ALL",
  "ASN": "/reference/1000G_phase3/ALL",
  "AFR": "/reference/1000G_phase3/ALL"
}' > genotypes.json

credtools prepare chunk_info.txt genotypes.json prepared/
```

### Example 3: High-Performance Setup

```bash
# Maximize parallel processing and keep intermediate files for debugging
credtools prepare chunk_info.txt genotypes.json prepared/ \
  --threads 8 \
  --keep-intermediate
```

## Genotype Data Setup

### Using 1000 Genomes Project Data

Download and prepare 1000 Genomes reference panels:

```bash
# Download 1000G Phase 3 data
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
tar -xzf 1000GP_Phase3.tgz

# Convert to PLINK format (example for chromosome 1)
plink --vcf 1000GP_Phase3_chr1.vcf.gz --make-bed --out 1000G_chr1

# Combine chromosomes and filter by ancestry
plink --bfile 1000G_merged --keep EUR_samples.txt --make-bed --out EUR_reference
```

### Population-Specific Panels

For best results, use ancestry-matched reference panels:

```json
{
  "EUR": "/reference/UKBB_EUR_50k",
  "ASN": "/reference/BBJ_ASN_10k", 
  "AFR": "/reference/H3Africa_AFR_5k"
}
```

## Performance Optimization

### Parallel Processing

The `--threads` option parallelizes across ancestries:

```bash
# Optimal thread count: number of ancestries or CPU cores, whichever is smaller
credtools prepare chunk_info.txt genotypes.json prepared/ --threads 6
```

### Memory Management

For large datasets:
- Use more threads to parallelize memory usage across cores
- Ensure sufficient RAM (recommended: 4-8GB per thread)
- Consider processing subsets of loci if memory-constrained

### Storage Considerations

Prepared files are optimized for size:
- LD matrices use float16 precision (sufficient for fine-mapping)
- Files are compressed for minimal storage footprint
- Intermediate files can be removed (default behavior)

## Integration with Workflow

Prepared files feed directly into fine-mapping:

```bash
# 1. Munge summary statistics
credtools munge ancestry_files.json munged/

# 2. Identify loci and chunk data
credtools chunk munged/ chunked/

# 3. Prepare LD matrices
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/

# 4. Run fine-mapping (uses final_loci_list.txt)
credtools finemap prepared/final_loci_list.txt results/

# 5. Results are saved in results/ directory
```

## Quality Control

### Checking Preparation Success

```bash
# Count successful preparations
grep -c "created" prepared/prepared_files.txt

# Check for failures
grep "failed\|error" prepared/prepared_files.txt

# Verify file completeness
ls prepared/*.ld.npz | wc -l
```

### Variant Intersection Stats

Review how many variants were successfully matched:

```bash
# Check intersection efficiency per locus
awk '{print $1, $7}' prepared/prepared_files.txt | head -10
```

## Troubleshooting

### Common Issues

**No variants found**: Check that reference genotype files cover the genomic regions of interest and use the correct chromosome encoding (1-22 vs chr1-chr22).

**Allele mismatches**: Verify that summary statistics and reference panels use consistent allele encoding (A/T/G/C vs A/T/C/G vs numeric).

**PLINK errors**: Ensure PLINK is installed and accessible in PATH. Check that genotype files are not corrupted.

**Memory errors**: Reduce the number of threads or process loci in smaller batches.

**File permission errors**: Verify read access to genotype files and write access to output directory.

### Performance Issues

**Slow processing**: Increase `--threads` up to the number of available CPU cores.

**Disk space**: Monitor storage usage, especially with `--keep-intermediate`. Clean up failed runs.

**Network storage**: If using network-mounted genotype files, consider copying to local storage first.

## Advanced Usage

### Custom LD Computation

For specialized reference panels or non-standard formats:

```bash
# Use VCF format (when supported)
credtools prepare chunk_info.txt genotypes.json prepared/ --ld-format vcf
```

### Debugging Failed Loci

```bash
# Keep intermediate files for troubleshooting
credtools prepare chunk_info.txt genotypes.json prepared/ \
  --keep-intermediate

# Check PLINK log files
ls prepared/*_temp.log
```

### Subset Processing

To process only specific loci:

```bash
# Filter chunk_info.txt to specific loci
grep "chr1_" chunk_info.txt > chr1_chunks.txt
credtools prepare chr1_chunks.txt genotypes.json prepared/
```

## Best Practices

1. **Use ancestry-matched references**: Match reference panels to summary statistics ancestry
2. **Verify coordinate systems**: Ensure consistent genome builds (GRCh37/hg19 vs GRCh38/hg38)
3. **Monitor resource usage**: Balance thread count with available memory
4. **Test with small datasets**: Validate setup with a few loci before processing all data
5. **Backup genotype files**: Keep copies of processed reference panels
6. **Document configurations**: Save genotype configurations for reproducibility
7. **Quality control**: Always review preparation summary before fine-mapping

## Tips for Success

- **Start simple**: Begin with single ancestry and small number of loci
- **Check early**: Verify first few loci process correctly before running all
- **Plan storage**: Ensure adequate disk space for output files
- **Use consistent naming**: Keep ancestry codes consistent across all steps
- **Monitor progress**: Watch console output for processing status and errors