# Summary Statistics Munging

The `credtools munge` command standardizes GWAS summary statistics from various formats into a consistent format suitable for fine-mapping analysis. This preprocessing step is essential for ensuring all your data follows the expected column names and data types.

## Overview

Different GWAS studies often use different column names, formats, and conventions for their summary statistics. The munging process:

- **Standardizes column names** to credtools format (CHR, BP, EA, NEA, BETA, SE, P, etc.)
- **Validates data quality** and identifies potential issues
- **Handles multiple file formats** including compressed files
- **Supports multi-ancestry datasets** with consistent processing across populations
- **Creates reusable configurations** for repeated analyses

## When to Use

Use `credtools munge` when you have:

- Raw GWAS summary statistics from different studies or consortiums
- Files with non-standard column names or formats
- Multiple ancestry-specific summary statistics that need consistent formatting
- Data that requires quality control and validation before fine-mapping

## Basic Usage

### Single File

```bash
credtools munge input_gwas.txt output_dir/
```

### Multiple Files

```bash
credtools munge "file1.txt,file2.txt,file3.txt" output_dir/
```

### Multi-Ancestry with Labels

Create a JSON configuration file mapping ancestry labels to file paths:

```json
{
  "EUR": "/path/to/european_gwas.txt",
  "ASN": "/path/to/asian_gwas.txt", 
  "AFR": "/path/to/african_gwas.txt"
}
```

Then run:

```bash
credtools munge ancestry_files.json output_dir/
```

## Command Options

| Option | Description | Default |
|--------|-------------|---------|
| `--config` / `-c` | Configuration file for column mappings | None |
| `--force` / `-f` | Overwrite existing output files | False |
| `--interactive` / `-i` | Create configuration interactively | False |

## Configuration Files

Configuration files specify how to map input columns to credtools standard format. They're JSON files with the following structure:

```json
{
  "EUR": {
    "CHR": "chromosome",
    "BP": "position", 
    "EA": "effect_allele",
    "NEA": "other_allele",
    "BETA": "beta",
    "SE": "se",
    "P": "pvalue",
    "EAF": "eaf"
  },
  "ASN": {
    "CHR": "chr",
    "BP": "pos",
    "EA": "a1", 
    "NEA": "a2",
    "BETA": "effect",
    "SE": "stderr",
    "P": "p"
  }
}
```

### Interactive Configuration

If your files have non-standard column names, use interactive mode to create a configuration:

```bash
credtools munge input_files.json output_dir/ --interactive
```

This will:
1. Examine the headers of your input files
2. Prompt you to map each required column
3. Save the configuration for future use
4. Apply the mapping to munge your files

## Expected Output

The munging process creates:

- **Munged files**: `{ancestry}.munged.txt.gz` - Standardized summary statistics
- **Validation report**: Console output showing number of variants and validation status per file
- **Configuration file**: If using `--interactive`, saves mapping for reuse

### Standard Column Format

Munged files will have these columns:

| Column | Description | Required |
|--------|-------------|----------|
| CHR | Chromosome | Yes |
| BP | Base pair position | Yes |
| SNPID | SNP identifier (auto-generated if missing) | Yes |
| EA | Effect allele | Yes |
| NEA | Non-effect allele | Yes |
| BETA | Effect size | Yes |
| SE | Standard error | Yes |
| P | P-value | Yes |
| EAF | Effect allele frequency | No |
| MAF | Minor allele frequency | No |

## Examples

### Example 1: Simple Single File

```bash
# Munge a single GWAS file
credtools munge gwas_height.txt output/

# Output:
# output/gwas_height.munged.txt.gz
```

### Example 2: Multi-Ancestry Study

```bash
# Create file mapping
echo '{
  "EUR": "gwas_eur.txt",
  "ASN": "gwas_asn.txt", 
  "AFR": "gwas_afr.txt"
}' > ancestry_files.json

# Munge all files
credtools munge ancestry_files.json output/

# Output:
# output/EUR.munged.txt.gz
# output/ASN.munged.txt.gz  
# output/AFR.munged.txt.gz
```

### Example 3: Custom Column Mapping

```bash
# Create configuration for non-standard columns
echo '{
  "study1": {
    "CHR": "chromosome_name",
    "BP": "genomic_position",
    "EA": "effect_allele", 
    "NEA": "reference_allele",
    "BETA": "effect_size",
    "SE": "standard_error",
    "P": "p_value"
  }
}' > column_config.json

# Munge with custom mapping
credtools munge study1_data.txt output/ --config column_config.json
```

### Example 4: Interactive Setup

```bash
# Let credtools help create the configuration
credtools munge unknown_format.txt output/ --interactive

# This will:
# 1. Show you the column headers found
# 2. Ask you to map each required column
# 3. Save the configuration
# 4. Apply the mapping
```

## Integration with Workflow

Munged files are ready for the next step in the credtools workflow:

```bash
# 1. Munge summary statistics
credtools munge ancestry_files.json munged/

# 2. Identify loci and chunk data  
credtools chunk munged/ chunked/ 

# 3. Prepare LD matrices
credtools prepare chunked/chunk_info.txt genotype_config.json prepared/

# 4. Run fine-mapping
credtools finemap prepared/final_loci_list.txt results/
```

## Troubleshooting

### Common Issues

**Missing required columns**: Ensure your input files contain at minimum: chromosome, position, alleles, effect size, standard error, and p-value.

**File not found errors**: Check file paths in JSON configuration files use absolute paths or paths relative to where you run the command.

**Memory issues**: For very large files, ensure sufficient RAM. Consider splitting large files by chromosome first.

**Inconsistent formats**: Use `--interactive` mode to examine and map columns if automatic detection fails.

### Validation Failures

The munge command validates output files and reports:
- ✓ Files that passed validation
- ✗ Files with missing required columns
- Number of variants per file

Check the console output to identify and fix any validation issues before proceeding to the next step.

## Tips for Success

1. **Keep original files**: Always work on copies of your original data
2. **Use meaningful ancestry labels**: Choose clear names like "EUR", "ASN", "AFR" rather than "study1", "study2"
3. **Save configurations**: Keep column mapping configurations for reuse with similar datasets
4. **Validate results**: Always check the validation output before proceeding
5. **Document your workflow**: Note which configuration was used for reproducibility