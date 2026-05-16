# `credtools qc`

Run quality control on locus-level summary statistics and LD matrices.

```bash
credtools qc INPUTS OUTDIR [OPTIONS]
```

## Common Use

```bash
credtools qc work/chunks/loci_list.txt work/qc --threads 4
```

With outlier removal:

```bash
credtools qc work/chunks/loci_list.txt work/qc_cleaned \
  --threads 4 \
  --remove-outlier
```

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--threads` | worker count | `1` |
| `--remove-outlier` | remove outliers and rerun QC on cleaned data | off |
| `--logLR-threshold` | LD mismatch threshold | `2.0` |
| `--z-threshold` | z-score threshold | `2.0` |
| `--z-std-diff-threshold` | marginal SNP outlier threshold | `3.0` |
| `--r-threshold` | correlation threshold with lead SNP | `0.8` |
| `--dentist-pvalue-threshold` | Dentist-S `-log10(p)` threshold | `4.0` |
| `--dentist-r2-threshold` | Dentist-S LD threshold | `0.6` |
| `--enable-c1b` | enable high-z residual rule | off |
| `--c1b-z-std-diff-threshold` | C1b residual threshold | `10.0` |
| `--adaptive-qc` | run two-stage adaptive QC | off |
| `--adaptive-lambda-threshold` | adaptive QC lambda-s trigger | `0.05` |
| `--log-file` | write logs to a file | none |

## Outputs

```text
OUTDIR/
- qc.txt.gz
- qc_run_summary.log
- {locus_id}/
- cleaned/
```

If cleaned files are created, use:

```text
OUTDIR/cleaned/cleaned_loci_info.txt.gz
```

as the next input for `finemap`.
