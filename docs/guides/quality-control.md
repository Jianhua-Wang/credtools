# Quality Control

QC asks a practical question: do the summary statistics and LD look consistent
enough to trust the fine-mapping result?

Run QC before interpreting a locus, especially for multi-ancestry analysis.

## Run QC

```bash
credtools qc loci_list.txt qc_results --threads 4
```

The main summary is:

```text
qc_results/qc.txt.gz
```

For each locus, CREDTOOLS also writes detailed tables under:

```text
qc_results/{locus_id}/
```

## Run QC With Outlier Removal

```bash
credtools qc loci_list.txt qc_cleaned \
  --threads 4 \
  --remove-outlier
```

If outliers are removed, CREDTOOLS writes cleaned files and a new input list:

```text
qc_cleaned/cleaned/cleaned_loci_info.txt.gz
```

You can pass that file to `finemap`:

```bash
credtools finemap qc_cleaned/cleaned/cleaned_loci_info.txt.gz finemap_cleaned \
  --tool susie
```

## What the Checks Mean

| Check | What it can reveal |
| --- | --- |
| Kriging RSS / lambda-s | summary statistics and LD disagree |
| Dentist-S | variants with suspicious association patterns |
| MAF comparison | allele frequency mismatch |
| SNP missingness | variants present in one study but missing in another |
| LD decay and LD moments | unusual LD structure |
| Cochran's Q | heterogeneity across studies |

## Kriging RSS and Dentist-S

Kriging RSS compares observed z-scores with LD-implied expected z-scores. The
main detailed file is:

```text
qc_results/{locus_id}/expected_z.txt.gz
```

Important columns:

| Column | Meaning |
| --- | --- |
| `z` | observed z-score after CREDTOOLS standardization |
| `condmean` | LD-implied expected z-score |
| `condvar` | conditional variance |
| `z_std_diff` | standardized residual |
| `logLR` | evidence for allele-switch-like LD mismatch |
| `lambda_s` | locus-level RSS regularization estimate |

Dentist-S checks whether each variant is consistent with the lead variant and
its LD to that lead variant. Its detailed file is:

```text
qc_results/{locus_id}/dentist_s.txt.gz
```

The default Dentist-S rule counts variants with
`-log10p_dentist_s >= 4` and `r2 >= 0.6`.

Use both checks together. A high `lambda_s` with many kriging RSS outliers often
points to broad LD/reference mismatch. A small number of Dentist-S outliers can
be local allele, build, or imputation problems around a strong lead signal.

## Useful Flags

| Flag | When to touch it |
| --- | --- |
| `--logLR-threshold` | LD mismatch calls are too strict or too loose |
| `--z-threshold` | z-score outlier calls need adjustment |
| `--dentist-pvalue-threshold` | Dentist-S calls need adjustment |
| `--enable-c1b` | you want the extra high-z residual rule |
| `--adaptive-qc` | you want a two-stage cleanup for elevated lambda-s |

Most users should start with defaults.

## How to Use QC in Practice

1. Run the pipeline or standalone QC.
2. Check `qc_run_summary.log` for failures.
3. Inspect `qc.txt.gz`.
4. Plot the locus if a metric looks odd.
5. Use outlier removal only when you understand what was removed.

!!! warning "Do not blindly remove outliers"
    Outlier removal can help with clear mismatches, but it can also remove real
    signal if thresholds are too aggressive. Keep the original run and compare.
