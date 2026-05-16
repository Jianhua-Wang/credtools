# QC Metrics Dictionary

CREDTOOLS QC is meant to answer one practical question: do the summary
statistics and LD reference look consistent enough to trust fine-mapping?

Start with `qc.txt.gz`, then open the detailed files only for suspicious loci.

## `qc.txt.gz`

| Column | Meaning | How to read it |
| --- | --- | --- |
| `popu` | population label | comes from `loci_list.txt` |
| `cohort` | cohort label | comes from `loci_list.txt` |
| `n_snps` | variants used in QC | small values can mean poor sumstats/LD overlap |
| `n_1e-5` | variants with p-value below `1e-5` | quick signal count |
| `n_5e-8` | variants with p-value below `5e-8` | genome-wide significant count |
| `maf_corr` | correlation between sumstats MAF and LD-reference MAF | low values suggest allele-frequency mismatch |
| `lambda_s` | RSS regularization estimate | high values suggest LD or summary-stat mismatch |
| `n_lambda_s_outlier` | count from LD-mismatch and marginal rules | variants flagged by the kriging RSS rules |
| `n_dentist_s_outlier` | count from Dentist-S | variants inconsistent with the lead-SNP LD pattern |
| `n_c1b_outlier` | optional C1b count | appears when C1b/adaptive QC is enabled |

!!! tip "Look for patterns, not one number"
    A single odd metric is a prompt to inspect the locus. A pattern across many
    loci or one cohort is a stronger sign that the input panel, build, or allele
    coding needs attention.

## `expected_z.txt.gz`

This file comes from the kriging RSS check.

| Column | Meaning |
| --- | --- |
| `SNPID` | variant identifier |
| `z` | transformed z-score |
| `condmean` | expected z-score given other variants and LD |
| `condvar` | conditional variance |
| `z_std_diff` | standardized residual between observed and expected z-score |
| `logLR` | log likelihood ratio for allele-switch-like behavior |
| `lambda_s` | locus-level regularization estimate |
| `cohort` | combined `popu_cohort` label |

Default C1 and C2 rules use:

| Rule | Default condition |
| --- | --- |
| C1 LD mismatch | `logLR > 2` and `abs(z) > 2` |
| C2 marginal | `abs(z) < 2`, `abs(z_std_diff) > 3`, and lead-SNP correlation above `0.8` |
| C1b high-z residual | `abs(z_std_diff) > 10` and `abs(z) > 2` |

C1b is only counted when enabled with `--enable-c1b` or through adaptive QC.

## `dentist_s.txt.gz`

Dentist-S asks whether a variant's marginal association is consistent with the
lead variant and LD.

| Column | Meaning |
| --- | --- |
| `SNPID` | variant identifier |
| `t_dentist_s` | Dentist-S test statistic |
| `-log10p_dentist_s` | evidence against consistency |
| `r2` | LD R2 with the lead variant |
| `cohort` | combined `popu_cohort` label |

The default Dentist-S outlier rule is:

```text
-log10p_dentist_s >= 4 and r2 >= 0.6
```

## `compare_maf.txt.gz`

| Column | Meaning |
| --- | --- |
| `SNPID` | variant identifier |
| `MAF_sumstats` | minor allele frequency from summary statistics |
| `MAF_ld` | minor allele frequency from the LD map `AF2` column |
| `cohort` | combined `popu_cohort` label |

If the LD map has no `AF2` column, CREDTOOLS cannot compute this comparison and
the detailed MAF file may be empty.

## Outlier Removal Files

When `--remove-outlier` or `--adaptive-qc` is used, CREDTOOLS writes a
`cleaned/` folder.

| File | Use |
| --- | --- |
| `cleaned/outlier_snps.txt.gz` | variant-level outlier calls |
| `cleaned/outlier_removal_summary.txt.gz` | counts removed per locus and cohort |
| `cleaned/cleaned_loci_info.txt.gz` | input file for downstream fine-mapping on cleaned data |

`--adaptive-qc` implies outlier removal. It first removes baseline C1, C2, and
C3 outliers. If the cleaned locus still has high `lambda_s`, it adds C1b
outliers from the original locus and recomputes QC.

## Heterogeneity Metrics

Meta-analysis writes heterogeneity outputs before combining cohorts.

| Metric | File | Meaning |
| --- | --- | --- |
| LD fourth moment | `ld_4th_moment.txt.gz` | local LD structure around each variant |
| LD decay | `ld_decay.txt.gz` | average LD by distance bin |
| SNP missingness | `snp_missingness.txt.gz` | which cohorts contain each variant |
| Cochran-Q | `cochran_q.txt.gz` | effect-size heterogeneity across cohorts |
| Summary | `heterogeneity.txt.gz` | per-cohort rollup |

Use `heterogeneity.txt.gz` for the first pass. Open detailed files when one
cohort has unusually high missingness, high LD differences, or many
heterogeneous SNPs.

## Practical Review Flow

1. Sort `qc.txt.gz` by `lambda_s`, `n_lambda_s_outlier`, and `maf_corr`.
2. Pick the worst locus for each suspicious cohort.
3. Open `expected_z.txt.gz`, `dentist_s.txt.gz`, and `compare_maf.txt.gz`.
4. Plot the locus with `credtools plot`.
5. Decide whether to remove outliers, change LD reference, or drop the locus.

