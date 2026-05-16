# Plot Gallery

Use plots after QC or fine-mapping, not as a replacement for checking the run
summary. The fastest habit is: read the table, plot the suspicious locus, then
decide whether to rerun.

## Quick Commands

=== "QC summary"

    ```bash
    credtools plot work/qc/qc.txt.gz \
      --type summary \
      --output plots/qc_summary.png
    ```

=== "Fine-mapping locus"

    ```bash
    credtools plot work/results/locus_1 \
      --type locusplot \
      --output plots/locus_1_finemap.png
    ```

=== "Locus QC"

    ```bash
    credtools plot work/qc/locus_1 \
      --type locus_qc \
      --output plots/locus_1_qc.png
    ```

## Auto-Detection

If you omit `--type`, CREDTOOLS tries to infer the plot:

| Input | Auto type |
| --- | --- |
| directory with `pips.txt.gz` | `locusplot` |
| directory without `pips.txt.gz` | `locus_qc` |
| file ending in `qc.txt.gz` or `qc.txt` | `summary` |
| file ending in `compare_maf.txt.gz` | `maf_corr` |

When auto-detection fails, pass `--type` explicitly.

## Available Plot Types

<div class="grid cards" markdown>

-   **`summary`**

    Input: `qc.txt.gz`

    Shows lambda-s, MAF correlation, lambda-s outlier counts, and Dentist-S
    outlier counts in one figure.

-   **`locus_qc`**

    Input: a per-locus QC directory

    Combines locus-level QC panels. Use this for the worst loci from the
    summary table.

-   **`locusplot`**

    Input: a fine-mapping locus directory with `pips.txt.gz`

    Shows p-values, LD-to-lead coloring, PIP-scaled points, and credible-set
    membership.

-   **`lambda_s`**

    Input: `qc.txt.gz`

    Boxplot of `lambda_s` by cohort.

-   **`maf_corr`**

    Input: `qc.txt.gz` or `compare_maf.txt.gz`

    Barplot of MAF correlation by cohort.

-   **`lambda_s_outliers`**

    Input: `qc.txt.gz`

    Count of kriging RSS outliers by cohort.

-   **`dentist_s_outliers`**

    Input: `qc.txt.gz`

    Count of Dentist-S outliers by cohort.

-   **`locus_pvalues`**

    Input: `expected_z.txt.gz`

    P-value style view from z-scores, optionally with credible-set annotations
    when called through Python.

-   **`zscore_qq`**

    Input: `expected_z.txt.gz`

    QQ plot for observed z-score behavior.

-   **`ld_decay`**

    Input: `ld_decay.txt.gz`

    LD decay by distance.

-   **`ld_4th_moment`**

    Input: `ld_4th_moment.txt.gz`

    Per-variant LD fourth-moment comparison.

-   **`snp_missingness`**

    Input: `snp_missingness.txt.gz`

    UpSet-style view of variant presence across cohorts.

</div>

## Size and Format

Use `--width`, `--height`, and `--dpi` for reports:

```bash
credtools plot work/qc/qc.txt.gz \
  --type summary \
  --width 14 \
  --height 10 \
  --dpi 300 \
  --output plots/qc_summary.pdf
```

The output extension controls the format. Common choices are `.png`, `.pdf`,
and `.svg`.

## Locus Plot Requirements

`locusplot` reads `pips.txt.gz` and expects:

| Column pattern | Why it matters |
| --- | --- |
| `BP` | x-axis position |
| `PIP` | point size |
| `CRED` | credible-set outline |
| columns ending with `_P` | cohort-specific p-values |
| matching columns ending with `_R2` | LD color scale |

If a locus plot fails with a missing `_R2` column, the PIP table was not created
from a fully matched summary-statistics and LD input.

## Suggested Review Set

For a large run, do not plot every locus at first. Plot:

- the locus with the highest `lambda_s`,
- the locus with the lowest `maf_corr`,
- the locus with the most Dentist-S outliers,
- two or three high-priority biological loci,
- one clean-looking locus as a baseline.

