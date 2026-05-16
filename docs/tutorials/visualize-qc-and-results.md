# Visualize QC and Results

Plots are a fast way to spot problems. Use them after a pipeline run or after
running `credtools qc`.

## Make a Summary QC Plot

For a pipeline locus directory:

```bash
credtools plot results/locus_1 \
  --type summary \
  --output results/locus_1/qc_summary.png
```

For a standalone QC directory:

```bash
credtools plot qc_results \
  --type summary \
  --output qc_results/qc_summary.png
```

## Plot One Locus

```bash
credtools plot results/locus_1 \
  --type locusplot \
  --output results/locus_1/locusplot.png
```

This plot helps you check whether high-PIP variants line up with the association
signal and LD pattern.

## Common Plot Types

| Plot type | Use it for |
| --- | --- |
| `summary` | quick QC overview |
| `locus_qc` | per-locus QC panels |
| `locusplot` | p-values, PIPs, and LD context |
| `lambda_s` | summary-statistic and LD consistency |
| `maf_corr` | allele-frequency agreement |
| `locus_pvalues` | p-value shape inside a locus |
| `zscore_qq` | z-score behavior |
| `ld_decay` | LD decay by population or cohort |
| `ld_4th_moment` | LD structure checks |
| `snp_missingness` | which studies are missing variants |

## Let CREDTOOLS Guess the Plot

If you do not pass `--type`, CREDTOOLS tries to infer the plot from the input:

```bash
credtools plot results/locus_1 --output results/locus_1/auto_plot.png
```

This is useful for quick checks. Use `--type` when you want reproducible reports.

## Size and Format

```bash
credtools plot results/locus_1 \
  --type locusplot \
  --output results/locus_1/locusplot.pdf \
  --width 14 \
  --height 8 \
  --dpi 300
```

The output extension controls the file type: PNG, PDF, or SVG.

!!! tip "Use plots to choose what to inspect, not as the final evidence"
    When a plot looks odd, open the matching table and check exact values. The
    plot is the map; the table is the record.
