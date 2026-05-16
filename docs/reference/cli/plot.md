# `credtools plot`

Create plots from QC or fine-mapping outputs.

```bash
credtools plot INPUT_PATH [OPTIONS]
```

## Common Use

```bash
credtools plot results/locus_1 \
  --type summary \
  --output results/locus_1/qc_summary.png
```

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--type` | plot type; auto-detected if omitted | none |
| `--output` | output path, such as PNG, PDF, SVG | none |
| `--width` | figure width in inches | `16` |
| `--height` | figure height in inches | `12` |
| `--dpi` | output DPI | `300` |
| `--include-upset` | include SNP missingness upset plot | on |

## Plot Types

```text
summary
locus_qc
locusplot
lambda_s
maf_corr
lambda_s_outliers
dentist_s_outliers
locus_pvalues
zscore_qq
ld_decay
ld_4th_moment
snp_missingness
```

## Examples

```bash
credtools plot qc_results --type summary --output qc_summary.png
credtools plot results/locus_1 --type locusplot --output locusplot.pdf
credtools plot results/locus_1 --type snp_missingness --output missingness.svg
```
