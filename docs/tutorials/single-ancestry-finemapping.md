# Single-Ancestry Fine-Mapping

Use this workflow when you have one population or one cohort per locus. It is a
good first run because there are fewer choices to make.

## Prepare the Loci List

A single-ancestry loci list still uses the same format:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_1	1	50000000	50500000	EUR	cohort1	10000	work/EUR_locus_1
```

Each prefix must have summary statistics, LD, and LD map files.

## Run QC First

QC catches many avoidable mistakes before fine-mapping.

```bash
credtools qc loci_list.txt qc_results --threads 4
```

Look at:

```text
qc_results/qc.txt.gz
qc_results/qc_run_summary.log
```

If you see many outliers, run QC with removal:

```bash
credtools qc loci_list.txt qc_cleaned --threads 4 --remove-outlier
```

When outlier removal creates cleaned files, the next input is:

```text
qc_cleaned/cleaned/cleaned_loci_info.txt.gz
```

## Run SuSiE

SuSiE is the safest first tool for most single-ancestry runs.

```bash
credtools finemap loci_list.txt finemap_susie \
  --tool susie \
  --max-causal 5 \
  --coverage 0.95
```

If you used cleaned QC output:

```bash
credtools finemap qc_cleaned/cleaned/cleaned_loci_info.txt.gz finemap_susie_cleaned \
  --tool susie \
  --max-causal 5
```

## Let CREDTOOLS Pick a Smaller `max_causal`

The pipeline enables COJO-based `L` setting by default. For `finemap`, you can
turn it on explicitly:

```bash
credtools finemap loci_list.txt finemap_susie_cojo \
  --tool susie \
  --set-L-by-cojo
```

This is useful when a locus has one or two clear signals and you do not want to
overfit with a large `max_causal`.

## Try ABF When You Do Not Have LD

Most tools require LD. ABF can run without an LD matrix, but the result is less
rich.

```bash
credtools finemap loci_list.txt finemap_abf --tool abf
```

Use this as a fallback, not as the main plan when good LD is available.

## Interpret the Output

The files you will read most often are:

| File | What to check |
| --- | --- |
| `pips.txt.gz` | variants sorted by fine-mapping probability |
| `credible_sets_summary.txt.gz` | one row per credible set |
| `causal_variants.txt.gz` | variants assigned to credible sets |
| `parameters.json` | tool and parameter record |

!!! tip "A clean single-ancestry result is not always a small credible set"
    If variants are in high LD, the tool may be unable to choose one. That is a
    data limitation, not necessarily a failed run.
