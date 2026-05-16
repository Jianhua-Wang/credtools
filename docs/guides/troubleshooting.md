# Troubleshooting

Most CREDTOOLS problems fall into a small number of buckets: paths, columns, LD
matching, external tools, or memory.

## Start Here

Check the run logs first:

```text
overall_run_summary.log
{locus_id}/run_summary.log
qc_run_summary.log
```

Then rerun the failing command with a log file:

```bash
credtools pipeline loci_list.txt results --log-file credtools.log
```

## File Not Found

Check the `prefix` in your loci list.

If the row says:

```text
prefix = data/EUR_locus_1
```

CREDTOOLS expects:

```text
data/EUR_locus_1.sumstats.gz
data/EUR_locus_1.ld.npz
data/EUR_locus_1.ldmap
```

Do not include extensions in `prefix`.

## Missing Columns

For raw summary statistics, run `munge` first:

```bash
credtools munge population_config.tsv work/munged --force
```

For loci lists, required columns are:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
```

Make sure the file is tab-separated, not comma-separated.

## LD Matrix and Summary Statistics Do Not Match

Symptoms:

- very few variants remain after intersection,
- QC flags many variants,
- fine-mapping fails on a locus,
- LD map row count does not match the matrix shape.

Check:

- same genome build,
- same chromosome and position system,
- alleles are normalized,
- LD map row order matches the LD matrix,
- summary statistics were munged.

## External Tool Not Found

Check the tool directly:

```bash
Rscript --version
finemap --help
```

For R packages:

```bash
Rscript -e "library(susieR)"
```

If you are on a cluster, load modules before running CREDTOOLS.

## Memory Problems

Large loci create large LD matrices. If a run is killed:

- reduce the locus size with `credtools chunk --distance`,
- run fewer workers,
- split a large job by chromosome or locus,
- avoid plotting huge loci until the analysis works.

Example:

```bash
credtools chunk population_config.tsv work/chunks --distance 250000 --threads 2
```

## No Credible Sets

This can be a valid result. Check:

- the smallest p-value in the locus,
- whether QC flagged problems,
- whether `--significant-threshold` is too strict,
- whether `--max-causal` is too low for a complex locus.

Try a focused rerun:

```bash
credtools pipeline loci_list.txt results_retry \
  --tool susie \
  --adaptive-max-causal
```

## Results Look Different Across Methods

That is normal. Fine-mapping tools make different modeling choices.

Compare:

- lead variant,
- credible set size,
- high-PIP variants,
- QC flags,
- heterogeneity.

If the locus is important, report the sensitivity instead of hiding it.
