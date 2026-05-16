# Run an Existing Loci List

Use this tutorial when you already have locus-level files. This is the shortest
path into CREDTOOLS.

## The Required Files

For each row in your loci list, CREDTOOLS reads files by `prefix`:

```text
{prefix}.sumstats.gz   # or {prefix}.sumstat
{prefix}.ld.npz        # or {prefix}.ld
{prefix}.ldmap         # or {prefix}.ldmap.gz
```

Your `loci_list.txt` should be tab-separated:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_9p21	9	21900000	22100000	EUR	UKBB	400000	data/EUR_UKBB_9p21
locus_9p21	9	21900000	22100000	AFR	MVP	90000	data/AFR_MVP_9p21
```

All rows with the same `locus_id` must use the same `chr`, `start`, and `end`.

!!! warning "Do not include file extensions in prefix"
    Use `data/EUR_UKBB_9p21`, not `data/EUR_UKBB_9p21.sumstats.gz`.

## Validate the Shape

Before a large run, inspect a few rows:

```bash
head loci_list.txt
ls data/EUR_UKBB_9p21.*
```

Then ask CREDTOOLS for command help:

```bash
credtools pipeline --help
```

## Run the Pipeline

```bash
credtools pipeline loci_list.txt results --tool susie
```

This creates one output directory per locus:

```text
results/
- locus_9p21/
  - pips.txt.gz
  - credible_sets_summary.txt.gz
  - causal_variants.txt.gz
  - run_summary.log
- overall_run_summary.log
```

## Choose the Meta Method

Use `--meta-method` to control how rows are combined before fine-mapping:

=== "`meta_all`"

    ```bash
    credtools pipeline loci_list.txt results_meta_all \
      --meta-method meta_all \
      --tool susie
    ```

    Use this when one combined analysis is the goal.

=== "`meta_by_population`"

    ```bash
    credtools pipeline loci_list.txt results_by_population \
      --meta-method meta_by_population \
      --tool susie
    ```

    Use this when you want population-level inputs before fine-mapping.

=== "`no_meta`"

    ```bash
    credtools pipeline loci_list.txt results_no_meta \
      --meta-method no_meta \
      --tool multisusie
    ```

    Use this when a multi-input tool should see the rows separately.

## Check the Main Result

Read the PIP table:

```bash
gzip -cd results/locus_9p21/pips.txt.gz | head
```

Read the credible set summary:

```bash
gzip -cd results/locus_9p21/credible_sets_summary.txt.gz | head
```

If no credible set is reported, check QC and the smallest p-values in the locus.
No credible set can be a valid result when the evidence is weak.
