# Multi-Ancestry Fine-Mapping

Multi-ancestry analysis is useful because LD patterns differ across populations.
That can make credible sets smaller, but it also adds choices. This tutorial
keeps the choices simple.

## Start With One Loci List

Rows with the same `locus_id` belong to the same region:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_1	1	50000000	50500000	EUR	UKBB	400000	data/EUR_UKBB_locus_1
locus_1	1	50000000	50500000	AFR	MVP	90000	data/AFR_MVP_locus_1
locus_1	1	50000000	50500000	EAS	BBJ	180000	data/EAS_BBJ_locus_1
```

Each row has its own summary statistics and LD files.

## Option A: Combine Everything First

```bash
credtools pipeline loci_list.txt results_meta_all \
  --meta-method meta_all \
  --tool susie
```

Use this when you want one combined result per locus.

This is a good first multi-ancestry run because the output is simple. The tradeoff
is that population-specific effects can be harder to see.

## Option B: Combine Within Each Population

```bash
credtools pipeline loci_list.txt results_by_population \
  --meta-method meta_by_population \
  --tool susie
```

Use this when you have multiple cohorts per population and want population-level
fine-mapping inputs.

## Option C: Keep Inputs Separate for a Multi-Input Tool

```bash
credtools pipeline loci_list.txt results_multisusie \
  --meta-method no_meta \
  --tool multisusie
```

This lets a multi-input tool use the rows together instead of collapsing them
first.

Other multi-input tools include:

```bash
credtools pipeline loci_list.txt results_susiex \
  --meta-method no_meta \
  --tool susiex
```

## Compare the Runs

For each run, compare:

```text
results_*/{locus_id}/credible_sets_summary.txt.gz
results_*/{locus_id}/causal_variants.txt.gz
results_*/{locus_id}/pips.txt.gz
```

Useful questions:

- Did the lead variant change?
- Did the credible set shrink?
- Are high-PIP variants shared across populations?
- Did QC flag one population more often than the others?

## Watch Heterogeneity

The pipeline computes heterogeneity before meta-analysis and writes summary
files in each locus output directory.

If a locus has strong heterogeneity, do not rush to interpret the combined
result. Run `meta_by_population` or `no_meta` and compare.

!!! tip "A practical default"
    Run `meta_all + susie` first for a quick baseline. Then run
    `no_meta + multisusie` for loci you care about or loci where QC suggests
    population differences.
