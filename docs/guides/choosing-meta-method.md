# Choosing a Meta Method

The meta method controls what happens before fine-mapping. It is one of the
most important choices in a multi-study run.

## Quick Choice

| If you want... | Use |
| --- | --- |
| one combined result | `meta_all` |
| one result per population | `meta_by_population` |
| a multi-input tool to see each row | `no_meta` |
| to debug cohort-specific issues | `no_meta` |

## `meta_all`

```bash
credtools pipeline loci_list.txt results \
  --meta-method meta_all \
  --tool susie
```

This combines all rows for a locus, regardless of population or cohort.

Use it when:

- you expect effects to be broadly shared,
- you want a simple first answer,
- your sample sizes are small and you need power.

Be careful when:

- QC shows strong heterogeneity,
- one population has a different lead signal,
- allele frequencies differ sharply and one cohort looks like an outlier.

## `meta_by_population`

```bash
credtools pipeline loci_list.txt results \
  --meta-method meta_by_population \
  --tool susie
```

This combines cohorts within each population, then keeps populations separate.

Use it when:

- you have several cohorts per population,
- you want ancestry-level interpretation,
- you want to compare populations without cohort noise.

## `no_meta`

```bash
credtools pipeline loci_list.txt results \
  --meta-method no_meta \
  --tool multisusie
```

This keeps every row separate.

Use it when:

- you are running `multisusie`, `susiex`, or another multi-input tool,
- you want to inspect cohort-specific behavior,
- you do not trust a combined result yet.

## A Practical Pattern

For important loci, run two passes:

```bash
credtools pipeline loci_list.txt results_meta_all \
  --meta-method meta_all \
  --tool susie

credtools pipeline loci_list.txt results_no_meta \
  --meta-method no_meta \
  --tool multisusie
```

Then compare the lead variants and credible sets.

!!! tip "Do not treat one method as always correct"
    Meta-analysis is a modeling choice. If the result matters, compare methods
    and read the QC tables before making a claim.
