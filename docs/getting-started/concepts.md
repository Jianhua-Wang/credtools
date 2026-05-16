# Core Concepts

This page explains the words used throughout the docs. You do not need to be a
fine-mapping expert to use CREDTOOLS, but these terms make the commands easier
to understand.

## Summary Statistics

Summary statistics are the GWAS results table. Each row is a variant. The
important columns are usually:

| Column | Meaning |
| --- | --- |
| `CHR` | chromosome |
| `BP` | base-pair position |
| `EA` | effect allele |
| `NEA` | non-effect allele |
| `BETA` | effect size |
| `SE` | standard error |
| `P` | p-value |
| `EAF` | effect allele frequency, if available |
| `N` | sample size, if available |

Different tools and cohorts often use different column names. `credtools munge`
turns them into one common format.

## LD Matrix

LD means linkage disequilibrium. In practical terms, it tells CREDTOOLS how
correlated nearby variants are.

CREDTOOLS expects an LD matrix plus a map file:

```text
my_locus.ld.npz
my_locus.ldmap
```

The `.ldmap` file tells CREDTOOLS which variant belongs to each row and column
of the matrix. The order must match.

!!! warning "LD and summary statistics must describe the same variants"
    Fine-mapping can fail or produce bad results when the LD matrix and summary
    statistics use different alleles, positions, or variant ordering. CREDTOOLS
    checks and intersects them, but the cleaner your inputs are, the better.

## Locus

A locus is one genomic region, such as `chr9:21900000-22100000`.

CREDTOOLS uses loci because fine-mapping is local. You do not usually fine-map
the whole genome as one huge matrix. You split the genome into regions, then
analyze each region.

## Locus Set

A locus set is the same locus measured in one or more studies.

For example, the same region may appear in:

| Population | Cohort | Sample size |
| --- | --- | --- |
| EUR | UKBB | 400000 |
| AFR | MVP | 90000 |
| EAS | BBJ | 180000 |

CREDTOOLS can meta-analyze these rows, run them separately, or pass them to a
multi-input fine-mapping tool.

## PIP

PIP means posterior inclusion probability. It is the model's estimate that a
variant is causal.

The value is between 0 and 1:

- `0.90` means strong evidence for that variant.
- `0.10` means weaker but still worth checking.
- `0.00` means little support in the current model.

PIP is not a p-value. It is a fine-mapping probability after considering the
variants in the locus and the LD pattern.

## Credible Set

A credible set is a small group of variants that should contain a causal variant
with a chosen coverage, often 95%.

If a credible set has 95% coverage, the model is saying: "given the data and
assumptions, this set should contain the causal variant with probability 0.95."

Small credible sets are easier to interpret. Large credible sets usually mean
the data or LD structure cannot separate the variants well.

## Meta-Analysis Method

The `--meta-method` flag controls how CREDTOOLS combines studies before
fine-mapping:

| Method | What it does | Use when |
| --- | --- | --- |
| `meta_all` | combine all rows into one analysis input | you want maximum power and expect shared effects |
| `meta_by_population` | combine cohorts within each population | you want population-level results |
| `no_meta` | keep every row separate | you want to preserve each cohort or use multi-input tools |

## Fine-Mapping Tool

The `--tool` flag chooses the statistical engine. A safe first choice is
`susie`. For multi-ancestry joint analysis, look at `multisusie` and `susiex`.

You do not need to pick the perfect tool on day one. Start with SuSiE, inspect
QC, then compare tools if the locus matters.

## QC

QC checks whether the summary statistics and LD make sense together. CREDTOOLS
looks for issues such as:

- allele or sign mismatches,
- outlier variants,
- missing variants across studies,
- unusual LD structure,
- heterogeneity across cohorts.

QC does not prove the result is correct. It tells you where to look before you
trust the result.
