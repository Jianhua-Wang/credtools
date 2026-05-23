# Input Files

CREDTOOLS has three main entry points:

1. a population config for raw genome-wide files,
2. a loci list for prepared locus-level files.
3. a genotype config when you run `credtools prepare` separately.

Use the population config when you want CREDTOOLS to create loci for you. Use
the loci list when your locus files already exist. Use the genotype config only
when `credtools prepare` needs to extract LD from PLINK references.

## Population Config

Use this with `credtools munge` and `credtools chunk`.

```text
popu	cohort	sample_size	path	ld_ref
EUR	UKBB	400000	/data/eur.sumstats	/ref/EUR
AFR	MVP	90000	/data/afr.sumstats	/ref/AFR
EAS	BBJ	180000	/data/eas.sumstats	/ref/EAS
```

| Column | Required for | Meaning |
| --- | --- | --- |
| `popu` | `munge`, `chunk` | population or ancestry label |
| `cohort` | `munge`, `chunk` | cohort or study label |
| `sample_size` | `munge`, `chunk` | sample size for the row |
| `path` | `munge`, `chunk` | summary statistics path |
| `ld_ref` | `chunk` with LD extraction | PLINK prefix for `.bed/.bim/.fam` |

`credtools munge` accepts extra columns, so it is fine to include `ld_ref` from
the beginning.

!!! tip "Use stable labels"
    Keep `popu` short (`EUR`, `AFR`, `EAS`) and keep `cohort` readable
    (`UKBB`, `MVP`, `BBJ`). These labels appear in output files.

## Summary Statistics Columns

After munging, CREDTOOLS expects standard columns:

```text
SNPID	CHR	BP	RSID	EA	NEA	EAF	MAF	BETA	SE	P	N
```

The raw input can use common aliases such as `CHROM`, `POS`, `A1`, `A2`,
`PVAL`, or `p_value`. `credtools munge` tries to map them automatically.

If your headers are unusual, create a mapping:

```json
{
  "column_mapping": {
    "chromosome_name": "CHR",
    "genomic_position": "BP",
    "tested_allele": "EA",
    "other_allele": "NEA",
    "effect": "BETA",
    "stderr": "SE",
    "p_value": "P"
  }
}
```

Then run:

```bash
credtools munge population_config.tsv work/munged --config column_mapping.json
```

## Loci List

Use this with `credtools meta`, `credtools qc`, `credtools finemap`, and
`credtools pipeline`.

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_1	1	50000000	50500000	EUR	UKBB	400000	data/EUR_UKBB_locus_1
locus_1	1	50000000	50500000	AFR	MVP	90000	data/AFR_MVP_locus_1
```

| Column | Meaning |
| --- | --- |
| `locus_id` | groups rows from the same genomic region |
| `chr` | chromosome |
| `start` | locus start position |
| `end` | locus end position |
| `popu` | population label |
| `cohort` | cohort label |
| `sample_size` | sample size |
| `prefix` | file prefix, without extension |

All rows with the same `locus_id` must have the same `chr`, `start`, and `end`.
Each `popu + cohort + locus_id` combination must be unique.

## Genotype Config for `prepare`

Use a genotype config when chunked summary statistics exist but LD files still
need to be extracted from genotype references.

JSON:

```json
{
  "EUR": "/ref/ukb_eur",
  "AFR": "/ref/1kg_afr",
  "EAS": "/ref/1kg_eas"
}
```

TSV:

```text
popu	ld_ref
EUR	/ref/ukb_eur
AFR	/ref/1kg_afr
EAS	/ref/1kg_eas
```

Each key must match a `popu` value in the input loci list. Each value is a
genotype prefix, usually a PLINK prefix without `.bed`, `.bim`, or `.fam`.

## Files Behind Each Prefix

CREDTOOLS searches for:

| Data | Accepted names |
| --- | --- |
| summary statistics | `{prefix}.sumstat`, `{prefix}.sumstats.gz` |
| LD matrix | `{prefix}.ld`, `{prefix}.ld.npz` |
| LD map | `{prefix}.ldmap`, `{prefix}.ldmap.gz` |

!!! warning "The prefix is not a directory"
    If `prefix` is `data/EUR_locus_1`, CREDTOOLS reads
    `data/EUR_locus_1.sumstats.gz`, not `data/EUR_locus_1/sumstats.gz`.

!!! warning "CLI workflows expect LD files"
    Even if you plan to use `--tool abf`, the current CLI loader expects an LD
    matrix and LD map behind each `prefix`. Use the Python API if you need a
    true no-LD ABF run.

## Quick Input Checklist

- Paths are relative to the directory where you run the command, or absolute.
- Summary statistics and LD map use the same genome build.
- Alleles are A/C/G/T after munging.
- `sample_size` is a positive integer.
- `prefix` does not include an extension.
- The LD matrix row order matches the LD map row order.
