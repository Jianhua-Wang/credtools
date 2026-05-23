# Pre-Generated LD With UKB and 1KG Panels

Use this tutorial when LD matrices have already been generated outside
CREDTOOLS. A common setup is UK Biobank EUR LD for European GWAS and 1000
Genomes LD for non-European populations.

This path skips `credtools prepare`. You provide a `loci_list.txt` whose
`prefix` points to existing summary-statistics, LD-matrix, and LD-map files.

## File Layout

For every population and locus, keep the three files behind one shared prefix:

```text
data/
- EUR_UKB_chr9_21900000_22100000.sumstats.gz
- EUR_UKB_chr9_21900000_22100000.ld.npz
- EUR_UKB_chr9_21900000_22100000.ldmap.gz
- AFR_1KG_chr9_21900000_22100000.sumstats.gz
- AFR_1KG_chr9_21900000_22100000.ld.npz
- AFR_1KG_chr9_21900000_22100000.ldmap.gz
```

The prefix for the first EUR row is:

```text
data/EUR_UKB_chr9_21900000_22100000
```

Do not include `.sumstats.gz`, `.ld.npz`, or `.ldmap.gz` in the `prefix`
column.

## Required Locus List

Create a tab-separated `loci_list.txt`:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
chr9_21900000_22100000	9	21900000	22100000	EUR	UKB	400000	data/EUR_UKB_chr9_21900000_22100000
chr9_21900000_22100000	9	21900000	22100000	AFR	1KG	90000	data/AFR_1KG_chr9_21900000_22100000
```

Use the GWAS sample size in `sample_size`. Do not use the reference-panel sample
size unless the reference panel is also the GWAS cohort.

## Check Variant Alignment

The LD map and LD matrix must use the same variant order. The summary
statistics do not need to include every LD variant, but the overlap should be
high enough for stable QC and fine-mapping.

The LD map should include:

```text
SNPID	CHR	BP	A1	A2	AF2
```

`AF2` is optional, but it enables the MAF comparison in QC. For UKB or 1KG LD,
`AF2` should be the allele frequency from that reference panel.

## Run QC First

```bash
credtools qc loci_list.txt qc_ukb_1kg --threads 4
```

Start with:

```text
qc_ukb_1kg/qc.txt.gz
qc_ukb_1kg/chr9_21900000_22100000/expected_z.txt.gz
qc_ukb_1kg/chr9_21900000_22100000/dentist_s.txt.gz
qc_ukb_1kg/chr9_21900000_22100000/compare_maf.txt.gz
```

High `lambda_s`, many kriging RSS outliers, many Dentist-S outliers, or low
`maf_corr` usually means the GWAS and LD reference disagree. Check genome build,
allele coding, population match, and variant order before fine-mapping.

## Run Fine-Mapping

For one combined meta-analysis:

```bash
credtools pipeline loci_list.txt results_meta_all \
  --meta-method meta_all \
  --tool susie
```

For a multi-input method that should see the UKB and 1KG-backed rows
separately:

```bash
credtools pipeline loci_list.txt results_no_meta \
  --meta-method no_meta \
  --tool multisusie
```

## When to Use `prepare` Instead

Use `credtools prepare` when you have PLINK genotype references but not
pre-generated `{prefix}.ld.npz` and `{prefix}.ldmap.gz` files. In that case,
prepare a genotype config:

```json
{
  "EUR": "/ref/ukb_eur",
  "AFR": "/ref/1kg_afr"
}
```

Then run:

```bash
credtools prepare work/chunks/loci_list.txt genotype_config.json work/prepared
```

After that, use `work/prepared/loci_list.txt` for QC and fine-mapping.
