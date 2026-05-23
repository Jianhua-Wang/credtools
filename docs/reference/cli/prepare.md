# `credtools prepare`

Build LD-backed locus files from chunked summary statistics and genotype
references.

```bash
credtools prepare INPUTS GENOTYPE_CONFIG OUTPUT_DIR [OPTIONS]
```

Use `prepare` when `chunk` created locus-sized summary-statistics files without
LD files, or when you want to rebuild LD from a different reference panel.

## Common Use

```bash
credtools prepare work/chunks/loci_list.txt genotype_config.json work/prepared \
  --threads 4 \
  --ld-format plink
```

## Inputs

`INPUTS` can be a standard CREDTOOLS loci list:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_1	1	50000000	50500000	EUR	UKB	400000	work/chunks/EUR.locus_1
```

It can also be the internal `chunk_info.txt` shape with `ancestry` and
`sumstats_file`; `prepare` normalizes those columns before processing.

`GENOTYPE_CONFIG` maps each `popu` value to a genotype prefix. JSON is the
simplest form:

```json
{
  "EUR": "/ref/ukb_eur",
  "AFR": "/ref/1kg_afr",
  "EAS": "/ref/1kg_eas"
}
```

For PLINK references, each prefix must point to `.bed`, `.bim`, and `.fam`
files.

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--threads` | worker count | `1` |
| `--ld-format` | genotype reference format | `plink` |
| `--keep-intermediate` | keep temporary PLINK files | off |
| `--log-file` | write logs to a file | none |

!!! warning "VCF extraction is not implemented"
    The CLI accepts `--ld-format vcf`, but the current implementation only has a
    working PLINK extraction path.

## Outputs

```text
OUTPUT_DIR/
- prepared_files.txt
- loci_list.txt
- {popu}.{locus_id}.sumstats.gz
- {popu}.{locus_id}.ld.npz
- {popu}.{locus_id}.ldmap.gz
```

Use `OUTPUT_DIR/loci_list.txt` as the input to `qc`, `meta`, `finemap`, or
`pipeline`.

## Relationship to `chunk`

If your population config already has `ld_ref`, `chunk` runs this preparation
step internally and writes `work/chunks/loci_list.txt`. Running `prepare`
separately is useful when:

- you passed direct summary-statistic paths to `chunk`,
- you changed reference panels,
- you want a separate, inspectable preparation step before QC and fine-mapping.
