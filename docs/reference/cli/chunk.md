# `credtools chunk`

Find independent loci, split summary statistics by locus, and optionally
extract LD matrices.

```bash
credtools chunk INPUT_CONFIG OUTPUT_DIR [OPTIONS]
```

## Common Use

```bash
credtools chunk work/munged/sumstat_info_updated.txt work/chunks --threads 4
```

## Inputs

Use a population config with:

```text
popu	cohort	sample_size	path	ld_ref
```

`path` should point to munged summary statistics. `ld_ref` should be a PLINK
prefix if you want LD extraction.

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--distance` | distance threshold for independent loci, in bp | `500000` |
| `--pvalue` | p-value threshold for lead variants | `5e-8` |
| `--merge-overlapping` | merge overlapping loci across inputs | on |
| `--use-most-sig` | use the best SNP if no SNP passes the threshold | on |
| `--min-variants` | minimum variants required in a locus | `10` |
| `--custom-chunks` | TSV with `chr`, `start`, `end` | none |
| `--ld-format` | LD reference format, usually `plink` | `plink` |
| `--keep-intermediate` | keep temporary LD extraction files | off |
| `--threads` | worker count | `1` |
| `--log-file` | write logs to a file | none |

## Outputs

```text
OUTPUT_DIR/
- identified_loci.txt
- loci_list.txt
- chunks/
- prepared/
```

`loci_list.txt` is the handoff file for `pipeline`, `meta`, `qc`, and `finemap`.

## Custom Regions

Use custom chunks when you already know the regions:

```text
chr	start	end
1	50000000	50500000
9	21900000	22100000
```

```bash
credtools chunk population_config.tsv work/chunks --custom-chunks regions.tsv
```
