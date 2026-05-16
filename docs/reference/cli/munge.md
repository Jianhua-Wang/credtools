# `credtools munge`

Clean and standardize GWAS summary statistics.

```bash
credtools munge INPUT_CONFIG OUTPUT_DIR [OPTIONS]
```

## Common Use

```bash
credtools munge population_config.tsv work/munged --force
```

## Inputs

`INPUT_CONFIG` can be:

- a population config TSV with `popu`, `cohort`, `sample_size`, and `path`,
- one summary statistics file path,
- comma-separated summary statistics file paths.

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--config` | JSON column mapping file | none |
| `--force` | overwrite existing output files | off |
| `--interactive` | create a mapping interactively | off |
| `--log-file` | write logs to a file | none |

## Outputs

```text
OUTPUT_DIR/
- {popu}_{cohort}.munged.txt.gz
- sumstat_info_updated.txt
```

`sumstat_info_updated.txt` is the usual input to `credtools chunk`.

## Notes

Munging creates the `SNPID` column and normalizes common columns such as
chromosome, position, alleles, beta, standard error, p-value, and sample size.
