# `credtools meta`

Run meta-analysis as a standalone step.

```bash
credtools meta INPUTS OUTDIR [OPTIONS]
```

## Common Use

```bash
credtools meta work/chunks/loci_list.txt work/meta \
  --meta-method meta_all \
  --threads 4
```

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--threads` | worker count | `1` |
| `--meta-method` | `meta_all`, `meta_by_population`, or `no_meta` | `meta_all` |
| `--calculate-lambda-s` | estimate lambda-s while loading loci | off |
| `--skip` | reuse completed loci from a previous run | off |
| `--log-file` | write logs to a file | none |

## Outputs

```text
OUTDIR/
- loci_info.txt
- heterogeneity.txt.gz
- {locus_id}/
```

Use `OUTDIR/loci_info.txt` as the next input for `qc` or `finemap`.

## Meta Methods

| Method | Behavior |
| --- | --- |
| `meta_all` | combine all rows for each locus |
| `meta_by_population` | combine rows within each population |
| `no_meta` | keep rows separate |
