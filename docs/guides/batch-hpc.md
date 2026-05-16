# Batch, HPC, and Performance

Fine-mapping is usually limited by LD matrix size, number of loci, and external
tool runtime. Treat a large run like a batch job: make it resumable, log every
step, and test one locus first.

## A Safe Batch Pattern

```bash
credtools pipeline loci_list.txt results \
  --tool susie \
  --meta-method meta_all \
  --max-causal 5 \
  --log-file results/pipeline.log
```

For a first production run:

1. Run one locus.
2. Run ten mixed loci.
3. Run the full list.

This catches schema, environment, and memory problems early.

## Parallelism

| Command | Option | What it parallelizes |
| --- | --- | --- |
| `chunk` | `--threads` | chunking and LD preparation work |
| `meta` | `--threads` | loci in meta-analysis |
| `qc` | `--threads` | loci in QC |
| `finemap` | `--processes` | per-locus fine-mapping tasks |
| `pipeline` | tool-specific options | use command logs to confirm actual tool behavior |
| FINEMAP tool | `--n-threads` | threads inside FINEMAP |

Do not set every knob to the full core count. If you run 16 CREDTOOLS
processes and each FINEMAP process uses 8 threads, the job can oversubscribe the
node.

## Split by Locus

For clusters, split the loci list and submit array jobs:

```bash
mkdir -p batches
split -l 100 -d --additional-suffix=.txt loci_list.txt batches/loci_
```

Make sure each split keeps the header. One simple approach is:

```bash
header=$(head -n 1 loci_list.txt)
for f in batches/loci_*.txt; do
  sed -i.bak "1i\\
$header" "$f"
  rm "${f}.bak"
done
```

Then run one batch per job:

```bash
credtools pipeline batches/loci_00.txt results_00 \
  --tool susie \
  --meta-method meta_all \
  --log-file results_00/pipeline.log
```

## Resume and Skip

`meta` supports `--skip` for completed loci:

```bash
credtools meta loci_list.txt work/meta --skip
```

For other steps, keep output directories separate by run attempt:

```text
results_susie_v1/
results_susie_v2/
results_finemap_timeout60/
```

This makes it easier to compare outputs and avoids mixing partial files from
different settings.

## Memory Rules of Thumb

LD matrices scale with the square of the variant count. A locus with twice as
many variants can need roughly four times the LD memory.

Before large runs:

- keep loci reasonably bounded,
- check `n_snps` in QC output,
- avoid very wide custom regions unless you need them,
- use fewer parallel workers for large loci,
- keep per-locus logs for failed regions.

## Logging

Use `--log-file` on every long run:

```bash
credtools qc loci_list.txt work/qc \
  --threads 8 \
  --log-file work/qc/qc.log
```

For external tools, also inspect per-tool logs in the locus output or temporary
directory when a run fails.

## FINEMAP Timeouts

FINEMAP has a per-locus timeout:

```bash
credtools finemap loci_list.txt results_finemap \
  --tool finemap \
  --timeout-minutes 60 \
  --processes 4
```

If many loci time out, reduce `--max-causal`, reduce the number of parallel
processes, or split large loci into smaller regions after reviewing the biology.

## Output Hygiene

For reproducible batch runs, save:

- the exact `loci_list.txt`,
- the command line,
- the CREDTOOLS version,
- `parameters.json`,
- `run_summary.log` or `overall_run_summary.log`,
- the environment module list or container tag.

