# Raw GWAS to Results

Use this tutorial when your starting point is whole-genome summary statistics.
You will clean the files, split them into loci, build LD inputs, and run
fine-mapping.

## What You Need

For each population or cohort:

- one GWAS summary statistics file,
- sample size,
- population label,
- cohort name,
- LD reference files in PLINK format (`.bed`, `.bim`, `.fam`).

The example below uses the small mock data in `exampledata/test_mock_data`.

## Step 1: Write a Population Config

Create a tab-separated file with one row per study.

```bash
cat > population_config.tsv <<'EOF'
popu	cohort	sample_size	path	ld_ref
EUR	cohort1	10000	exampledata/test_mock_data/EUR_all_loci.sumstats	exampledata/test_mock_data/EUR_all_loci
AFR	cohort1	8000	exampledata/test_mock_data/AFR_all_loci.sumstats	exampledata/test_mock_data/AFR_all_loci
EAS	cohort1	12000	exampledata/test_mock_data/EAS_all_loci.sumstats	exampledata/test_mock_data/EAS_all_loci
EOF
```

The `path` column points to raw summary statistics. The `ld_ref` column is the
PLINK prefix, without `.bed`, `.bim`, or `.fam`.

!!! info "Why this file matters"
    CREDTOOLS carries this metadata forward. Later outputs can still tell which
    row came from EUR, AFR, EAS, or any cohort labels you used.

## Step 2: Munge the Summary Statistics

```bash
credtools munge population_config.tsv work/munged --force
```

Munging does three things:

1. renames common columns into the CREDTOOLS format,
2. removes obvious bad rows,
3. creates a stable `SNPID` from chromosome, position, and alleles.

After this step, check:

```text
work/munged/sumstat_info_updated.txt
```

This updated config points to the munged files and keeps the original metadata.

## Step 3: Chunk the Genome Into Loci

```bash
credtools chunk \
  work/munged/sumstat_info_updated.txt \
  work/chunks \
  --threads 4
```

Chunking finds significant regions, cuts the summary statistics to each region,
and extracts LD matrices from the PLINK reference panels.

The handoff file is:

```text
work/chunks/loci_list.txt
```

Open it once. It should have columns like:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
```

The `prefix` column is important. CREDTOOLS uses it to find:

```text
{prefix}.sumstats.gz
{prefix}.ld.npz
{prefix}.ldmap
```

## Step 4: Run the Full Pipeline

```bash
credtools pipeline \
  work/chunks/loci_list.txt \
  work/results \
  --tool susie \
  --meta-method meta_all
```

For each `locus_id`, the pipeline creates a subdirectory:

```text
work/results/
- locus_1/
- locus_2/
- overall_run_summary.log
```

Inside each locus directory, look for:

```text
pips.txt.gz
credible_sets_summary.txt.gz
causal_variants.txt.gz
parameters.json
run_summary.log
```

## Step 5: Make a Quick Plot

```bash
credtools plot \
  work/results/locus_1 \
  --type summary \
  --output work/results/locus_1/qc_summary.png
```

Use plots as a fast check. Use the tables when you need exact values.

## When Something Fails

Start with `run_summary.log` and `overall_run_summary.log`. Most failures are
path or input-shape problems:

- a `prefix` points to files that do not exist,
- LD and summary statistics do not share enough variants,
- an external tool is not installed,
- a locus is too large for available memory.

See [Troubleshooting](../guides/troubleshooting.md) for the common fixes.
