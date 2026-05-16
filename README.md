# CREDTOOLS

[![pypi](https://img.shields.io/pypi/v/credtools.svg)](https://pypi.org/project/credtools/)
[![python](https://img.shields.io/pypi/pyversions/credtools.svg)](https://pypi.org/project/credtools/)
[![Build Status](https://github.com/Jianhua-Wang/credtools/actions/workflows/dev.yml/badge.svg)](https://github.com/Jianhua-Wang/credtools/actions/workflows/dev.yml)
[![codecov](https://codecov.io/gh/Jianhua-Wang/credtools/branch/main/graphs/badge.svg)](https://codecov.io/github/Jianhua-Wang/credtools)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

CREDTOOLS is a command-line and Python toolkit for GWAS fine-mapping workflows.
It helps you standardize summary statistics, define loci, prepare LD-backed
locus files, run QC, combine cohorts, and run fine-mapping tools from one
repeatable interface.

Documentation: <https://Jianhua-Wang.github.io/credtools>

## What It Does

- Standardizes raw GWAS summary statistics with `credtools munge`.
- Finds and chunks loci from genome-wide summary statistics with
  `credtools chunk`.
- Prepares locus-level summary statistics, LD matrices, and LD maps.
- Runs meta-analysis with `meta_all`, `meta_by_population`, or `no_meta`.
- Runs QC for LD and summary-statistic consistency.
- Runs fine-mapping with SuSiE, FINEMAP, ABF, CARMA, SuSiEx, MultiSuSiE,
  MESuSiE, and RSparsePro wrappers.
- Writes PIP tables, credible set summaries, causal-variant tables, logs, and
  plots.

## Install

```bash
python -m pip install credtools
credtools --version
```

CREDTOOLS supports Python `>=3.9,<3.13`.

Some workflows need external tools:

| Workflow | Extra dependency |
| --- | --- |
| LD extraction in `chunk` | PLINK |
| `--tool finemap` | FINEMAP executable |
| `--tool susiex` | SuSiEx executable |
| `--tool susie_ash` or `susie_inf` | Rscript and `susieR >= 0.16.1` |
| `--tool carma` | Rscript and `CARMA` |
| `--tool mesusie` | Rscript and `MESuSiE` |

See the [external tools guide](https://Jianhua-Wang.github.io/credtools/guides/external-tools/)
for setup checks.

## Quick Start

Start from a population config:

```text
popu	cohort	sample_size	path	ld_ref
EUR	UKBB	400000	/data/EUR.sumstats.gz	/ref/EUR
AFR	MVP	90000	/data/AFR.sumstats.gz	/ref/AFR
```

Run the workflow:

```bash
credtools munge population_config.tsv work/munged
credtools chunk work/munged/sumstat_info_updated.txt work/chunks
credtools pipeline work/chunks/loci_list.txt work/results \
  --tool susie \
  --meta-method meta_all \
  --max-causal 5
```

Before a full run, inspect `work/chunks/loci_list.txt`. If `chunk` wrote a
placeholder `sample_size` or cohort label, edit it before fine-mapping.

## Main Commands

| Command | Use it for |
| --- | --- |
| `credtools munge` | clean and standardize summary statistics |
| `credtools chunk` | find loci, create chunks, and prepare LD-backed files |
| `credtools meta` | combine cohorts before QC or fine-mapping |
| `credtools qc` | check summary-statistic and LD consistency |
| `credtools finemap` | run fine-mapping only |
| `credtools pipeline` | run meta-analysis, QC, and fine-mapping together |
| `credtools plot` | create QC and fine-mapping plots |

## Key Input

Most downstream commands use a `loci_list.txt`:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_1	1	50000000	50500000	EUR	UKBB	400000	data/EUR_UKBB_locus_1
```

Each `prefix` should point to:

```text
{prefix}.sumstats.gz
{prefix}.ld.npz
{prefix}.ldmap.gz
```

See the [file schema reference](https://Jianhua-Wang.github.io/credtools/reference/file-schemas/)
for exact columns and accepted file names.

## Development

```bash
git clone https://github.com/Jianhua-Wang/credtools.git
cd credtools
uv sync
uv run pytest
uv run mkdocs serve
```

## Links

- Documentation: <https://Jianhua-Wang.github.io/credtools>
- GitHub: <https://github.com/Jianhua-Wang/credtools>
- PyPI: <https://pypi.org/project/credtools/>
- License: MIT
