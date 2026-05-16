# External Tools and Environment

CREDTOOLS is a Python package, but some workflows call external programs. You
only need the tools used by your chosen path.

## What You Need

| Workflow | Extra tools |
| --- | --- |
| `munge` only | none beyond Python dependencies |
| `chunk` with LD extraction | `plink` on `PATH` |
| `finemap --tool susie` | none beyond Python dependencies |
| `finemap --tool finemap` | `finemap` executable |
| `finemap --tool susiex` | `SuSiEx` executable |
| `finemap --tool susie_ash` or `susie_inf` | `Rscript` and `susieR >= 0.16.1` |
| `finemap --tool carma` | `Rscript` and R package `CARMA` |
| `finemap --tool mesusie` | `Rscript` and R package `MESuSiE` |

!!! tip "Install only what you use"
    A basic SuSiE run does not need R, FINEMAP, or SuSiEx. Add specialized
    tools after the core workflow works on your machine.

## Python Environment

CREDTOOLS supports Python `>=3.9,<3.13`.

```bash
python -m pip install credtools
credtools --version
credtools --help
```

For local development:

```bash
uv sync
uv run pytest
```

## PLINK

`credtools chunk` uses PLINK when `ld_ref` is present in the population config.
The reference prefix must point to `.bed`, `.bim`, and `.fam` files:

```text
popu	cohort	sample_size	path	ld_ref
EUR	UKBB	400000	work/munged/EUR_UKBB.munged.txt.gz	/ref/1kg/EUR
```

Check PLINK first:

```bash
plink --version
ls /ref/1kg/EUR.bed /ref/1kg/EUR.bim /ref/1kg/EUR.fam
```

CREDTOOLS currently supports PLINK LD extraction. VCF LD extraction is not
implemented.

## FINEMAP and SuSiEx Executables

CREDTOOLS resolves executable paths in this order:

1. a custom path set in Python with `tool_manager.set_tool_path(...)`,
2. the system `PATH`,
3. the bundled `credtools/bin/<tool>` path, if present.

For CLI usage, the most practical option is to put the executable on `PATH`:

```bash
which finemap
finemap --help

which SuSiEx
SuSiEx --help
```

On a cluster, load the module before running CREDTOOLS:

```bash
module load finemap
module load plink
credtools pipeline loci_list.txt results --tool finemap
```

## R-Based Tools

Check R:

```bash
Rscript --version
```

For SuSiE-ash and SuSiE-inf:

```bash
Rscript -e 'library(susieR); packageVersion("susieR")'
```

The wrapper expects a SuSiE 2.0 capable `susieR` with support for
`unmappable_effects`. If the installed package is too old, upgrade `susieR`.

For CARMA:

```r
install.packages("devtools")
devtools::install_github("ZikunY/CARMA")
```

For MESuSiE:

```r
install.packages("devtools")
devtools::install_github("borangao/MESuSiE")
```

## Temporary Files

Some wrappers write intermediate files under `./tmp/`:

| Tool | Temporary parent directory |
| --- | --- |
| FINEMAP | `./tmp/FINEMAP` |
| SuSiEx | `./tmp/SuSiEx` |
| MESuSiE | `./tmp/MESuSiE` |

Normal runs clean these folders after each tool call. Debug-level logging may
keep temporary directories so you can inspect tool inputs and logs.

## Quick Environment Test

Run a tiny smoke test before a large job:

```bash
credtools pipeline test_loci_list.txt smoke_results \
  --tool susie \
  --max-causal 2 \
  --log-file smoke.log
```

If the smoke test fails, fix the environment or schema first. Do not submit the
full batch and hope the cluster logs are easier to read.

