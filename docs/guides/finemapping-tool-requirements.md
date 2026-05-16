# Fine-Mapping Tool Requirements

Use this page when you are deciding whether a prepared locus has enough data for
a specific fine-mapping tool.

## Short Version

| Tool | Input mode | LD needed by method | Extra software | Notes |
| --- | --- | --- | --- | --- |
| `abf` | single input | no | none | ABF assumes one causal signal; CLI still expects LD files when loading from `loci_list.txt`. |
| `abf_cojo` | single input | yes | none | Uses COJO-style conditioning before ABF. |
| `susie` | single input | yes | none | Good default for most first-pass runs. |
| `susie_ash` | single input | yes | `Rscript`, `susieR >= 0.16.1` | SuSiE with adaptive-shrinkage background. |
| `susie_inf` | single input | yes | `Rscript`, `susieR >= 0.16.1` | SuSiE with infinitesimal background. |
| `finemap` | single input | yes | `finemap` executable | Requires `MAF`; has a per-locus timeout. |
| `rsparsepro` | single input | yes | none | Robust sparse model for LD mismatch sensitivity checks. |
| `carma` | single input | yes | `Rscript`, `CARMA` | Includes CARMA outlier modeling. |
| `multisusie` | multi input | yes | none | Joint multi-population SuSiE-style model. |
| `susiex` | multi input | yes | `SuSiEx` executable | Cross-ancestry fine-mapping wrapper. |
| `mesusie` | multi input | yes | `Rscript`, `MESuSiE` | Joint multi-ancestry model with shared and ancestry-specific signals. |

!!! note "Single input versus multi input"
    Single-input tools run once per row in a locus set and CREDTOOLS combines
    the results. Multi-input tools analyze all rows for the same `locus_id`
    together.

## Required Columns

All CLI fine-mapping workflows start from `loci_list.txt`:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
```

Each `prefix` should resolve to:

```text
{prefix}.sumstats.gz
{prefix}.ld.npz
{prefix}.ldmap.gz
```

The summary statistics should contain:

| Column | Used for |
| --- | --- |
| `SNPID` | matching summary statistics to LD |
| `CHR`, `BP` | location and plotting |
| `EA`, `NEA` | allele alignment |
| `EAF` | allele-frequency checks and derived `MAF` |
| `MAF` | required by FINEMAP; usually derived from `EAF` when loaded |
| `BETA`, `SE`, `P` | all fine-mapping models |
| `N` | useful in prepared files; CLI also uses `sample_size` from `loci_list.txt` |
| `RSID` | optional reporting field |

## No-Signal Behavior

Most single-input wrappers check `--significant-threshold` before running the
model. If no variant passes the threshold, CREDTOOLS returns an empty result
instead of forcing a weak credible set.

```bash
credtools finemap loci_list.txt results \
  --tool susie \
  --significant-threshold 5e-8
```

An empty result usually has:

- `n_cs = 0`
- `PIP = 0` for every variant
- no rows in `causal_variants.txt.gz`

If that is not what you want, loosen `--significant-threshold` for exploratory
runs and label the output accordingly.

## Convergence Options

Iterative tools expose `--max-iter` and convergence tolerance options. When
using adaptive search:

```bash
credtools pipeline loci_list.txt results \
  --tool susie \
  --adaptive-max-causal
```

CREDTOOLS starts from the requested `--max-causal`, tries larger values if the
model appears capped, and can retry smaller values when convergence fails.

For SuSiEx, keep `--mult-step` off when using `--adaptive-max-causal`; both
features try to refine the model and can make debugging harder.

## Choosing the First Tool

Start with:

```bash
credtools pipeline loci_list.txt results \
  --tool susie \
  --meta-method meta_all \
  --max-causal 5
```

Then compare a smaller set of important loci with another tool:

```bash
credtools finemap important_loci.txt results_finemap \
  --tool finemap \
  --max-causal 5 \
  --timeout-minutes 30
```

Use comparison runs to answer a specific question, such as whether a top
credible set is stable under a different model.

