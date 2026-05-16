# `credtools finemap`

Run fine-mapping without the full pipeline.

```bash
credtools finemap INPUTS OUTDIR [OPTIONS]
```

## Common Use

```bash
credtools finemap work/chunks/loci_list.txt work/finemap \
  --tool susie \
  --max-causal 5
```

## Tool Choice

| Tool | Type |
| --- | --- |
| `abf` | single-input |
| `abf_cojo` | single-input |
| `carma` | single-input |
| `finemap` | single-input |
| `rsparsepro` | single-input |
| `susie` | single-input |
| `susie_ash` | single-input |
| `susie_inf` | single-input |
| `multisusie` | multi-input |
| `susiex` | multi-input |
| `mesusie` | multi-input |

## Common Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--tool` | fine-mapping engine | `susie` |
| `--max-causal` | maximum causal variants | `5` |
| `--adaptive-max-causal` | retry with adjusted `max_causal` | off |
| `--set-L-by-cojo` | estimate signal count with COJO | off |
| `--coverage` | credible set coverage | `0.95` |
| `--timeout-minutes` | per-locus FINEMAP timeout | `30` for FINEMAP |
| `--processes` | per-locus worker processes | `1` |
| `--combine-cred` | `union`, `intersection`, or `cluster` | `union` |
| `--combine-pip` | `max`, `min`, `mean`, or `meta` | `max` |
| `--significant-threshold` | p-value threshold for considering a signal | `5e-8` |
| `--jaccard-threshold` | clustering threshold for credible sets | `0.1` |
| `--purity` | minimum credible set purity | `0.0` |
| `--calculate-lambda-s` | estimate lambda-s while loading loci | off |
| `--log-file` | write logs to a file | none |

## SuSiE Options

| Option | Default |
| --- | --- |
| `--max-iter` | `100` |
| `--estimate-residual-variance` | off |
| `--convergence-tol` | `1e-3` |

## COJO Options

Used with `--set-L-by-cojo`.

| Option | Default |
| --- | --- |
| `--p-cutoff` | `5e-8` |
| `--collinear-cutoff` | `0.9` |
| `--window-size` | `10000000` |
| `--maf-cutoff` | `0.01` |
| `--diff-freq-cutoff` | `0.2` |

## Outputs

```text
OUTDIR/
- {locus_id}/
  - pips.txt.gz
  - credible_sets_summary.txt.gz
  - causal_variants.txt.gz
  - parameters.json
- overall_run_summary.log
```
