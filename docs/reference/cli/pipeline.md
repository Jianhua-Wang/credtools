# `credtools pipeline`

Run meta-analysis, QC, and fine-mapping in one command.

```bash
credtools pipeline INPUTS OUTDIR [OPTIONS]
```

## Common Use

```bash
credtools pipeline work/chunks/loci_list.txt work/results \
  --tool susie \
  --meta-method meta_all
```

## Workflow

For each `locus_id`, the pipeline:

1. loads all rows for the locus,
2. computes heterogeneity,
3. runs the selected meta-analysis method,
4. runs QC unless `--skip-qc` is set,
5. runs fine-mapping,
6. writes result tables and logs.

## Main Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--meta-method` | `meta_all`, `meta_by_population`, or `no_meta` | `meta_all` |
| `--skip-qc` | skip the QC phase | off |
| `--tool` | fine-mapping engine | `susie` |
| `--max-causal` | maximum causal variants | `5` |
| `--adaptive-max-causal` | retry with adjusted `max_causal` | off |
| `--set-L-by-cojo` | estimate signal count with COJO | on |
| `--coverage` | credible set coverage | `0.95` |
| `--significant-threshold` | p-value threshold for considering a signal | `5e-8` |
| `--combine-cred` | combine credible sets from separate runs | `union` |
| `--combine-pip` | combine PIPs from separate runs | `max` |
| `--jaccard-threshold` | cluster threshold for credible sets | `0.1` |
| `--calculate-lambda-s` | estimate lambda-s while loading loci | off |
| `--log-file` | write logs to a file | none |

For exact adaptive-L retry rules, see
[Adaptive L](../../guides/finemapping-tool-requirements.md#adaptive-l).

## Selected Tool-Specific Options

| Tool family | Options |
| --- | --- |
| ABF | `--var-prior` |
| FINEMAP | `--n-iter`, `--n-threads` |
| SuSiE | `--max-iter`, `--estimate-residual-variance`, `--purity`, `--convergence-tol` |
| RSparsePro | `--eps`, `--ubound`, `--cthres`, `--eincre`, `--minldthres`, `--maxldthres`, `--varemax`, `--varemin` |
| SuSiEx | `--mult-step`, `--keep-ambig`, `--min-purity`, `--tol` |
| CARMA | `--outlier-switch`, `--effect-size-prior` |
| MULTISUSIE | `--rho`, `--scaled-prior-variance`, `--standardize`, `--pop-spec-standardization`, `--estimate-prior-variance`, `--estimate-prior-method`, `--pop-spec-effect-priors`, `--iter-before-zeroing-effects`, `--prior-tol` |

Use `credtools pipeline --help` for the full current flag list.

## Outputs

```text
OUTDIR/
- overall_run_summary.log
- {locus_id}/
  - pips.txt.gz
  - credible_sets_summary.txt.gz
  - causal_variants.txt.gz
  - parameters.json
  - run_summary.log
  - QC tables
```
