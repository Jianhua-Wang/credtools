# Choosing a Fine-Mapping Tool

Start simple. Use SuSiE for the first pass, then compare other tools for loci
that matter.

## Quick Choice

| Situation | Good first tool |
| --- | --- |
| standard single-input locus | `susie` |
| no LD matrix available in Python | `abf` |
| many cohorts or populations kept separate | `multisusie` |
| multi-ancestry comparison with separate inputs | `susiex` |
| sensitivity analysis against another method | `finemap`, `carma`, `rsparsepro` |

!!! warning "CLI note"
    The ABF model itself does not need LD, but the current CLI loader still
    expects LD files behind each `prefix` in `loci_list.txt`. Use the Python API
    for true no-LD ABF runs.

## Single-Input Tools

These tools analyze one input at a time. If your `loci_list.txt` has several
rows for the same locus, CREDTOOLS runs the tool per row and combines results.

| Tool | Notes |
| --- | --- |
| `susie` | strong default; use first |
| `abf` | simple and can work without LD |
| `abf_cojo` | ABF with COJO-style conditioning |
| `finemap` | external FINEMAP engine |
| `rsparsepro` | alternative sparse model |
| `carma` | R-based model with outlier handling |
| `susie_ash` | SuSiE with adaptive shrinkage |
| `susie_inf` | SuSiE-inf wrapper |

Example:

```bash
credtools finemap loci_list.txt results_susie --tool susie --max-causal 5
```

## Multi-Input Tools

These tools analyze a full locus set together.

| Tool | Notes |
| --- | --- |
| `multisusie` | good for joint multi-population analysis |
| `susiex` | designed for cross-population fine-mapping |
| `mesusie` | multi-ethnic SuSiE-style method |

Example:

```bash
credtools pipeline loci_list.txt results_multisusie \
  --meta-method no_meta \
  --tool multisusie
```

## Important Parameters

| Parameter | Plain meaning |
| --- | --- |
| `--max-causal` | maximum number of signals to allow |
| `--coverage` | credible set coverage, usually `0.95` |
| `--adaptive-max-causal` | let CREDTOOLS retry with different `max_causal` values |
| `--set-L-by-cojo` | use COJO to estimate the number of signals |
| `--purity` | filter credible sets with weak internal LD |

For the full adaptive-L retry flow, see
[Adaptive L](finemapping-tool-requirements.md#adaptive-l).

## A Safe First Run

```bash
credtools pipeline loci_list.txt results \
  --tool susie \
  --meta-method meta_all \
  --max-causal 5 \
  --coverage 0.95
```

If many loci hit the maximum number of credible sets, try:

```bash
credtools pipeline loci_list.txt results_adaptive \
  --tool susie \
  --adaptive-max-causal
```

## When to Compare Tools

Compare tools when:

- the credible set drives a downstream experiment,
- QC is clean but the credible set is large,
- two populations point to different lead variants,
- the result is sensitive to `--max-causal`.

Keep the comparison small. Run several tools on a few important loci before
expanding to the full genome.

For exact requirements by tool, see
[Fine-Mapping Tool Requirements](finemapping-tool-requirements.md).
