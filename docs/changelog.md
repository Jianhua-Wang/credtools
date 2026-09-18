# Changelog

This page gives the short release story. For the full raw changelog, see
[`CHANGELOG.md`](https://github.com/Jianhua-Wang/credtools/blob/main/CHANGELOG.md).

## 0.9.5

Added `--ld-weighting ess|se` to `meta` and `pipeline`, with ESS as the default.
Both modes match GWAS and LD within cohorts, recompute matched IVW statistics,
and merge LD with geometric per-SNP normalization. ESS is no longer the legacy
pairwise average, so use fresh output directories. Mode manifests and matching
audits make the change explicit and protect resume from incompatible outputs.

Also fixed meta output metadata when the input TSV uses a different column
order. Adaptive L, float16 LD storage, and the absence of default PSD repair
are unchanged. See [matched LD weighting](reference/cli/meta.md#matched-ld-weighting).

Known limitation: INS validation found population-order sensitivity in MESuSiE
and MultiSuSiE. Keep population order fixed when comparing results; this
release does not resolve that sensitivity or change their fitting algorithms.

## 0.9.4

Fixed the adaptive L (`--adaptive-max-causal`) fallback logic.

When a larger L failed, the loop used to restart from `initial - 1` and discard
earlier successful runs (5 → 10 → 15 fails → 4). It now falls back one step at
a time from the failed L (15 fails → 14 → 13 → …), reuses an earlier valid
result when it reaches that L, never exceeds L = 20, and never accepts a
non-converged result as success. See
[Adaptive L](guides/finemapping-tool-requirements.md#adaptive-l).

## 0.9.3

Documented adaptive L (`--adaptive-max-causal`) in more detail.

The guide now explains how CREDTOOLS changes `max_causal`, when it increases or
decreases L, how non-convergence and purity filtering affect the retry loop, and
which tools support the mechanism.

## 0.9.0

Added SuSiE-inf support through `--tool susie_inf`.

Use it when you want a SuSiE 2.0-style model that can account for a broad
background component in addition to sparse signals.

## 0.8.0

Added SuSiE-ash support through `--tool susie_ash`.

This is useful for loci where a few strong effects sit on top of weaker
background effects.

## 0.7.0

Added CARMA support through `--tool carma`.

CARMA is R-based and includes its own outlier-aware modeling behavior.

## 0.6.x

Improved convergence handling and adaptive `max_causal` behavior.

The main user-facing change: adaptive runs can retry more cleanly when a tool
does not converge.

## 0.5.x

Improved QC, outlier handling, heterogeneity output, and credible set
post-processing.

Highlights:

- C1b and adaptive QC options,
- cleaned loci outputs after outlier removal,
- better per-locus credible set columns,
- heterogeneity summaries from meta-analysis.

## Earlier Releases

Earlier releases added MESuSiE, broad wrapper test coverage, preprocessing
tests, and the current CLI structure.

Use GitHub releases or the root changelog when you need exact patch-level
details.
