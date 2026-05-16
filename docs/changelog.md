# Changelog

This page gives the short release story. For the full raw changelog, see
[`CHANGELOG.md`](https://github.com/Jianhua-Wang/credtools/blob/main/CHANGELOG.md).

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
