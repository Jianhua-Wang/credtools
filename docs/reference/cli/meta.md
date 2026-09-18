# `credtools meta`

Run meta-analysis as a standalone step.

```bash
credtools meta INPUTS OUTDIR [OPTIONS]
```

## Common Use

```bash
credtools meta work/chunks/loci_list.txt work/meta \
  --meta-method meta_all \
  --threads 4
```

## Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--threads` | worker count | `1` |
| `--meta-method` | `meta_all`, `meta_by_population`, or `no_meta` | `meta_all` |
| `--ld-weighting` | geometric LD weighting: `ess` or `se` | `ess` |
| `--calculate-lambda-s` | estimate lambda-s while loading loci | off |
| `--skip` | reuse completed loci from a previous run | off |
| `--log-file` | write logs to a file | none |

## Outputs

```text
OUTDIR/
- meta_config.json
- loci_info.txt
- heterogeneity.txt.gz
- {locus_id}/
```

Use `OUTDIR/loci_info.txt` as the next input for `qc` or `finemap`.

## Meta Methods

| Method | Behavior |
| --- | --- |
| `meta_all` | combine all rows for each locus |
| `meta_by_population` | combine rows within each population |
| `no_meta` | keep rows separate |

## Matched LD Weighting

Both `ess` (default) and `se` first match GWAS and LD **within each cohort**,
align allele directions, then use the **union** of those matches. A variant
does not have to occur in every cohort. Only finite BETA and positive finite SE
contribute. A cohort with no eligible matches contributes zero; an entirely
empty output fails explicitly.

In both modes, meta BETA, SE and P are recomputed by inverse-variance weighting
from exactly those matched contributions. The option changes **LD weighting**,
not the summary-statistics estimator. Original discovery GWAS inputs are not
overwritten; these are LD-supported fine-mapping inputs and can lose power.

For cohort `k` and SNP `j`, define `w[k,j]` as:

- `ess`: the supplied cohort `sample_size` for an eligible contribution;
- `se`: `1 / SE[k,j]^2` from the input summary statistics;
- zero when that cohort/SNP contribution is unavailable.

With `a[k,j] = sqrt(w[k,j] / sum_k w[k,j])`, merge LD as
`R[i,j] = sum_k a[k,i] * a[k,j] * R_k[i,j]`. There is **no pairwise
renormalization**. Supply effective sample sizes in `sample_size` for ESS
weighting; credtools does not infer case/control ESS from that field.
Supply any desired LDSC-corrected SE upstream; the option does not apply a new
correction. SE weights describe the covariance of IVW meta Z under independent
cohorts and appropriate LD; ESS weights remain an approximation to that covariance.
Neither option corrects sample overlap or poor reference LD.

```bash
# Default: matched geometric ESS LD
credtools meta work/chunks/loci_list.txt work/meta_ess --meta-method meta_all

# Matched SE LD, with the same matched IVW summary statistics
credtools meta work/chunks/loci_list.txt work/meta_se --ld-weighting se

# The same option works within each population
credtools meta work/chunks/loci_list.txt work/popmeta_se \
  --meta-method meta_by_population --ld-weighting se
```

`no_meta` retains its existing intersection behavior and ignores the weighting
choice. A single-cohort population still uses the matched input contract, with
its BETA/SE and LD unchanged apart from allele alignment/order and numeric rounding;
P is computed from BETA/SE.

### Migration and Audits

This changes the default ESS **algorithm** compared with releases through
v0.9.4: matching occurs before aggregation, and geometric normalization replaces
the old pairwise ESS average. Use a **new output directory**, even for default
ESS. Keep the old release and old outputs for exact historical reproduction.

Output prefixes include `.matched_ess` or `.matched_se`. Each output has
`*.meta_audit.tsv.gz` (per-SNP original valid information, matched information,
fraction retained, and contributing-cohort counts) and `*.meta_cohorts.tsv.gz`
(per-cohort valid/matched SNP counts). Information is `1/SE^2` in these audits
even in ESS mode. The scalar output sample size is the sum of input cohort
sample sizes, **not a per-SNP observed sample size**.

`meta_config.json` records the algorithm and weighting. A different mode or a
legacy output without a manifest cannot be reused or overwritten in place;
`--skip` only reuses compatible output. As before, use a new output directory
when changing the source input files; the manifest is not an input-content hash.

Information and summary-statistics arithmetic use float64; LD accumulation uses
float32 and `.ld.npz` storage remains float16. There is no automatic PSD repair.
The existing file loader's NaN-to-zero LD policy is unchanged: matching an LD
map entry does not establish that every pairwise LD cell was observed.
