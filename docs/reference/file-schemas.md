# File Schemas

This page is the strict version of the input and output file formats. Use it
when you are creating files outside CREDTOOLS or checking a failed run.

## Population Config

Used by `credtools munge`:

```text
popu	cohort	sample_size	path
EUR	UKBB	400000	/data/EUR.sumstats.gz
```

Used by `credtools chunk` when LD extraction is needed:

```text
popu	cohort	sample_size	path	ld_ref
EUR	UKBB	400000	work/munged/EUR_UKBB.munged.txt.gz	/ref/EUR
```

| Column | Type | Required | Meaning |
| --- | --- | --- | --- |
| `popu` | string | yes | population or ancestry label |
| `cohort` | string | yes | cohort or study label |
| `sample_size` | integer | yes | cohort sample size |
| `path` | path | yes | summary statistics file |
| `ld_ref` | PLINK prefix | for `chunk` LD extraction | path prefix for `.bed/.bim/.fam` |

## Raw Summary Statistics Aliases

`credtools munge` can recognize common raw headers.

| CREDTOOLS column | Common aliases |
| --- | --- |
| `CHR` | `CHROM`, `#CHROM`, `chromosome`, `Chromosome` |
| `BP` | `POS`, `Position`, `position`, `base_pair_location`, `pos` |
| `SNPID` | `SNP`, `MarkerName`, `variant`, `ID` |
| `EA` | `A1`, `effect_allele`, `ALT`, `Allele1` |
| `NEA` | `A2`, `other_allele`, `REF`, `Allele2` |
| `EAF` | `FRQ`, `FREQ`, `frequency`, `Freq1` |
| `MAF` | `MAF` |
| `BETA` | `beta`, `Beta`, `effect`, `Effect` |
| `SE` | `StdErr`, `stderr`, `standard_error` |
| `P` | `PVAL`, `P_BOLT_LMM`, `pvalue`, `P-value`, `p_value` |
| `N` | `n`, `sample_size`, `NMISS` |
| `INFO` | `info`, `imputation_quality` |
| `Z` | `STAT`, `zscore`, `z_score` |
| `RSID` | `rsid`, `rs` |

When aliases are not enough, pass a JSON mapping to `munge`.

## Munged Summary Statistics

`credtools munge` writes:

```text
CHR	BP	SNPID	EA	NEA	EAF	BETA	SE	P	N	RSID
```

Prepared locus files loaded by fine-mapping may also include `MAF`, which is
derived from `EAF` when CREDTOOLS loads and munges the locus file.

| Column | Required for fine-mapping | Notes |
| --- | --- | --- |
| `SNPID` | yes | generated as `chr-bp-sortedAllele1-sortedAllele2` |
| `CHR`, `BP` | yes | chromosome and base-pair position |
| `EA`, `NEA` | yes | effect and non-effect allele |
| `EAF` | strongly recommended | used to derive `MAF` |
| `MAF` | required by FINEMAP | derived during loading when `EAF` exists |
| `BETA`, `SE`, `P` | yes | model inputs |
| `N` | recommended | sample size also comes from `loci_list.txt` |
| `RSID` | optional | carried into reports when present |

## Loci List

Used by `meta`, `qc`, `finemap`, and `pipeline`.

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
locus_1	1	50000000	50500000	EUR	UKBB	400000	data/EUR_UKBB_locus_1
locus_1	1	50000000	50500000	AFR	MVP	90000	data/AFR_MVP_locus_1
```

| Column | Type | Rule |
| --- | --- | --- |
| `locus_id` | string | rows with the same value are analyzed together |
| `chr` | integer | must be the same for all rows in one `locus_id` |
| `start` | integer | must be positive |
| `end` | integer | must be greater than `start` |
| `popu` | string | population label |
| `cohort` | string | cohort label |
| `sample_size` | integer | must be positive |
| `prefix` | path prefix | no file extension |

Each `popu + cohort + locus_id` combination must be unique.

## Files Behind `prefix`

CREDTOOLS checks these names:

| Data | Accepted names |
| --- | --- |
| summary statistics | `{prefix}.sumstat`, `{prefix}.sumstats.gz` |
| LD matrix | `{prefix}.ld`, `{prefix}.ld.npz` |
| LD map | `{prefix}.ldmap`, `{prefix}.ldmap.gz` |

!!! warning "Prefix is not a folder"
    If `prefix` is `data/EUR_locus_1`, CREDTOOLS reads
    `data/EUR_locus_1.sumstats.gz`, not
    `data/EUR_locus_1/sumstats.gz`.

## LD Matrix

Text LD files use lower-triangular rows:

```text
1
0.12	1
-0.03	0.25	1
```

`.npz` LD files store a square matrix. CREDTOOLS loads the first array in the
archive and replaces missing values with zero.

## LD Map

Minimum columns:

```text
CHR	BP	A1	A2
1	50000123	A	G
1	50000456	C	T
```

Optional but useful:

| Column | Meaning |
| --- | --- |
| `SNPID` | if absent, CREDTOOLS creates one |
| `AF2` | allele frequency used by MAF comparison and SuSiEx preparation |

The number of LD map rows must match the number of LD matrix rows.

## Fine-Mapping Outputs

`pips.txt.gz` always includes:

| Column | Meaning |
| --- | --- |
| `SNPID` | variant identifier |
| `PIP` | posterior inclusion probability |
| `CRED` | credible set index; `0` means not assigned |

For one input row, it also includes available summary-statistic columns such as
`CHR`, `BP`, `RSID`, `EA`, `NEA`, `EAF`, `MAF`, `BETA`, `SE`, `P`, and `R2`.

For multiple input rows, study-specific columns are prefixed:

```text
EUR_UKBB_P
EUR_UKBB_R2
AFR_MVP_P
AFR_MVP_R2
```

Other common result files:

| File | Contents |
| --- | --- |
| `credible_sets_summary.txt.gz` | one row per credible set |
| `causal_variants.txt.gz` | variants with `CRED != 0` |
| `parameters.json` | tool, settings, and run metadata |
| `run_summary.log` | success, failure, and parameter summary |

## QC Outputs

Global and per-locus QC summaries use the same schema:

```text
popu	cohort	n_snps	n_1e-5	n_5e-8	maf_corr	lambda_s	n_lambda_s_outlier	n_dentist_s_outlier
```

When C1b is enabled, `n_c1b_outlier` is appended.

Detailed QC files:

| File | Key columns |
| --- | --- |
| `expected_z.txt.gz` | `SNPID`, `z`, `condmean`, `condvar`, `z_std_diff`, `logLR`, `lambda_s`, `cohort` |
| `dentist_s.txt.gz` | `SNPID`, `t_dentist_s`, `-log10p_dentist_s`, `r2`, `cohort` |
| `compare_maf.txt.gz` | `SNPID`, `MAF_sumstats`, `MAF_ld`, `cohort` |
| `cleaned/outlier_snps.txt.gz` | `SNPID`, `C1_ld_mismatch`, `C2_marginal`, `C3_dentist_s`, optional `C1b_high_z_residual` |
| `cleaned/cleaned_loci_info.txt.gz` | loci list pointing to cleaned files |

## Heterogeneity Outputs

`meta` and `pipeline` can write heterogeneity summaries:

| File | Key columns |
| --- | --- |
| `heterogeneity.txt.gz` | `popu`, `cohort`, `ld_4th_moment_mean`, `ld_decay_rate`, `missing_rate`, `cochran_q_median`, `i_squared_median`, `n_het_snps` |
| `ld_4th_moment.txt.gz` | per-variant LD fourth-moment values by cohort |
| `ld_decay.txt.gz` | `distance_kb`, `r2_avg`, `decay_rate`, `cohort` |
| `cochran_q.txt.gz` | `SNPID`, `Q`, `Q_pvalue`, `I_squared`, `k` |
| `snp_missingness.txt.gz` | variant presence or absence by cohort |

