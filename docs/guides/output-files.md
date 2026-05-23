# Output Files

CREDTOOLS writes a lot of files because each step keeps its handoff files. This
page explains the files you will actually open most often.

## After `munge`

```text
work/munged/
- EUR_UKBB.munged.txt.gz
- AFR_MVP.munged.txt.gz
- sumstat_info_updated.txt
```

| File | Use it for |
| --- | --- |
| `*.munged.txt.gz` | cleaned summary statistics |
| `sumstat_info_updated.txt` | input to `credtools chunk` |

## After `chunk`

```text
work/chunks/
- identified_loci.txt
- loci_list.txt
- sumstat_info_updated.txt
- chunks/
- prepared/
```

| File or directory | Use it for |
| --- | --- |
| `identified_loci.txt` | see which regions were found |
| `chunks/` | locus-sized summary statistics |
| `prepared/` | prepared summary statistics, LD matrices, and LD maps |
| `loci_list.txt` | input to `pipeline`, `meta`, `qc`, or `finemap` |

If `ld_ref` was present, `chunk` also tries to build LD files.

## After `prepare`

```text
work/prepared/
- prepared_files.txt
- loci_list.txt
- EUR.locus_1.sumstats.gz
- EUR.locus_1.ld.npz
- EUR.locus_1.ldmap.gz
```

| File | Use it for |
| --- | --- |
| `prepared_files.txt` | status table for prepared locus file sets |
| `loci_list.txt` | input to `pipeline`, `meta`, `qc`, or `finemap` |
| `*.sumstats.gz` | locus-level munged summary statistics |
| `*.ld.npz` | compressed LD matrix |
| `*.ldmap.gz` | LD variant map and reference allele frequencies |

## After `meta`

```text
work/meta/
- loci_info.txt
- heterogeneity.txt.gz
- {locus_id}/
  - {prefix}.sumstats.gz
  - {prefix}.ld.npz
  - {prefix}.ldmap.gz
```

| File | Use it for |
| --- | --- |
| `loci_info.txt` | input to QC or fine-mapping after meta-analysis |
| `heterogeneity.txt.gz` | cross-study heterogeneity summary |
| `{locus_id}/...` | meta-analyzed locus files |

## After `qc`

```text
work/qc/
- qc.txt.gz
- qc_run_summary.log
- {locus_id}/
- cleaned/
```

Common files:

| File | Use it for |
| --- | --- |
| `qc.txt.gz` | one summary table across loci |
| `{locus_id}/qc.txt.gz` | per-locus QC summary |
| `{locus_id}/expected_z.txt.gz` | kriging RSS LD and z-score consistency |
| `{locus_id}/dentist_s.txt.gz` | Dentist-S outlier statistics |
| `{locus_id}/compare_maf.txt.gz` | allele-frequency checks |
| `cleaned/cleaned_loci_info.txt.gz` | input after `--remove-outlier` |

The `cleaned/` directory appears only when outlier removal is requested and
outliers are found.

## After `pipeline`

```text
work/results/
- overall_run_summary.log
- {locus_id}/
  - pips.txt.gz
  - credible_sets_summary.txt.gz
  - causal_variants.txt.gz
  - parameters.json
  - run_summary.log
  - expected_z.txt.gz
  - dentist_s.txt.gz
  - compare_maf.txt.gz
```

Start with these:

| File | First question it answers |
| --- | --- |
| `overall_run_summary.log` | did all loci finish? |
| `{locus_id}/run_summary.log` | did this locus finish? |
| `{locus_id}/pips.txt.gz` | which variants have high PIP? |
| `{locus_id}/credible_sets_summary.txt.gz` | how many credible sets were found? |
| `{locus_id}/causal_variants.txt.gz` | which variants are in credible sets? |
| `{locus_id}/parameters.json` | which tool and settings were used? |

## Reading PIP Tables

`pips.txt.gz` includes `SNPID`, `PIP`, `CRED`, and available summary-statistic
columns. For multi-input results, study-specific columns are prefixed by
population and cohort.

```bash
gzip -cd work/results/locus_1/pips.txt.gz | head
```

`CRED = 0` means the variant is not assigned to a credible set. `CRED = 1`
means it belongs to the first credible set.

## Reading Credible Set Summaries

```bash
gzip -cd work/results/locus_1/credible_sets_summary.txt.gz | head
```

Use this file for reports. Use `pips.txt.gz` when you need variant-level detail.
