# Known Limitations and Gotchas

Most CREDTOOLS runs fail for simple reasons: a placeholder in `loci_list.txt`,
an LD file that was never created, or an external tool that is not on `PATH`.
This page lists the current traps so you can check them before starting a long
run.

!!! warning "Read this before a genome-wide run"
    Run the first locus end to end before launching hundreds of loci. It is much
    cheaper to find a schema or environment problem on one locus.

## Current Limits

| Area | What happens now | What to do |
| --- | --- | --- |
| VCF LD extraction | `credtools chunk --ld-format vcf` is accepted by the CLI, but VCF LD extraction is not implemented. | Use PLINK `.bed/.bim/.fam` references and keep `--ld-format plink`. |
| Custom chunks | `--custom-chunks` reads `chr`, `start`, and `end`, then assigns the internal ancestry label `custom`. This may not match your real population labels. | For now, prefer an explicit `loci_list.txt` when you already know the regions. |
| Auto-created sample size | `chunk` can write `sample_size=50000` as a placeholder in `loci_list.txt`. | Replace it with the real cohort sample size before `meta`, `qc`, `finemap`, or `pipeline`. |
| Auto-created cohort label | `chunk` may set `cohort` equal to the ancestry label. | Edit `cohort` if you need study-level labels such as `UKBB`, `MVP`, or `BBJ`. |
| Direct chunk input | Passing raw file paths to `chunk` skips LD extraction because no `ld_ref` is available. | Use a population config with `ld_ref` or prepare LD files yourself. |
| ABF without LD | The ABF method itself can run without LD in Python, but the current CLI locus loader expects `{prefix}.ld` or `{prefix}.ld.npz` and a matching `{prefix}.ldmap`. | Use the Python API for a true no-LD ABF run, or provide LD files for CLI workflows. |
| FINEMAP MAF | FINEMAP requires a `MAF` column after CREDTOOLS loads the locus. | Make sure `EAF` is present so CREDTOOLS can derive `MAF`, or provide `MAF` in prepared inputs. |
| Multi-input tools | `susiex`, `multisusie`, and `mesusie` analyze all rows in a locus together. | Keep rows for the same `locus_id` aligned to the same chromosome, start, and end. |

## The `loci_list.txt` Check

Before running `pipeline`, open the generated loci list:

```bash
head -n 5 work/chunks/loci_list.txt
```

Check these columns first:

| Column | Check |
| --- | --- |
| `sample_size` | not the placeholder `50000` unless that is really correct |
| `popu` | matches the population label you want in output files |
| `cohort` | matches the study or cohort name you want in reports |
| `prefix` | points to files that actually exist |

Use a quick file check:

```bash
prefix=$(awk 'NR==2 {print $8}' work/chunks/loci_list.txt)
ls "${prefix}.sumstats.gz" "${prefix}.ld.npz" "${prefix}.ldmap.gz"
```

If this fails, fix the input files before running the full pipeline.

## Custom Regions

Custom region files use this shape:

```text
chr	start	end
1	1000000	1500000
1	5000000	5600000
```

At the moment, this path is best for region discovery and manual inspection,
not for a fully automatic LD-prepared pipeline. If you already know the regions
and want a reliable production run, make a `loci_list.txt` directly:

```text
locus_id	chr	start	end	popu	cohort	sample_size	prefix
chr1_1000000_1500000	1	1000000	1500000	EUR	UKBB	400000	data/EUR_UKBB_chr1_1000000_1500000
```

That direct loci list is the most explicit handoff into `qc`, `finemap`, and
`pipeline`.

## Empty Results Are Not Always Errors

Most fine-mapping wrappers check whether any variant passes
`--significant-threshold` before doing expensive work. If no variant passes,
the result can be a valid empty credible set:

- `n_cs = 0`
- all PIPs set to zero
- no lead SNPs

This usually means the locus did not pass the significance threshold used for
fine-mapping. It is different from a tool crash.

## First Run Checklist

- Run one locus with `--log-file first_locus.log`.
- Confirm `pips.txt.gz` and `credible_sets_summary.txt.gz` are written.
- Confirm `run_summary.log` has `Failed: 0`.
- Plot one locus with `credtools plot`.
- Only then scale to all loci.

