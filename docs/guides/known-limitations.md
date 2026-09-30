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
| Direct chunk input | Passing raw file paths to `chunk` skips LD extraction because no `ld_ref` is available. | Use a population config with `ld_ref`, run `credtools prepare` with a genotype config, or provide pre-generated LD files yourself. |
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

## Native convergence metadata

The multi-ancestry wrappers expose a Boolean in their native `CredibleSet` result:

| Tool | Source of `converged` |
| --- | --- |
| Python SuSiE | Finite objective and `0 <= new - old < tol`. |
| MultiSuSiE | The fit's `converged` field: guarded ELBO criterion, after the initial effect-zeroing period. The wrapper uses float64 inference. |
| SuSiEx | The native `.cs` file: `FAIL` means False; `NULL` or a valid, nonempty CS table means True. A missing, empty or malformed file raises an error. |
| MESuSiE | Actual adapter stopping plus two finite computed ELBO values satisfying `0 <= new - old < tol`. The requested tolerance is honored. |

MESuSiE's initial `-Inf` and unused trailing `NA` entries are not computed
observations. Fewer than two computed values, internal missing values,
NaN/infinite values or a nonfinite difference produce False rather than a
default True. From v0.9.7, negative increments no longer count as convergence
in Python SuSiE, MultiSuSiE or the MESuSiE adapter. A substantial decrease is
warned about; a nonfinite computed objective fails the fit. This matches the
nonnegative ELBO convergence rule in recent susieR, not an absolute-difference
rule. It does not promise that all earlier ELBO steps were monotone.

MESuSiE defaults to covariance EM and retains the original Gaussian model.
The development EM adapter evaluates the marginal likelihood with Cholesky
log determinants and computes multivariate posterior moments by latent
square-root conditioning, avoiding direct determinants and subtractive
posterior covariance cancellation. It rejects substantive negative prior
eigenvalues; only eigensolver-scale negative roundoff may be clipped to zero.
This does not project or repair LD. The `native` optimizer mode retains the
installed package's arithmetic as an explicit sensitivity control.
The R wrapper writes `mesusie_status.tsv` and `mesusie_elbo.tsv`
in its temporary directory; iteration counts and status are carried in the
native `CredibleSet` and its `mesusie_runtime` parameters. Inner EM updates
that exhaust the budget are counted separately and do not imply inner
convergence. The validated five-population INS comparison used fixed L=5;
order stability is not a claim about the true number of causal signals,
global optimality or validity under LD mismatch.

The development EM kernel caches data invariant within each single-effect
update and streams configuration-level sufficient moments using log-scaled
weights. It reuses a Cholesky factor for the small covariance solves and log
determinant; dimensions above five use generic Armadillo operations. No SNPs
or configurations are pruned, and float64, the inner iteration/tolerance
settings, eigenvalue floor, independent likelihood guard (absolute tolerance
1e-6) and outer stopping policy
are unchanged. Algebraic equivalence does not imply bitwise equality: tests
compare covariance/likelihood against the frozen v0.9.7 kernel, and application
validation must also check PIPs, CS membership and convergence. Compilation
overhead can dominate tiny fits, so compute-only benchmarks must not be read
as cold-start wall-time guarantees.

For cross-host reproducibility on heterogeneous OpenBLAS CPUs, controlled
comparisons use one thread and a common supported `OPENBLAS_CORETYPE` (for
example, `Nehalem`). Thread count alone does not fix CPU-specific dispatch.
This is a testing/launch setting, not a global environment mutation or an
automatic backend setting imposed on users.

CS/PIP parsing and the existing `empty_on_nonconvergence` behavior are unchanged.
A corrected False may therefore activate that pre-existing setting for MESuSiE;
with the setting disabled, the native PIPs and CS remain available. These are
native-wrapper flags: this change does not add CLI output columns or repair
metadata propagation through later purity filtering/combination.

## First Run Checklist

- Run one locus with `--log-file first_locus.log`.
- Confirm `pips.txt.gz` and `credible_sets_summary.txt.gz` are written.
- Confirm `run_summary.log` has `Failed: 0`.
- Plot one locus with `credtools plot`.
- Only then scale to all loci.
