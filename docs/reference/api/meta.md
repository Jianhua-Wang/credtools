# Meta-Analysis API

`meta_all(inputs, ld_weighting="ess")`, `meta_by_population(...)`,
`meta(inputs, meta_method="meta_all", ld_weighting="ess")` and `meta_loci(...)`
share the [matched geometric ESS/SE contract](../cli/meta.md#matched-ld-weighting).
`meta_lds(inputs, ld_weighting="ess")` returns only the LD half of that contract;
use `meta_all` when preparing paired fine-mapping inputs. The standalone
`meta_sumstats` function remains an LD-independent IVW helper and must not be
used on unmatched cohorts to supply statistics for matched LD.

::: credtools.meta
