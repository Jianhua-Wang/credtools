# Getting Started

CREDTOOLS is a command-line workflow for fine-mapping. You give it GWAS
summary statistics and LD information. It gives you QC reports, posterior
inclusion probabilities, and credible sets.

The docs are arranged around the question you probably have right now:

## Which Starting Point Matches You?

=== "I have raw GWAS files"

    Start here if you have whole-genome summary statistics and genotype
    reference files such as PLINK `.bed/.bim/.fam`.

    1. Create a population config file.
    2. Run `credtools munge`.
    3. Run `credtools chunk`.
    4. Run `credtools pipeline`.

    [:octicons-arrow-right-24: Raw GWAS to Results](../tutorials/from-raw-gwas-to-results.md)

=== "I already have locus files"

    Start here if each locus already has:

    - `{prefix}.sumstats.gz` or `{prefix}.sumstat`
    - `{prefix}.ld` or `{prefix}.ld.npz`
    - `{prefix}.ldmap` or `{prefix}.ldmap.gz`

    [:octicons-arrow-right-24: Existing Loci List](../tutorials/run-an-existing-loci-list.md)

=== "I only want one command"

    Use `credtools pipeline` once your `loci_list.txt` is ready.

    ```bash
    credtools pipeline loci_list.txt results --tool susie
    ```

    The pipeline runs meta-analysis, QC, and fine-mapping.

## The Mental Model

Think of the workflow as a series of small handoffs:

| Step | What you have | What CREDTOOLS makes |
| --- | --- | --- |
| `munge` | messy GWAS files | clean, standard summary statistics |
| `chunk` | clean genome-wide files + LD reference | locus-level files and a `loci_list.txt` |
| `meta` | one locus across cohorts or ancestries | combined data for fine-mapping |
| `qc` | summary statistics + LD | QC tables and outlier flags |
| `finemap` | one prepared locus set | PIPs and credible sets |
| `pipeline` | `loci_list.txt` | meta-analysis, QC, and fine-mapping in one run |
| `plot` | QC or result files | quick visual checks |

!!! info "A locus is the unit of work"
    CREDTOOLS fine-maps one genomic region at a time. If you start with
    whole-genome data, `chunk` creates those regions for you.

## A Useful First Rule

Run a small example before a big batch. It catches path issues, missing LD files,
and external tool problems before you spend hours on a full run.

```bash
credtools --help
credtools pipeline --help
```

Then continue to [Installation](installation.md) or [Quickstart](quickstart.md).
