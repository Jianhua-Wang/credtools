# Python API Cookbook

The CLI is the main path for routine runs. Use the Python API when you want to
inspect objects, run one step in a notebook, or connect CREDTOOLS to your own
workflow.

## Load a Locus Set

```python
import pandas as pd
from credtools.locus import load_locus_set

loci_info = pd.read_csv("loci_list.txt", sep="\t")
locus_info = loci_info[loci_info["locus_id"] == "locus_1"]

locus_set = load_locus_set(locus_info, if_intersect=True)
print(locus_set)
```

`load_locus_set` expects each row to have `prefix`, `popu`, `cohort`,
`sample_size`, `chr`, `start`, `end`, and `locus_id`.

## Run QC for One Locus

```python
from credtools.qc import locus_qc, locus_qc_summary

metrics = locus_qc(locus_set, out_dir="work/qc/locus_1")
summary = locus_qc_summary(metrics)
print(summary)
```

Open the detailed data frames directly:

```python
expected_z = metrics["expected_z"]
dentist_s = metrics["dentist_s"]
compare_maf = metrics["compare_maf"]
```

## Run Meta-Analysis on One Locus

```python
from credtools.meta import meta

meta_locus_set = meta(locus_set, meta_method="meta_all")
print(meta_locus_set)
```

Common `meta_method` values:

| Method | Use |
| --- | --- |
| `meta_all` | combine all input rows into one meta-analyzed locus |
| `meta_by_population` | combine cohorts within each population |
| `no_meta` | keep rows separate |

## Run Fine-Mapping

```python
from credtools.credtools import fine_map

creds = fine_map(
    locus_set,
    tool="susie",
    max_causal=5,
    coverage=0.95,
)

print(creds.n_cs)
print(creds.pips.sort_values(ascending=False).head())
```

Create the same enhanced PIP table that the CLI writes:

```python
pips = creds.create_enhanced_pips_df(locus_set)
pips.to_csv("pips.txt.gz", sep="\t", index=False, compression="gzip")
```

## Run ABF Without LD

The CLI loader expects LD files for `loci_list.txt` workflows. If you want to
run true no-LD ABF, construct a `Locus` object directly.

```python
import pandas as pd
from credtools.locus import Locus
from credtools.wrappers.abf import run_abf

sumstats = pd.read_csv("single_locus.sumstats.gz", sep="\t")

locus = Locus(
    popu="EUR",
    cohort="UKBB",
    sample_size=400000,
    sumstats=sumstats,
    locus_start=50_000_000,
    locus_end=50_500_000,
    ld=None,
)

creds = run_abf(locus, coverage=0.95)
```

## Plot from Python

```python
from credtools.plot import plot_summary_qc, plot_locusplot

plot_summary_qc("work/qc/qc.txt.gz", output_file="qc_summary.png")
plot_locusplot("work/results/locus_1", output_file="locus_1.png")
```

For notebook work, omit `output_file` and display the returned Matplotlib
figure.

## Set an External Tool Path

The CLI uses executables on `PATH`. In Python, you can set a tool path directly:

```python
from credtools.utils import tool_manager

tool_manager.set_tool_path("finemap", "/opt/finemap/finemap")
```

This affects the current Python process.

