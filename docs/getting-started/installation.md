# Installation

Use Python 3.9, 3.10, 3.11, or 3.12.

## Install CREDTOOLS

=== "pip"

    ```bash
    python -m pip install credtools
    ```

=== "uv"

    ```bash
    uv pip install credtools
    ```

=== "From source"

    ```bash
    git clone https://github.com/Jianhua-Wang/credtools.git
    cd credtools
    python -m pip install -e ".[dev]"
    ```

Check that the command is available:

```bash
credtools --version
credtools --help
```

## External Tools

Some fine-mapping engines need extra software. You only need to install the
tools you plan to use.

| Tool | Extra dependency | Notes |
| --- | --- | --- |
| `susie` | none beyond Python dependencies | common default |
| `susie_ash` | R package `susieR >= 0.16.1` | adaptive shrinkage background |
| `susie_inf` | R package `susieR >= 0.16.1` | infinitesimal background |
| `finemap` | FINEMAP binary | binary is bundled in `credtools/bin`, but system installs are also supported |
| `carma` | R package `CARMA` | needs R |
| `mesusie` | R package for MESuSiE | needs R |
| `multisusie` | none beyond Python dependencies | best checked with a small run |
| `susiex` | SuSiEx binary/wrapper | needs method-specific setup |

!!! tip "Start with the default"
    If you only want to verify the workflow, install CREDTOOLS and try
    `--tool susie` first. Add the more specialized tools after the basic path
    works on your machine.

For the full environment checklist, see
[External Tools and Environment](../guides/external-tools.md).

## R Setup

Several tools call R through `Rscript`. Check that it is available:

```bash
Rscript --version
```

Install common R packages:

```r
install.packages(c("susieR", "ashr", "devtools"))
```

For CARMA:

```r
devtools::install_github("ZikunY/CARMA")
```

## Development Setup

For local development:

```bash
git clone https://github.com/Jianhua-Wang/credtools.git
cd credtools
uv sync
uv run pytest
```

If you are using `pip` instead of `uv`:

```bash
python -m pip install -e ".[dev]"
pytest
```

## Common Install Problems

!!! warning "R package is installed but CREDTOOLS cannot find it"
    Make sure the same `Rscript` on your `PATH` can load the package:

    ```bash
    Rscript -e "library(susieR); packageVersion('susieR')"
    ```

!!! warning "External binary is not found"
    Run the command directly from your shell first. If your cluster uses modules,
    load the module before running CREDTOOLS.

!!! warning "Package installs but command is missing"
    Your Python environment may not be active. Run:

    ```bash
    python -m pip show credtools
    python -m pip install --upgrade credtools
    ```
