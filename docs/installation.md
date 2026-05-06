# Installation

CREDTOOLS requires Python 3.9 or higher. The base installation includes all dependencies needed for fine-mapping analysis.

## Basic Installation

To install the base CREDTOOLS package, run this command in your terminal:

```bash
$ pip install credtools
```

This is the preferred method to install CREDTOOLS, as it will always install the most recent stable release.


## Development Installation

For development, you may want to install CREDTOOLS in "editable" mode:

```bash
$ git clone https://github.com/Jianhua-Wang/credtools.git
$ cd credtools
$ pip install -e .
```

Check that CREDTOOLS is installed correctly:

```bash
$ credtools --help
```


## Conda Installation

You can also install CREDTOOLS using conda:

```bash
$ conda install -c conda-forge credtools
```

## Source

The source for CREDTOOLS can be downloaded from the [Github repo][].

You can either clone the public repository:

```bash
$ git clone git://github.com/Jianhua-Wang/credtools
```

Or download the [tarball][]:

```bash
$ curl -OJL https://github.com/Jianhua-Wang/credtools/tarball/master
```

Once you have a copy of the source, you can install it with:

```bash
$ cd credtools
$ pip install -e .
```


## External Tool Dependencies

Some fine-mapping tools in CREDTOOLS require external software to be installed separately. The core tools (SuSiE, ABF, MultiSuSiE, RSparsePro) are pure Python and work out of the box. The following tools need additional setup:

### SuSiEx

SuSiEx requires a pre-compiled binary. See the [SuSiEx documentation](https://github.com/getian107/SuSiEx) for installation instructions.

### MESuSiE

MESuSiE (Multiple Ancestry Sum of Single Effects) is an R-based multi-ancestry fine-mapping tool. It requires **R** and the **MESuSiE R package** to be installed on your system.

#### Option 1: Native R Installation

**Step 1: Install R**

=== "macOS"

    ```bash
    # Using Homebrew
    brew install r
    ```

=== "Ubuntu/Debian"

    ```bash
    sudo apt-get update
    sudo apt-get install r-base r-base-dev
    ```

=== "CentOS/RHEL"

    ```bash
    sudo yum install R R-devel
    ```

**Step 2: Install MESuSiE R package**

MESuSiE depends on `RcppArmadillo`, which requires a Fortran compiler. Make sure `gfortran` is available:

=== "macOS"

    ```bash
    # Install gfortran via gcc
    brew install gcc

    # Create ~/.R/Makevars to tell R where to find gfortran
    mkdir -p ~/.R
    cat > ~/.R/Makevars << 'EOF'
    FC = /opt/homebrew/bin/gfortran
    F77 = /opt/homebrew/bin/gfortran
    FLIBS = -L/opt/homebrew/lib/gcc/current -lgfortran -lquadmath
    EOF
    ```

=== "Ubuntu/Debian"

    ```bash
    sudo apt-get install gfortran
    ```

=== "CentOS/RHEL"

    ```bash
    sudo yum install gcc-gfortran
    ```

Then install MESuSiE from GitHub:

```r
# In R console
install.packages("devtools")
devtools::install_github("borangao/MESuSiE")

# Verify installation
library(MESuSiE)
```

Or from the command line:

```bash
Rscript -e 'install.packages("devtools", repos="https://cran.r-project.org")'
Rscript -e 'devtools::install_github("borangao/MESuSiE")'

# Verify
Rscript -e 'library(MESuSiE); cat("MESuSiE installed successfully\n")'
```

#### Option 2: Conda Environment

If you prefer using conda to manage R and its dependencies:

```bash
# Create a conda environment with R
conda create -n credtools_env python=3.12 r-base r-devtools r-rcpp r-rcpparmadillo -c conda-forge

# Activate the environment
conda activate credtools_env

# Install credtools
pip install credtools

# Install MESuSiE R package
Rscript -e 'devtools::install_github("borangao/MESuSiE")'

# Verify both are available
credtools --version
Rscript -e 'library(MESuSiE); cat("MESuSiE installed successfully\n")'
```

!!! tip "Conda Tips"
    - Make sure to install `r-rcpparmadillo` from conda-forge to avoid Fortran compiler issues.
    - If `devtools::install_github` fails inside conda, try installing additional R packages first:
      ```bash
      conda install -c conda-forge r-rcpp r-rcpparmadillo r-nloptr r-ggplot2 r-cowplot r-ggrepel r-progress r-tidyr
      ```

#### Option 3: HPC / Module System

On HPC clusters with a module system:

```bash
# Load R module (exact name varies by cluster)
module load R/4.3.0
# or
module load r/4.3.0

# Install MESuSiE to a user library
Rscript -e 'devtools::install_github("borangao/MESuSiE")'
```

!!! note "HPC Note"
    On some HPC systems, you may need to load additional modules (e.g., `gcc`, `openblas`) before compiling R packages with C++/Fortran code. Check your cluster's documentation or contact your sysadmin.

#### Verifying MESuSiE Setup

After installation, verify that CREDTOOLS can find MESuSiE:

```bash
# Check Rscript is on PATH
which Rscript

# Check MESuSiE is loadable
Rscript -e 'library(MESuSiE); cat("OK\n")'
```

If either command fails, CREDTOOLS will provide a clear error message with installation instructions when you try to run MESuSiE.

### CARMA

CARMA (CAusal Robust Mapping method with Annotations) is an R-based fine-mapping
tool that explicitly models LD / summary-statistic outliers. It requires **R**
and the **CARMA R package**. The Intel **MKL** library is recommended by the
upstream authors for performance.

```bash
# Verify Rscript exists
which Rscript

# Install CARMA from GitHub (after installing R + devtools as above)
Rscript -e 'devtools::install_github("ZikunY/CARMA")'

# Verify the install
Rscript -e 'library(CARMA); cat("CARMA installed successfully\n")'
```

If you do not have Intel MKL available, CARMA still runs against the default
BLAS/LAPACK shipped with R, but expect longer per-locus runtimes.

CREDTOOLS calls CARMA via `Rscript`. If either Rscript or the CARMA package is
missing, the wrapper raises a `FileNotFoundError` with the install command
above.

### SuSiE-ash

SuSiE-ash is the SuSiE 2.0 model with an adaptive-shrinkage prior on the
"unmappable" background effects. It is exposed by the **susieR** R package
(version ≥ 0.16.1) through the `unmappable_effects = "ash"` argument of
`susie_rss`.

```bash
# Verify Rscript exists
which Rscript

# Install susieR from CRAN (preferred)
Rscript -e 'install.packages("susieR", repos = "https://cloud.r-project.org")'

# Or the development version
Rscript -e 'remotes::install_github("stephenslab/susieR")'

# Verify the install
Rscript -e 'library(susieR); cat(as.character(packageVersion("susieR")), "\n")'
```

CREDTOOLS calls susieR via `Rscript`. If either Rscript or the susieR
package is missing (or the installed susieR is too old to expose
`unmappable_effects`), the wrapper raises a clear `FileNotFoundError` /
`RuntimeError` pointing at this section.

### SuSiE-inf

SuSiE-inf is the SuSiE 2.0 model with a single-Gaussian "infinitesimal" prior
on the unmappable background effects. It uses the same **susieR** R package
(version ≥ 0.16.1) entry point as SuSiE-ash, switched on by
`unmappable_effects = "inf"` in `susie_rss`. Follow the SuSiE-ash install
recipe above — a single susieR installation provides both extensions.

CREDTOOLS calls susieR via `Rscript`. If either Rscript or the susieR
package is missing (or the installed susieR is too old to expose
`unmappable_effects`), the wrapper raises a clear `FileNotFoundError` /
`RuntimeError` pointing at this section.


## Troubleshooting

If you encounter any issues during installation, please check:

1. Python version (3.9+ required)
2. pip is up to date
3. You have write permissions for the installation directory

### MESuSiE-Specific Issues

**`Rscript not found`**: Ensure R is installed and `Rscript` is on your `PATH`. Run `which Rscript` to verify.

**`library 'gfortran' not found` (macOS)**: Install gfortran via `brew install gcc` and create `~/.R/Makevars` as described above.

**`ggrepel` not available for this version of R**: Install an older compatible version:
```r
devtools::install_version("ggrepel", version = "0.9.6", repos = "https://cran.r-project.org")
```

**`RcppArmadillo` compilation errors**: Ensure you have a C++ compiler and Fortran compiler available. On macOS, install Xcode Command Line Tools (`xcode-select --install`) and gcc (`brew install gcc`).


## Links

* [Github repo]: https://github.com/Jianhua-Wang/credtools
* [tarball]: https://github.com/Jianhua-Wang/credtools/tarball/master
