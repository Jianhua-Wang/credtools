"""Preprocessing package for credtools.

This package provides functionality to preprocess whole genome summary statistics
and genotype data for multi-ancestry fine-mapping analysis.
"""

from credtools.preprocessing.chunk import chunk_sumstats, identify_independent_loci
from credtools.preprocessing.munge import munge_sumstats
from credtools.preprocessing.prepare import prepare_finemap_inputs

__all__ = [
    "munge_sumstats",
    "chunk_sumstats",
    "identify_independent_loci",
    "prepare_finemap_inputs",
]
