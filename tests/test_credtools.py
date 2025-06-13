#!/usr/bin/env python
"""Tests for `credtools` package."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from credtools import __version__
from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet, load_locus


def test_version():
    """Test that version is a string."""
    assert isinstance(__version__, str)


@pytest.fixture
def sample_locus():
    """Create a sample locus for testing."""
    # Create sample data with all mandatory columns
    sumstats = pd.DataFrame(
        {
            ColName.SNPID: ["rs1", "rs2", "rs3"],
            ColName.CHR: [1, 1, 1],
            ColName.BP: [1000, 2000, 3000],
            ColName.EA: ["A", "C", "G"],  # Effect allele
            ColName.NEA: ["T", "G", "T"],  # Non-effect allele
            ColName.EAF: [0.1, 0.2, 0.3],  # Effect allele frequency
            ColName.A1: ["A", "C", "G"],
            ColName.A2: ["T", "G", "T"],
            ColName.BETA: [0.1, 0.2, 0.3],
            ColName.SE: [0.01, 0.02, 0.03],
            ColName.P: [0.001, 0.002, 0.003],
            ColName.MAF: [0.1, 0.2, 0.3],
        }
    )

    # Create sample LD matrix
    ld_matrix = np.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.3], [0.2, 0.3, 1.0]])

    # Create LDMatrix object
    ld = LDMatrix(sumstats, ld_matrix)

    return Locus(popu="EUR", cohort="test", sample_size=1000, sumstats=sumstats, ld=ld)


def test_locus_creation(sample_locus):
    """Test that a Locus object can be created."""
    assert isinstance(sample_locus, Locus)
    assert len(sample_locus.sumstats) == 3
    assert sample_locus.ld.r.shape == (3, 3)


def test_credible_set_creation(sample_locus):
    """Test that a CredibleSet can be created."""
    # Create a sample credible set
    cs = CredibleSet(
        tool=Method.SUSIE,
        parameters={"max_causal": 1, "coverage": 0.95},
        coverage=0.95,
        n_cs=1,
        cs_sizes=[2],
        lead_snps=["rs1"],
        snps=[["rs1", "rs2"]],
        pips=pd.Series({"rs1": 0.8, "rs2": 0.2, "rs3": 0.0}),
    )

    assert isinstance(cs, CredibleSet)
    assert len(cs.snps) == 1
    assert cs.coverage == 0.95


def test_locus_set_creation(sample_locus):
    """Test that a LocusSet can be created."""
    locus_set = LocusSet([sample_locus])
    assert isinstance(locus_set, LocusSet)
    assert locus_set.n_loci == 1


@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory for file operations."""
    return tmp_path


def test_file_io(temp_dir, sample_locus):
    """Test that loci can be saved and loaded."""
    # Save locus
    output_prefix = str(temp_dir / "test_locus")
    sample_locus.sumstats.to_csv(
        f"{output_prefix}.sumstats.gz", sep="\t", index=False, compression="gzip"
    )
    np.savez_compressed(
        f"{output_prefix}.ld.npz", ld=sample_locus.ld.r.astype(np.float16)
    )
    sample_locus.ld.map.to_csv(
        f"{output_prefix}.ldmap.gz", sep="\t", index=False, compression="gzip"
    )

    # Load locus
    loaded_locus = load_locus(
        prefix=output_prefix, popu="EUR", cohort="test", sample_size=1000
    )

    assert isinstance(loaded_locus, Locus)
    assert len(loaded_locus.sumstats) == len(sample_locus.sumstats)
    assert loaded_locus.ld.r.shape == sample_locus.ld.r.shape
