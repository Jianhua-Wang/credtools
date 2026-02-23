"""Shared fixtures for wrapper tests."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet


def _make_locus(
    popu: str = "EUR",
    cohort: str = "test",
    seed: int = 42,
    sample_size: int = 10000,
    locus_start: int = 1000,
    locus_end: int = 3000,
    n_snps: int = 20,
    p_range: tuple = (1e-10, 1e-9),
    add_maf: bool = True,
    add_af2: bool = False,
    matched: bool = True,
) -> Locus:
    """Create a Locus fixture with configurable p-value range.

    Parameters
    ----------
    p_range : tuple
        (low, high) for uniform p-value generation.
        Use (1e-10, 1e-9) for significant SNPs (all pass 5e-8).
        Use (0.1, 0.9) for non-significant SNPs (none pass 5e-8).
    matched : bool
        If True, LD map SNPIDs match sumstats SNPIDs exactly.
        If False, LD map has a different order.
    """
    rng = np.random.default_rng(seed)

    bps = [locus_start + i * 100 for i in range(n_snps)]
    snpids = [f"1-{bp}-A-G" for bp in bps]

    p_vals = rng.uniform(p_range[0], p_range[1], size=n_snps)
    betas = rng.normal(0, 0.1, size=n_snps)
    ses = np.abs(rng.normal(0.02, 0.005, size=n_snps))
    ses = np.maximum(ses, 1e-6)
    eafs = rng.uniform(0.1, 0.9, size=n_snps)

    sumstats_data = {
        ColName.SNPID: snpids,
        ColName.CHR: [1] * n_snps,
        ColName.BP: bps,
        ColName.RSID: snpids,
        ColName.EA: ["A"] * n_snps,
        ColName.NEA: ["G"] * n_snps,
        ColName.EAF: eafs,
        ColName.BETA: betas,
        ColName.SE: ses,
        ColName.P: p_vals,
    }
    if add_maf:
        sumstats_data[ColName.MAF] = np.minimum(eafs, 1 - eafs)

    sumstats = pd.DataFrame(sumstats_data)

    # Create a positive-definite LD matrix
    A = rng.normal(0, 1, size=(n_snps, n_snps))
    cov = A @ A.T
    d = np.sqrt(np.diag(cov))
    r = cov / np.outer(d, d)
    np.fill_diagonal(r, 1.0)

    ld_map_data = {
        ColName.SNPID: snpids,
        ColName.CHR: [1] * n_snps,
        ColName.BP: bps,
        ColName.A1: ["A"] * n_snps,
        ColName.A2: ["G"] * n_snps,
    }
    if add_af2:
        ld_map_data["AF2"] = rng.uniform(0.1, 0.9, size=n_snps)

    ld_map = pd.DataFrame(ld_map_data)

    if not matched:
        # Shuffle LD map order to make it unmatched
        perm = rng.permutation(n_snps)
        ld_map = ld_map.iloc[perm].reset_index(drop=True)
        r = r[perm, :][:, perm]

    ld = LDMatrix(ld_map, r)
    return Locus(popu, cohort, sample_size, sumstats, locus_start, locus_end, ld=ld)


def _make_locus_set(
    n_pop: int = 2,
    seed: int = 42,
    **kwargs,
) -> LocusSet:
    """Create a LocusSet with n_pop populations sharing the same locus region."""
    pop_names = ["EUR", "EAS", "AFR", "SAS", "AMR"][:n_pop]
    loci = []
    for i, popu in enumerate(pop_names):
        loci.append(
            _make_locus(
                popu=popu,
                cohort=f"cohort{i}",
                seed=seed + i,
                add_af2=True,
                **kwargs,
            )
        )
    return LocusSet(loci)


# --------------- Fixtures ---------------


@pytest.fixture
def locus_significant():
    """Locus with all SNPs passing 5e-8 threshold."""
    return _make_locus(p_range=(1e-12, 1e-9))


@pytest.fixture
def locus_no_significant():
    """Locus with no SNPs passing 5e-8 threshold."""
    return _make_locus(p_range=(0.1, 0.9))


@pytest.fixture
def locus_with_maf():
    """Locus with MAF column present."""
    return _make_locus(add_maf=True, p_range=(1e-12, 1e-9))


@pytest.fixture
def locus_without_maf():
    """Locus without MAF column."""
    return _make_locus(add_maf=False, p_range=(1e-12, 1e-9))


@pytest.fixture
def locus_with_af2():
    """Locus with AF2 column in LD map (for COJO)."""
    return _make_locus(add_af2=True, p_range=(1e-12, 1e-9))


@pytest.fixture
def locus_unmatched():
    """Locus where sumstats and LD map have different SNPID order."""
    return _make_locus(matched=False, p_range=(1e-12, 1e-9))


@pytest.fixture
def locus_set_two_pop():
    """Create LocusSet with 2 populations."""
    return _make_locus_set(n_pop=2, p_range=(1e-12, 1e-9))


@pytest.fixture
def locus_set_three_pop():
    """Create LocusSet with 3 populations."""
    return _make_locus_set(n_pop=3, p_range=(1e-12, 1e-9))
