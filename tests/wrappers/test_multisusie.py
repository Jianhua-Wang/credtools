"""Tests for the MultiSuSiE multi-ancestry fine-mapping wrapper."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.wrappers.multisusie import run_multisusie

from .conftest import _make_locus, _make_locus_set


def _make_mock_multisusie_rss(n_variants, cs_indices=None, purity_values=None):
    """Factory to create mock multisusie_rss functions.

    Parameters
    ----------
    n_variants : int
        Number of variants in the union set.
    cs_indices : list of list of int or None
        Credible set variant indices. None for no CS.
    purity_values : list of float or None
        Purity values for each CS.
    """

    def mock_fn(**kwargs):
        pip = np.zeros(n_variants)
        if cs_indices:
            for cs in cs_indices:
                for idx in cs:
                    pip[idx] = max(pip[idx], 0.5 / len(cs))
            # Set higher PIP for first variant in first CS
            if cs_indices[0]:
                pip[cs_indices[0][0]] = 0.8

        cs_list = cs_indices if cs_indices else []
        purity_arr = np.array(purity_values) if purity_values else np.full(
            len(cs_list), np.nan
        )
        # include_mask: True for valid CS, False otherwise
        include_mask = [True] * len(cs_list)

        return SimpleNamespace(
            pip=pip,
            sets=(cs_list, purity_arr, [0.95] * len(cs_list), include_mask),
        )

    return mock_fn


class TestRunMultisusieBasic:
    """Basic MultiSuSiE functionality tests."""

    def test_basic_call_two_pop(self, locus_set_two_pop, monkeypatch):
        """MultiSuSiE returns a CredibleSet for 2-population LocusSet."""
        n = 20  # expected union size (same SNPs across pops)
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(n, cs_indices=[[0, 1, 2]], purity_values=[0.85]),
        )
        result = run_multisusie(locus_set_two_pop)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.MULTISUSIE

    def test_output_structure(self, locus_set_two_pop, monkeypatch):
        """MultiSuSiE result has correct structure."""
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(n, cs_indices=[[0, 1]], purity_values=[0.9]),
        )
        result = run_multisusie(locus_set_two_pop)
        assert result.n_cs >= 0
        assert isinstance(result.pips, pd.Series)
        assert result.coverage == 0.95

    def test_parameters_stored(self, locus_set_two_pop, monkeypatch):
        """Parameters should be stored in the result."""
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(n, cs_indices=[[0]], purity_values=[0.9]),
        )
        result = run_multisusie(
            locus_set_two_pop, max_causal=5, coverage=0.99, rho=0.8
        )
        assert result.parameters["max_causal"] == 5
        assert result.parameters["coverage"] == 0.99
        assert result.parameters["rho"] == 0.8


class TestRunMultisusieVariantUnion:
    """Tests for variant union merging across populations."""

    def test_variant_union_dedup(self, locus_set_two_pop, monkeypatch):
        """Variant union should be deduplicated."""
        # Track what the mock receives to validate z_list/R_list sizes
        captured_kwargs = {}

        def capturing_mock(**kwargs):
            captured_kwargs.update(kwargs)
            n = len(kwargs["z_list"][0])
            pip = np.zeros(n)
            pip[0] = 0.5
            return SimpleNamespace(
                pip=pip,
                sets=([], np.array([]), [], []),
            )

        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss", capturing_mock
        )
        run_multisusie(locus_set_two_pop)
        # z_list and R_list should match the number of populations
        assert len(captured_kwargs["z_list"]) == 2
        assert len(captured_kwargs["R_list"]) == 2

    def test_rho_matrix_construction(self, locus_set_two_pop, monkeypatch):
        """Rho matrix should be correctly constructed."""
        captured_kwargs = {}

        def capturing_mock(**kwargs):
            captured_kwargs.update(kwargs)
            n = len(kwargs["z_list"][0])
            return SimpleNamespace(
                pip=np.zeros(n),
                sets=([], np.array([]), [], []),
            )

        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss", capturing_mock
        )
        run_multisusie(locus_set_two_pop, rho=0.6)
        rho_mat = captured_kwargs["rho"]
        assert rho_mat.shape == (2, 2)
        assert rho_mat[0, 0] == 1.0
        assert rho_mat[1, 1] == 1.0
        assert rho_mat[0, 1] == 0.6
        assert rho_mat[1, 0] == 0.6


class TestRunMultisusieCS:
    """Tests for credible set handling."""

    def test_no_cs(self, locus_set_two_pop, monkeypatch):
        """Should handle no credible sets."""
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(n, cs_indices=[], purity_values=[]),
        )
        result = run_multisusie(locus_set_two_pop)
        assert result.n_cs == 0
        assert result.snps == []

    def test_purity_extraction(self, locus_set_two_pop, monkeypatch):
        """Purity values should be extracted."""
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(
                n, cs_indices=[[0, 1], [5, 6]], purity_values=[0.9, 0.7]
            ),
        )
        result = run_multisusie(locus_set_two_pop)
        assert result.purity is not None
        assert len(result.purity) == 2

    def test_include_mask_filtering(self, locus_set_two_pop, monkeypatch):
        """CS with include_mask=False should be excluded."""
        n = 20

        def mock_fn(**kwargs):
            return SimpleNamespace(
                pip=np.zeros(n),
                sets=(
                    [[0, 1], [5, 6]],
                    np.array([0.9, 0.3]),
                    [0.95, 0.95],
                    [True, False],  # Second CS filtered out
                ),
            )

        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss", mock_fn
        )
        result = run_multisusie(locus_set_two_pop)
        assert result.n_cs == 1

    def test_cs_maps_back_to_snpids(self, locus_set_two_pop, monkeypatch):
        """CS indices should map back to actual SNPIDs."""
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(n, cs_indices=[[0, 1]], purity_values=[0.9]),
        )
        result = run_multisusie(locus_set_two_pop)
        if result.n_cs > 0:
            for snp_list in result.snps:
                for snp in snp_list:
                    assert isinstance(snp, str)
                    assert "-" in snp  # SNPID format is "1-BP-A-G"


class TestRunMultisusieThreePop:
    """Tests for 3-population analysis."""

    def test_three_pop(self, locus_set_three_pop, monkeypatch):
        """Should work with 3 populations."""
        captured_kwargs = {}

        def capturing_mock(**kwargs):
            captured_kwargs.update(kwargs)
            n = len(kwargs["z_list"][0])
            pip = np.zeros(n)
            pip[0] = 0.7
            return SimpleNamespace(
                pip=pip,
                sets=([[0, 1]], np.array([0.8]), [0.95], [True]),
            )

        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss", capturing_mock
        )
        result = run_multisusie(locus_set_three_pop)
        assert len(captured_kwargs["z_list"]) == 3
        assert len(captured_kwargs["R_list"]) == 3
        assert captured_kwargs["rho"].shape == (3, 3)
        assert result.n_cs == 1
