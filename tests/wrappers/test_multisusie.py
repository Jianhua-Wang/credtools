"""Tests for the MultiSuSiE multi-ancestry fine-mapping wrapper."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet
from credtools.wrappers.multisusie import run_multisusie

from .conftest import _make_locus, _make_locus_set


def _make_mock_multisusie_rss(
    n_variants,
    cs_indices=None,
    purity_values=None,
    converged=True,
    niter=10,
):
    """Factory to create mock multisusie_rss functions.

    Parameters
    ----------
    n_variants : int
        Number of variants in the union set.
    cs_indices : list of list of int or None
        Credible set variant indices. None for no CS.
    purity_values : list of float or None
        Purity values for each CS.
    converged : bool
        Convergence status to expose on the returned object.
    niter : int
        Number of iterations to expose on the returned object.
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
        purity_arr = (
            np.array(purity_values) if purity_values else np.full(len(cs_list), np.nan)
        )
        # include_mask: True for valid CS, False otherwise
        include_mask = [True] * len(cs_list)

        return SimpleNamespace(
            pip=pip,
            sets=(cs_list, purity_arr, [0.95] * len(cs_list), include_mask),
            converged=converged,
            niter=niter,
        )

    return mock_fn


class TestRunMultisusieBasic:
    """Basic MultiSuSiE functionality tests."""

    def test_basic_call_two_pop(self, locus_set_two_pop, monkeypatch):
        """Verify MultiSuSiE returns a CredibleSet for 2-population LocusSet."""
        n = 20  # expected union size (same SNPs across pops)
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(n, cs_indices=[[0, 1, 2]], purity_values=[0.85]),
        )
        result = run_multisusie(locus_set_two_pop)
        assert isinstance(result, CredibleSet)
        assert result.tool == Method.MULTISUSIE

    def test_output_structure(self, locus_set_two_pop, monkeypatch):
        """Verify MultiSuSiE result has correct structure."""
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
        result = run_multisusie(locus_set_two_pop, max_causal=5, coverage=0.99, rho=0.8)
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

        monkeypatch.setattr("credtools.wrappers.multisusie.multisusie_rss", mock_fn)
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


class TestRunMultisusieConvergence:
    """Tests for convergence reporting and empty-on-nonconvergence."""

    def test_converged_propagated_when_true(self, locus_set_two_pop, monkeypatch):
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(
                n,
                cs_indices=[[0, 1]],
                purity_values=[0.9],
                converged=True,
                niter=20,
            ),
        )
        result = run_multisusie(locus_set_two_pop)
        assert result.converged is True
        assert result.n_iter == 20

    def test_converged_propagated_when_false(self, locus_set_two_pop, monkeypatch):
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(
                n,
                cs_indices=[[0, 1]],
                purity_values=[0.5],
                converged=False,
                niter=100,
            ),
        )
        result = run_multisusie(locus_set_two_pop)
        assert result.converged is False
        assert result.n_iter == 100
        assert result.n_cs == 1  # CS preserved by default

    def test_empty_on_nonconvergence_zeros_cs(self, locus_set_two_pop, monkeypatch):
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(
                n,
                cs_indices=[[0, 1]],
                purity_values=[0.5],
                converged=False,
                niter=100,
            ),
        )
        result = run_multisusie(locus_set_two_pop, empty_on_nonconvergence=True)
        assert result.converged is False
        assert result.n_cs == 0
        assert result.snps == []
        assert (result.pips == 0).all()

    def test_empty_on_nonconvergence_keeps_converged_cs(
        self, locus_set_two_pop, monkeypatch
    ):
        n = 20
        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss",
            _make_mock_multisusie_rss(
                n,
                cs_indices=[[0, 1]],
                purity_values=[0.9],
                converged=True,
                niter=10,
            ),
        )
        result = run_multisusie(locus_set_two_pop, empty_on_nonconvergence=True)
        assert result.converged is True
        assert result.n_cs == 1


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


class TestMissingVariantInputs:
    """Variants missing from a population must not enter as observed zeros."""

    @staticmethod
    def _capture(monkeypatch):
        captured = {}

        def capturing_mock(**kwargs):
            captured.update(kwargs)
            n = len(kwargs["z_list"][0])
            return SimpleNamespace(pip=np.zeros(n), sets=([], np.array([]), [], []))

        monkeypatch.setattr(
            "credtools.wrappers.multisusie.multisusie_rss", capturing_mock
        )
        return captured

    @staticmethod
    def _subset(locus, keep):
        """Locus restricted to the variants at positions ``keep`` (sumstats and LD)."""
        keep = np.asarray(keep)
        ld_map = locus.ld.map.iloc[keep].reset_index(drop=True)
        return Locus(
            locus.popu,
            locus.cohort,
            locus.sample_size,
            locus.sumstats.iloc[keep].reset_index(drop=True),
            locus._locus_start,
            locus._locus_end,
            ld=LDMatrix(ld_map, locus.ld.r[np.ix_(keep, keep)]),
        )

    def test_partial_overlap_masks_missing_variants(self, monkeypatch):
        """A variant absent from one population has NaN z and an empty LD row."""
        locus_a = self._subset(_make_locus(popu="EUR", seed=1, add_af2=True), range(15))
        locus_b = self._subset(
            _make_locus(popu="EAS", seed=2, add_af2=True), range(5, 20)
        )
        captured = self._capture(monkeypatch)
        run_multisusie(LocusSet([locus_a, locus_b]))

        snps_a = set(locus_a.sumstats[ColName.SNPID])
        snps_b = set(locus_b.sumstats[ColName.SNPID])
        union = sorted(snps_a | snps_b, key=lambda s: int(s.split("-")[1]))
        assert len(union) == 20
        for k, own in enumerate([snps_a, snps_b]):
            z, R = captured["z_list"][k], captured["R_list"][k]
            observed = np.array([s in own for s in union])
            assert observed.sum() == 15
            assert np.all(np.isfinite(z[observed]))
            assert np.all(np.isnan(z[~observed]))
            assert np.allclose(np.diag(R)[observed], 1.0)
            assert np.all(R[~observed, :] == 0)
            assert np.all(R[:, ~observed] == 0)

    def test_sumstats_variants_missing_from_ld_are_masked(self, monkeypatch):
        """Variants in the sumstats but not in the LD panel are dropped, not kept with z."""
        full = _make_locus(popu="EUR", seed=1, add_af2=True)
        n_keep = 17
        ld_map = full.ld.map.iloc[:n_keep].reset_index(drop=True)
        short = Locus(
            full.popu,
            full.cohort,
            full.sample_size,
            full.sumstats,
            full._locus_start,
            full._locus_end,
            ld=LDMatrix(ld_map, full.ld.r[:n_keep, :n_keep]),
        )
        other = _make_locus(popu="EAS", seed=2, add_af2=True)
        captured = self._capture(monkeypatch)
        run_multisusie(LocusSet([short, other]))

        dropped = set(full.sumstats[ColName.SNPID]) - set(ld_map[ColName.SNPID])
        union = sorted(
            set(full.sumstats[ColName.SNPID]), key=lambda s: int(s.split("-")[1])
        )
        miss = np.array([s in dropped for s in union])
        assert miss.sum() == len(full.sumstats) - n_keep
        z, R = captured["z_list"][0], captured["R_list"][0]
        assert np.all(np.isnan(z[miss]))
        assert np.all(np.diag(R)[miss] == 0)
        assert np.allclose(np.diag(R)[~miss], 1.0)
