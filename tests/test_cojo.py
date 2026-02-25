"""Tests for credtools.cojo module."""

from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus


def _make_cojo_locus(n_snps=20, p_range=(1e-10, 1e-9), add_af2=False, seed=42):
    """Create a Locus for COJO testing."""
    rng = np.random.default_rng(seed)
    bps = [1000 + i * 100 for i in range(n_snps)]
    snpids = [f"1-{bp}-A-G" for bp in bps]

    sumstats = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: [1] * n_snps,
            ColName.BP: bps,
            ColName.RSID: snpids,
            ColName.EA: ["A"] * n_snps,
            ColName.NEA: ["G"] * n_snps,
            ColName.EAF: rng.uniform(0.1, 0.9, size=n_snps),
            ColName.BETA: rng.normal(0, 0.1, size=n_snps),
            ColName.SE: np.abs(rng.normal(0.02, 0.005, size=n_snps)),
            ColName.P: rng.uniform(*p_range, size=n_snps),
        }
    )

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
    ld = LDMatrix(ld_map, r)
    return Locus("EUR", "test", 10000, sumstats, 1000, 3000, ld=ld)


class TestConditionalSelection:
    """Tests for conditional_selection function."""

    def test_normal_path(self):
        """Test that conditional_selection returns a DataFrame with mock COJO."""
        locus = _make_cojo_locus(add_af2=True)
        fake_result = pd.DataFrame({"SNP": ["1-1000-A-G"], "b": [0.1], "p": [1e-9]})

        fake_cojo_instance = MagicMock()
        fake_cojo_instance.conditional_selection.return_value = fake_result

        with patch("credtools.cojo.COJO", return_value=fake_cojo_instance):
            from credtools.cojo import conditional_selection

            result = conditional_selection(locus)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        fake_cojo_instance.load_sumstats.assert_called_once()

    def test_p_cutoff_relaxed_when_no_snps_pass(self):
        """When no SNPs pass p_cutoff, it should be relaxed to 1e-5."""
        locus = _make_cojo_locus(p_range=(0.01, 0.1))
        fake_result = pd.DataFrame({"SNP": [], "b": [], "p": []})

        fake_cojo_instance = MagicMock()
        fake_cojo_instance.conditional_selection.return_value = fake_result

        with patch("credtools.cojo.COJO", return_value=fake_cojo_instance) as MockCOJO:
            from credtools.cojo import conditional_selection

            result = conditional_selection(locus, p_cutoff=5e-8)

        # COJO should have been created with relaxed p_cutoff=1e-5
        call_kwargs = MockCOJO.call_args[1]
        assert call_kwargs["p_cutoff"] == 1e-5

    def test_af2_missing_sets_ld_freq_none(self):
        """When AF2 is not in LD map, ld_freq should be set to None."""
        locus = _make_cojo_locus(add_af2=False)
        fake_result = pd.DataFrame({"SNP": ["1-1000-A-G"], "b": [0.1], "p": [1e-9]})

        fake_cojo_instance = MagicMock()
        fake_cojo_instance.conditional_selection.return_value = fake_result

        with patch("credtools.cojo.COJO", return_value=fake_cojo_instance):
            from credtools.cojo import conditional_selection

            result = conditional_selection(locus)

        # Check that ld_freq was passed as None
        call_kwargs = fake_cojo_instance.load_sumstats.call_args[1]
        assert call_kwargs["ld_freq"] is None

    def test_af2_present_creates_freq_df(self):
        """When AF2 is in LD map, ld_freq should be a DataFrame."""
        locus = _make_cojo_locus(add_af2=True)
        fake_result = pd.DataFrame({"SNP": ["1-1000-A-G"], "b": [0.1], "p": [1e-9]})

        fake_cojo_instance = MagicMock()
        fake_cojo_instance.conditional_selection.return_value = fake_result

        with patch("credtools.cojo.COJO", return_value=fake_cojo_instance):
            from credtools.cojo import conditional_selection

            result = conditional_selection(locus)

        call_kwargs = fake_cojo_instance.load_sumstats.call_args[1]
        ld_freq = call_kwargs["ld_freq"]
        assert isinstance(ld_freq, pd.DataFrame)
        assert list(ld_freq.columns) == ["SNP", "freq"]
