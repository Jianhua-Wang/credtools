"""Tests for credtools.preprocessing.prepare module."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName


# ---------------------------------------------------------------------------
# TestExtractLdMatrix
# ---------------------------------------------------------------------------
class TestExtractLdMatrix:
    """Tests for _extract_ld_matrix function."""

    def test_unsupported_format_raises(self):
        from credtools.preprocessing.prepare import _extract_ld_matrix

        with pytest.raises(ValueError, match="Unsupported LD format"):
            _extract_ld_matrix(
                genotype_prefix="test",
                chrom=1,
                start=1000,
                end=2000,
                output_prefix="out",
                ld_format="unknown_format",
                keep_intermediate=False,
            )

    def test_vcf_not_implemented(self):
        from credtools.preprocessing.prepare import _extract_ld_vcf

        with pytest.raises(NotImplementedError, match="VCF format"):
            _extract_ld_vcf(
                genotype_prefix="test",
                chrom=1,
                start=1000,
                end=2000,
                output_prefix="out",
                keep_intermediate=False,
            )


# ---------------------------------------------------------------------------
# TestHandleAlleleFlipping
# ---------------------------------------------------------------------------
class TestHandleAlleleFlipping:
    """Tests for _handle_allele_flipping function."""

    def test_no_flip_needed(self):
        """When alleles are already sorted, no flip should occur."""
        from credtools.preprocessing.prepare import _handle_allele_flipping

        ldmap = pd.DataFrame(
            {
                "A1": ["A", "C", "A"],
                "A2": ["G", "T", "T"],
            }
        )
        ld_matrix = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.2], [0.3, 0.2, 1.0]])
        result_map, result_ld = _handle_allele_flipping(ldmap.copy(), ld_matrix.copy())
        np.testing.assert_array_almost_equal(result_ld, ld_matrix)
        assert list(result_map["A1"]) == ["A", "C", "A"]

    def test_flip_occurs(self):
        """When A1 > A2 alphabetically, flip should happen."""
        from credtools.preprocessing.prepare import _handle_allele_flipping

        ldmap = pd.DataFrame(
            {
                "A1": ["G", "A"],  # First SNP needs flip (G > A)
                "A2": ["A", "T"],
            }
        )
        ld_matrix = np.array([[1.0, 0.5], [0.5, 1.0]])
        result_map, result_ld = _handle_allele_flipping(ldmap.copy(), ld_matrix.copy())
        # After flip, A1 should be sorted alphabetically
        assert result_map["A1"].iloc[0] == "A"
        assert result_map["A2"].iloc[0] == "G"
        # Diagonal should remain 1
        assert result_ld[0, 0] == 1.0
        assert result_ld[1, 1] == 1.0
        # Off-diagonal: row 0 flipped -> ld[0,:] *= -1 -> ld[0,1] = -0.5
        # Then col 0 flipped -> ld[:,0] *= -1 -> ld[1,0] = -0.5
        # But ld[0,1] stays -0.5 (only row was flipped for that element)
        # Actually: row flip makes ld[0,1]=-0.5, then col flip makes ld[0,1]=0.5
        # Wait - the code does: ld_matrix[swapped_indices] *= -1 (row),
        # then ld_matrix[:, swapped_indices] *= -1 (col).
        # For element [0,1]: row flip -> -0.5, col flip of col 0 doesn't affect [0,1]
        # So [0,1] = -0.5, then diagonal fill makes [0,0] = 1.0 again
        assert result_ld[0, 1] == -0.5

    def test_diagonal_preserved(self):
        """Diagonal should always be 1 after flipping."""
        from credtools.preprocessing.prepare import _handle_allele_flipping

        n = 4
        ldmap = pd.DataFrame(
            {
                "A1": ["T", "G", "C", "A"],
                "A2": ["A", "A", "A", "G"],
            }
        )
        rng = np.random.default_rng(42)
        ld = rng.uniform(-1, 1, (n, n))
        ld = (ld + ld.T) / 2
        np.fill_diagonal(ld, 1.0)
        _, result_ld = _handle_allele_flipping(ldmap.copy(), ld.copy())
        np.testing.assert_array_almost_equal(np.diag(result_ld), np.ones(n))


# ---------------------------------------------------------------------------
# TestIntersectSumstatsLd
# ---------------------------------------------------------------------------
class TestIntersectSumstatsLd:
    """Tests for _intersect_sumstats_ld function."""

    def _make_sumstats_and_ld(self, common_snps=3, extra_sumstats=2, extra_ld=1):
        """Create matching sumstats and LD data with some overlap."""
        bps = list(
            range(1000, 1000 + (common_snps + extra_sumstats + extra_ld) * 100, 100)
        )

        common_ids = [f"1-{bps[i]}-A-G" for i in range(common_snps)]
        extra_sum_ids = [f"1-{bps[common_snps + i]}-A-G" for i in range(extra_sumstats)]
        extra_ld_ids = [
            f"1-{bps[common_snps + extra_sumstats + i]}-C-T" for i in range(extra_ld)
        ]

        all_sum_ids = common_ids + extra_sum_ids
        all_ld_ids = common_ids + extra_ld_ids
        n_sum = len(all_sum_ids)
        n_ld = len(all_ld_ids)

        sumstats = pd.DataFrame(
            {
                ColName.SNPID: all_sum_ids,
                ColName.CHR: [1] * n_sum,
                ColName.BP: bps[:n_sum],
                ColName.EA: ["A"] * n_sum,
                ColName.NEA: ["G"] * n_sum,
                ColName.EAF: [0.3] * n_sum,
                ColName.BETA: [0.1] * n_sum,
                ColName.SE: [0.01] * n_sum,
                ColName.P: [1e-8] * n_sum,
            }
        )

        ld_matrix = np.eye(n_ld)
        ldmap = pd.DataFrame(
            {
                ColName.SNPID: all_ld_ids,
                ColName.CHR: [1] * n_ld,
                ColName.BP: bps[:common_snps] + bps[n_sum : n_sum + extra_ld],
                "A1": ["A"] * common_snps + ["C"] * extra_ld,
                "A2": ["G"] * common_snps + ["T"] * extra_ld,
            }
        )
        return sumstats, ld_matrix, ldmap

    def test_common_variants_intersected(self):
        from credtools.preprocessing.prepare import _intersect_sumstats_ld

        sumstats, ld_matrix, ldmap = self._make_sumstats_and_ld(
            common_snps=5, extra_sumstats=3, extra_ld=2
        )
        result_ss, result_ld, result_map = _intersect_sumstats_ld(
            sumstats, ld_matrix, ldmap
        )
        assert len(result_ss) == 5
        assert result_ld.shape == (5, 5)

    def test_no_common_variants(self):
        from credtools.preprocessing.prepare import _intersect_sumstats_ld

        sumstats = pd.DataFrame(
            {
                ColName.SNPID: ["1-1000-A-G"],
                ColName.CHR: [1],
                ColName.BP: [1000],
                ColName.EA: ["A"],
                ColName.NEA: ["G"],
                ColName.EAF: [0.3],
                ColName.BETA: [0.1],
                ColName.SE: [0.01],
                ColName.P: [1e-8],
            }
        )
        ldmap = pd.DataFrame(
            {
                ColName.SNPID: ["1-2000-C-T"],
                ColName.CHR: [1],
                ColName.BP: [2000],
                "A1": ["C"],
                "A2": ["T"],
            }
        )
        ld_matrix = np.eye(1)
        result_ss, result_ld, result_map = _intersect_sumstats_ld(
            sumstats, ld_matrix, ldmap
        )
        assert len(result_ss) == 0
        assert result_ld.shape == (0,)
