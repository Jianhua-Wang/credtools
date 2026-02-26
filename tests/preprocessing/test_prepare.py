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


# ---------------------------------------------------------------------------
# Imports for new tests
# ---------------------------------------------------------------------------
import os
from unittest.mock import MagicMock, call, patch


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_chunk_info(**overrides):
    """Create a pd.Series mimicking a chunk_info row."""
    defaults = {
        "locus_id": "locus_1",
        "chr": 1,
        "start": 1000,
        "end": 50000,
        "cohort": "cohort_A",
        "sample_size": 5000,
        "prefix": "/tmp/fake_prefix",
        "popu": "EUR",
    }
    defaults.update(overrides)
    return pd.Series(defaults)


def _make_chunk_df(n=2, ancestry="EUR"):
    """Create a small DataFrame resembling chunk_info_df."""
    rows = []
    for i in range(n):
        rows.append(
            {
                "locus_id": f"locus_{i}",
                "chr": 1,
                "start": 1000 + i * 1000,
                "end": 2000 + i * 1000,
                "cohort": "cohort_A",
                "sample_size": 5000,
                "prefix": f"/tmp/prefix_{i}",
                "popu": ancestry,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# TestPrepareFinemapInputs
# ---------------------------------------------------------------------------
class TestPrepareFinemapInputs:
    """Tests for prepare_finemap_inputs function."""

    @patch("credtools.preprocessing.prepare._prepare_ancestry_files")
    def test_single_threaded_processing(self, mock_prepare, tmp_path):
        """Single-threaded path delegates to _prepare_ancestry_files sequentially."""
        from credtools.preprocessing.prepare import prepare_finemap_inputs

        mock_prepare.return_value = [
            {"locus_id": "locus_0", "popu": "EUR", "status": "created"}
        ]

        chunk_df = _make_chunk_df(n=2, ancestry="EUR")
        genotype_files = {"EUR": "/data/eur_geno"}

        result = prepare_finemap_inputs(
            chunk_info_df=chunk_df,
            genotype_files=genotype_files,
            output_dir=str(tmp_path / "out"),
            threads=1,
        )

        mock_prepare.assert_called_once()
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        # Verify output file was written
        assert os.path.exists(os.path.join(str(tmp_path / "out"), "prepared_files.txt"))

    @patch("credtools.preprocessing.prepare._prepare_ancestry_files")
    def test_multi_threaded_processing(self, mock_prepare, tmp_path):
        """Multi-threaded path uses Pool.starmap."""
        from credtools.preprocessing.prepare import prepare_finemap_inputs

        # Pool.starmap will call the function directly in the subprocess,
        # but we need to patch at the module level. Instead, we patch Pool.
        mock_prepare.return_value = [
            {"locus_id": "locus_0", "popu": "EUR", "status": "created"}
        ]

        chunk_df = _make_chunk_df(n=2, ancestry="EUR")
        genotype_files = {"EUR": "/data/eur_geno"}

        # Patch Pool to avoid real multiprocessing
        with patch("credtools.preprocessing.prepare.Pool") as mock_pool_cls:
            mock_pool = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_pool)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
            mock_pool.starmap.return_value = [
                [{"locus_id": "locus_0", "popu": "EUR", "status": "created"}]
            ]

            result = prepare_finemap_inputs(
                chunk_info_df=chunk_df,
                genotype_files=genotype_files,
                output_dir=str(tmp_path / "out"),
                threads=2,
            )

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    @patch("credtools.preprocessing.prepare._prepare_ancestry_files")
    def test_ancestry_not_in_genotype_files_skipped(self, mock_prepare, tmp_path):
        """Ancestry without a genotype file entry is warned and skipped."""
        from credtools.preprocessing.prepare import prepare_finemap_inputs

        chunk_df = _make_chunk_df(n=1, ancestry="AFR")
        genotype_files = {"EUR": "/data/eur_geno"}  # AFR is missing

        result = prepare_finemap_inputs(
            chunk_info_df=chunk_df,
            genotype_files=genotype_files,
            output_dir=str(tmp_path / "out"),
            threads=1,
        )

        mock_prepare.assert_not_called()
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# TestPrepareAncestryFiles
# ---------------------------------------------------------------------------
class TestPrepareAncestryFiles:
    """Tests for _prepare_ancestry_files function."""

    @patch("credtools.preprocessing.prepare._prepare_single_locus")
    def test_normal_path_returns_results(self, mock_single):
        """Successful _prepare_single_locus results are collected."""
        from credtools.preprocessing.prepare import _prepare_ancestry_files

        mock_single.return_value = {
            "locus_id": "locus_0",
            "popu": "EUR",
            "status": "created",
        }

        chunk_group = _make_chunk_df(n=2, ancestry="EUR")

        result = _prepare_ancestry_files(
            ancestry="EUR",
            genotype_prefix="/data/eur",
            chunk_group=chunk_group,
            output_dir="/out",
            ld_format="plink",
            keep_intermediate=False,
            kwargs={},
        )

        assert len(result) == 2
        assert mock_single.call_count == 2

    @patch("credtools.preprocessing.prepare._prepare_single_locus")
    def test_exception_logged_and_continues(self, mock_single):
        """When _prepare_single_locus raises, the error is logged and processing continues."""
        from credtools.preprocessing.prepare import _prepare_ancestry_files

        mock_single.side_effect = [
            RuntimeError("boom"),
            {"locus_id": "locus_1", "popu": "EUR", "status": "created"},
        ]

        chunk_group = _make_chunk_df(n=2, ancestry="EUR")

        result = _prepare_ancestry_files(
            ancestry="EUR",
            genotype_prefix="/data/eur",
            chunk_group=chunk_group,
            output_dir="/out",
            ld_format="plink",
            keep_intermediate=False,
            kwargs={},
        )

        assert len(result) == 1
        assert result[0]["locus_id"] == "locus_1"

    @patch("credtools.preprocessing.prepare._prepare_single_locus")
    def test_none_result_not_appended(self, mock_single):
        """When _prepare_single_locus returns None, it is not included in the results."""
        from credtools.preprocessing.prepare import _prepare_ancestry_files

        mock_single.return_value = None

        chunk_group = _make_chunk_df(n=3, ancestry="EUR")

        result = _prepare_ancestry_files(
            ancestry="EUR",
            genotype_prefix="/data/eur",
            chunk_group=chunk_group,
            output_dir="/out",
            ld_format="plink",
            keep_intermediate=False,
            kwargs={},
        )

        assert len(result) == 0
        assert mock_single.call_count == 3


# ---------------------------------------------------------------------------
# TestPrepareSingleLocus
# ---------------------------------------------------------------------------
class TestPrepareSingleLocus:
    """Tests for _prepare_single_locus function."""

    def test_output_files_exist_returns_existed(self, tmp_path):
        """When all expected output files already exist, return dict with status='existed'."""
        from credtools.preprocessing.prepare import _prepare_single_locus

        chunk_info = _make_chunk_info()
        ancestry = "EUR"
        output_dir = str(tmp_path)
        output_prefix = os.path.join(output_dir, f"{ancestry}.{chunk_info['locus_id']}")

        # Create the expected output files
        for suffix in [".sumstats.gz", ".ld.npz", ".ldmap.gz"]:
            with open(output_prefix + suffix, "w") as f:
                f.write("placeholder")

        result = _prepare_single_locus(
            chunk_info=chunk_info,
            ancestry=ancestry,
            genotype_prefix="/data/eur",
            output_dir=output_dir,
            ld_format="plink",
            keep_intermediate=False,
        )

        assert result is not None
        assert result["status"] == "existed"
        assert result["locus_id"] == "locus_1"
        assert result["popu"] == "EUR"
        assert result["chr"] == 1

    @patch("credtools.preprocessing.prepare.subprocess")
    @patch("credtools.preprocessing.prepare._handle_allele_flipping")
    @patch("credtools.preprocessing.prepare._extract_ld_matrix")
    @patch("credtools.preprocessing.prepare.make_SNPID_unique")
    @patch("credtools.preprocessing.prepare.munge")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    def test_normal_flow_returns_created(
        self,
        mock_read_csv,
        mock_munge,
        mock_unique,
        mock_extract_ld,
        mock_flip,
        mock_subprocess,
        tmp_path,
    ):
        """Normal processing flow reads sumstats, extracts LD, and saves files."""
        from credtools.preprocessing.prepare import _prepare_single_locus

        # Set up mock chain
        fake_sumstats = pd.DataFrame(
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
        mock_read_csv.return_value = fake_sumstats.copy()
        mock_munge.return_value = fake_sumstats.copy()
        mock_unique.return_value = fake_sumstats.copy()

        fake_ldmap = pd.DataFrame(
            {
                ColName.CHR: [1],
                ColName.BP: [1000],
                "A1": ["A"],
                "A2": ["G"],
                ColName.SNPID: ["1-1000-A-G"],
            }
        )
        fake_ld = np.array([[1.0]])
        mock_extract_ld.return_value = (fake_ldmap, fake_ld)
        mock_flip.return_value = (fake_ldmap, fake_ld)

        mock_subprocess.run.return_value = MagicMock(returncode=0)

        chunk_info = _make_chunk_info()
        output_dir = str(tmp_path)

        result = _prepare_single_locus(
            chunk_info=chunk_info,
            ancestry="EUR",
            genotype_prefix="/data/eur",
            output_dir=output_dir,
            ld_format="plink",
            keep_intermediate=False,
        )

        assert result is not None
        assert result["status"] == "created"
        assert result["n_variants"] == 1
        assert result["locus_id"] == "locus_1"
        mock_munge.assert_called_once()
        mock_extract_ld.assert_called_once()

    @patch("credtools.preprocessing.prepare._extract_ld_matrix")
    @patch("credtools.preprocessing.prepare.make_SNPID_unique")
    @patch("credtools.preprocessing.prepare.munge")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    def test_ld_extraction_returns_none(
        self, mock_read_csv, mock_munge, mock_unique, mock_extract_ld, tmp_path
    ):
        """When _extract_ld_matrix returns None, function returns None."""
        from credtools.preprocessing.prepare import _prepare_single_locus

        fake_sumstats = pd.DataFrame(
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
        mock_read_csv.return_value = fake_sumstats.copy()
        mock_munge.return_value = fake_sumstats.copy()
        mock_unique.return_value = fake_sumstats.copy()
        mock_extract_ld.return_value = None

        chunk_info = _make_chunk_info()

        result = _prepare_single_locus(
            chunk_info=chunk_info,
            ancestry="EUR",
            genotype_prefix="/data/eur",
            output_dir=str(tmp_path),
            ld_format="plink",
            keep_intermediate=False,
        )

        assert result is None

    @patch("credtools.preprocessing.prepare._extract_ld_matrix")
    @patch("credtools.preprocessing.prepare.make_SNPID_unique")
    @patch("credtools.preprocessing.prepare.munge")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    def test_exception_cleans_up_partial_files(
        self, mock_read_csv, mock_munge, mock_unique, mock_extract_ld, tmp_path
    ):
        """When an exception occurs, partial output files are removed."""
        from credtools.preprocessing.prepare import _prepare_single_locus

        mock_read_csv.side_effect = FileNotFoundError("sumstats not found")

        chunk_info = _make_chunk_info()
        output_dir = str(tmp_path)
        output_prefix = os.path.join(output_dir, f"EUR.{chunk_info['locus_id']}")

        # Create a partial file that should be cleaned up
        partial_file = f"{output_prefix}.sumstats.gz"
        with open(partial_file, "w") as f:
            f.write("partial data")

        result = _prepare_single_locus(
            chunk_info=chunk_info,
            ancestry="EUR",
            genotype_prefix="/data/eur",
            output_dir=output_dir,
            ld_format="plink",
            keep_intermediate=False,
        )

        assert result is None
        # The partial file should have been cleaned up
        assert not os.path.exists(partial_file)


# ---------------------------------------------------------------------------
# TestExtractLdPlink
# ---------------------------------------------------------------------------
class TestExtractLdPlink:
    """Tests for _extract_ld_plink function."""

    @patch("credtools.preprocessing.prepare.os.path.exists")
    @patch("credtools.preprocessing.prepare.os.remove")
    @patch("credtools.preprocessing.prepare.np.loadtxt")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_plink_extraction_succeeds(
        self, mock_run, mock_read_csv, mock_loadtxt, mock_remove, mock_exists
    ):
        """Successful plink extraction returns (ldmap, ld_matrix)."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        # All plink commands succeed
        mock_run.return_value = MagicMock(returncode=0, stderr="", stdout="")

        # bim file exists, ld file exists, frq file exists
        mock_exists.return_value = True

        # bim data
        bim_df = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "RSID": ["rs1", "rs2", "rs3"],
                "CM": [0, 0, 0],
                "BP": [1000, 2000, 3000],
                "A1": ["A", "C", "G"],
                "A2": ["G", "T", "A"],
            }
        )

        # freq data
        freq_df = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "SNP": ["rs1", "rs2", "rs3"],
                "A1": ["A", "C", "G"],
                "A2": ["G", "T", "A"],
                "MAF": [0.1, 0.2, 0.15],
                "NCHROBS": [10000, 10000, 10000],
            }
        )

        # read_csv is called twice: once for bim, once for frq
        mock_read_csv.side_effect = [bim_df, freq_df]

        # LD matrix
        ld_data = np.eye(3)
        mock_loadtxt.return_value = ld_data

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=False,
        )

        assert result is not None
        ldmap, ld_matrix = result
        assert len(ldmap) == 3
        assert ld_matrix.shape == (3, 3)
        assert "A1" in ldmap.columns
        assert "A2" in ldmap.columns
        assert "AF2" in ldmap.columns
        # Verify plink was called 3 times (extract, ld, freq)
        assert mock_run.call_count == 3

    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_plink_command_fails_returns_none(self, mock_run):
        """When plink extraction command fails (returncode != 0), returns None."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        mock_run.return_value = MagicMock(
            returncode=1, stderr="Error: file not found", stdout=""
        )

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=False,
        )

        assert result is None

    @patch("credtools.preprocessing.prepare.os.path.exists")
    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_bim_file_not_exists_returns_none(self, mock_run, mock_exists):
        """When bim file does not exist after extraction, returns None."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        mock_run.return_value = MagicMock(returncode=0, stderr="", stdout="")
        # bim file does not exist
        mock_exists.return_value = False

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=False,
        )

        assert result is None

    @patch("credtools.preprocessing.prepare.os.path.exists")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_fewer_than_two_variants_returns_none(
        self, mock_run, mock_read_csv, mock_exists
    ):
        """When fewer than 2 variants are found, returns None."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        mock_run.return_value = MagicMock(returncode=0, stderr="", stdout="")
        mock_exists.return_value = True

        # Only one variant in bim
        bim_df = pd.DataFrame(
            {
                "CHR": [1],
                "RSID": ["rs1"],
                "CM": [0],
                "BP": [1000],
                "A1": ["A"],
                "A2": ["G"],
            }
        )
        mock_read_csv.return_value = bim_df

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=False,
        )

        assert result is None

    @patch("credtools.preprocessing.prepare.os.path.exists")
    @patch("credtools.preprocessing.prepare.os.remove")
    @patch("credtools.preprocessing.prepare.np.loadtxt")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_ld_computation_fails_returns_none(
        self, mock_run, mock_read_csv, mock_loadtxt, mock_remove, mock_exists
    ):
        """When the second plink call (LD computation) fails, returns None."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        # First call succeeds (extract), second fails (LD)
        mock_run.side_effect = [
            MagicMock(returncode=0, stderr="", stdout=""),
            MagicMock(returncode=1, stderr="LD computation error", stdout=""),
        ]
        mock_exists.return_value = True

        bim_df = pd.DataFrame(
            {
                "CHR": [1, 1],
                "RSID": ["rs1", "rs2"],
                "CM": [0, 0],
                "BP": [1000, 2000],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        mock_read_csv.return_value = bim_df

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=False,
        )

        assert result is None

    @patch("credtools.preprocessing.prepare.os.path.exists")
    @patch("credtools.preprocessing.prepare.os.remove")
    @patch("credtools.preprocessing.prepare.np.loadtxt")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_freq_computation_fails_returns_none(
        self, mock_run, mock_read_csv, mock_loadtxt, mock_remove, mock_exists
    ):
        """When the third plink call (freq) fails, returns None."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        # First call succeeds (extract), second succeeds (LD), third fails (freq)
        mock_run.side_effect = [
            MagicMock(returncode=0, stderr="", stdout=""),
            MagicMock(returncode=0, stderr="", stdout=""),
            MagicMock(returncode=1, stderr="freq error", stdout=""),
        ]

        def exists_side_effect(path):
            return True

        mock_exists.side_effect = exists_side_effect

        bim_df = pd.DataFrame(
            {
                "CHR": [1, 1],
                "RSID": ["rs1", "rs2"],
                "CM": [0, 0],
                "BP": [1000, 2000],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        mock_read_csv.return_value = bim_df
        mock_loadtxt.return_value = np.eye(2)

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=False,
        )

        assert result is None

    @patch("credtools.preprocessing.prepare.os.path.exists")
    @patch("credtools.preprocessing.prepare.os.remove")
    @patch("credtools.preprocessing.prepare.np.loadtxt")
    @patch("credtools.preprocessing.prepare.pd.read_csv")
    @patch("credtools.preprocessing.prepare.subprocess.run")
    def test_keep_intermediate_preserves_files(
        self, mock_run, mock_read_csv, mock_loadtxt, mock_remove, mock_exists
    ):
        """When keep_intermediate=True, temp files are not removed."""
        from credtools.preprocessing.prepare import _extract_ld_plink

        mock_run.return_value = MagicMock(returncode=0, stderr="", stdout="")
        mock_exists.return_value = True

        bim_df = pd.DataFrame(
            {
                "CHR": [1, 1],
                "RSID": ["rs1", "rs2"],
                "CM": [0, 0],
                "BP": [1000, 2000],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        freq_df = pd.DataFrame(
            {
                "CHR": [1, 1],
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
                "MAF": [0.1, 0.2],
                "NCHROBS": [10000, 10000],
            }
        )
        mock_read_csv.side_effect = [bim_df, freq_df]
        mock_loadtxt.return_value = np.eye(2)

        result = _extract_ld_plink(
            genotype_prefix="/data/eur",
            chrom=1,
            start=1000,
            end=3000,
            output_prefix="/out/test",
            keep_intermediate=True,
        )

        assert result is not None
        # os.remove should not have been called for temp files
        mock_remove.assert_not_called()
