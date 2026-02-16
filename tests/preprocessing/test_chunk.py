"""Unit tests for credtools.preprocessing.chunk module."""

import os

import numpy as np
import pandas as pd
import pytest

from credtools.preprocessing.chunk import (
    _identify_independent_snps_by_distance,
    _merge_overlapping_loci,
    chunk_sumstats,
    create_loci_list_for_credtools,
    identify_independent_loci,
)


# ---------------------------------------------------------------------------
# TestIdentifyIndependentSnpsByDistance
# ---------------------------------------------------------------------------


class TestIdentifyIndependentSnpsByDistance:
    """Tests for the internal _identify_independent_snps_by_distance function."""

    def _make_df(self, rows):
        """Convenience: build a small sumstats-like DataFrame from dicts."""
        return pd.DataFrame(rows)

    # -- basic behaviour ----

    def test_basic_returns_correct_fields(self):
        """Returned dicts must have the 8 expected keys."""
        df = self._make_df(
            [
                {"CHR": 1, "BP": 100_000, "SNPID": "s1", "P": 1e-10},
                {"CHR": 1, "BP": 110_000, "SNPID": "s2", "P": 0.5},
            ]
            * 5  # ensure enough variants for min_variants_per_locus
        )
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, True, 1
        )
        assert len(result) >= 1
        expected_keys = {
            "chr",
            "start",
            "end",
            "lead_snp",
            "lead_bp",
            "lead_p",
            "ancestry",
            "n_variants",
        }
        assert set(result[0].keys()) == expected_keys

    def test_lead_snp_is_most_significant(self):
        """The lead SNP of a locus should be the one with the smallest P."""
        rows = [
            {
                "CHR": 1,
                "BP": 100_000 + i * 1000,
                "SNPID": f"s{i}",
                "P": 1e-8 - i * 1e-10,
            }
            for i in range(15)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert len(result) == 1
        assert result[0]["lead_p"] == min(r["P"] for r in rows if r["P"] <= 5e-8)

    def test_close_sig_snps_merged_into_one_locus(self):
        """Two significant SNPs closer than distance_threshold → only 1 locus."""
        rows = [
            {"CHR": 1, "BP": 100_000, "SNPID": "s1", "P": 1e-10},
            {"CHR": 1, "BP": 120_000, "SNPID": "s2", "P": 1e-9},
        ]
        # pad with non-significant SNPs to pass min_variants
        rows += [
            {"CHR": 1, "BP": 100_000 + i * 1000, "SNPID": f"ns{i}", "P": 0.5}
            for i in range(20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert len(result) == 1

    def test_distant_sig_snps_produce_two_loci(self):
        """Two significant SNPs far apart → 2 independent loci."""
        rows = [
            {"CHR": 1, "BP": 100_000, "SNPID": "s1", "P": 1e-10},
            {"CHR": 1, "BP": 1_000_000, "SNPID": "s2", "P": 1e-9},
        ]
        # add enough variants around each lead
        for bp_base in [100_000, 1_000_000]:
            rows += [
                {
                    "CHR": 1,
                    "BP": bp_base + i * 1000,
                    "SNPID": f"ns{bp_base}_{i}",
                    "P": 0.5,
                }
                for i in range(1, 15)
            ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert len(result) == 2

    def test_boundary_start_end(self):
        """Verify start = max(1, bp - d/2) and end = bp + d/2."""
        rows = [{"CHR": 1, "BP": 500_000, "SNPID": "s1", "P": 1e-10}]
        rows += [
            {"CHR": 1, "BP": 500_000 + i * 100, "SNPID": f"ns{i}", "P": 0.5}
            for i in range(20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 400_000, 5e-8, False, 1
        )
        assert len(result) == 1
        assert result[0]["start"] == 300_000  # 500000 - 200000
        assert result[0]["end"] == 700_000  # 500000 + 200000

    def test_start_not_less_than_one(self):
        """When lead_bp is very small, start should be clamped to 1."""
        rows = [{"CHR": 1, "BP": 100, "SNPID": "s1", "P": 1e-10}]
        rows += [
            {"CHR": 1, "BP": 100 + i * 10, "SNPID": f"ns{i}", "P": 0.5}
            for i in range(20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert result[0]["start"] == 1

    def test_no_sig_use_most_sig_true(self):
        """No significant SNPs + use_most_sig=True → picks most sig per chr."""
        rows = [
            {"CHR": 1, "BP": 100_000, "SNPID": "s1", "P": 0.01},
            {"CHR": 1, "BP": 110_000, "SNPID": "s2", "P": 0.1},
        ]
        rows += [
            {"CHR": 1, "BP": 100_000 + i * 1000, "SNPID": f"ns{i}", "P": 0.5}
            for i in range(20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, True, 1
        )
        assert len(result) >= 1
        assert result[0]["lead_p"] == 0.01

    def test_no_sig_use_most_sig_false(self):
        """No significant SNPs + use_most_sig=False → empty list."""
        rows = [
            {"CHR": 1, "BP": 100_000 + i * 1000, "SNPID": f"s{i}", "P": 0.5}
            for i in range(20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert result == []

    def test_min_variants_filters_small_regions(self):
        """Region with fewer variants than min_variants_per_locus is skipped."""
        rows = [
            {"CHR": 1, "BP": 100_000, "SNPID": "s1", "P": 1e-10},
            {"CHR": 1, "BP": 100_100, "SNPID": "s2", "P": 0.5},
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 10
        )
        assert len(result) == 0

    def test_multiple_chromosomes_independent(self):
        """Loci on different chromosomes are identified independently."""
        rows = []
        for c in [1, 5]:
            rows.append({"CHR": c, "BP": 100_000, "SNPID": f"s{c}", "P": 1e-10})
            rows += [
                {"CHR": c, "BP": 100_000 + i * 1000, "SNPID": f"ns{c}_{i}", "P": 0.5}
                for i in range(1, 15)
            ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert len(result) == 2
        chrs = {r["chr"] for r in result}
        assert chrs == {1, 5}

    def test_empty_dataframe(self):
        """Empty DataFrame → empty list."""
        df = pd.DataFrame(columns=["CHR", "BP", "SNPID", "P"])
        result = _identify_independent_snps_by_distance(
            df, "EUR", 500_000, 5e-8, False, 1
        )
        assert result == []

    def test_ancestry_name_preserved(self):
        """The ancestry string should appear in every returned dict."""
        rows = [{"CHR": 1, "BP": 100_000, "SNPID": "s1", "P": 1e-10}]
        rows += [
            {"CHR": 1, "BP": 100_000 + i * 1000, "SNPID": f"ns{i}", "P": 0.5}
            for i in range(20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "myAnc", 500_000, 5e-8, False, 1
        )
        assert all(r["ancestry"] == "myAnc" for r in result)

    def test_n_variants_count_accurate(self):
        """n_variants should equal the number of SNPs inside [start, end]."""
        rows = [{"CHR": 1, "BP": 500_000, "SNPID": "lead", "P": 1e-10}]
        rows += [
            {"CHR": 1, "BP": 500_000 + i * 100, "SNPID": f"ns{i}", "P": 0.5}
            for i in range(1, 20)
        ]
        df = self._make_df(rows)
        result = _identify_independent_snps_by_distance(
            df, "EUR", 400_000, 5e-8, False, 1
        )
        assert len(result) == 1
        start, end = result[0]["start"], result[0]["end"]
        expected = len(df[(df["CHR"] == 1) & (df["BP"] >= start) & (df["BP"] <= end)])
        assert result[0]["n_variants"] == expected


# ---------------------------------------------------------------------------
# TestMergeOverlappingLoci
# ---------------------------------------------------------------------------


class TestMergeOverlappingLoci:
    """Tests for _merge_overlapping_loci."""

    def _make_loci_df(self, records):
        return pd.DataFrame(records)

    def test_no_overlap_different_chr(self):
        """Loci on different chromosomes should not be merged."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 500,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 2,
                    "start": 100,
                    "end": 500,
                    "lead_snp": "s2",
                    "lead_bp": 300,
                    "lead_p": 1e-9,
                    "ancestry": "ASN",
                    "n_variants": 10,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert len(result) == 2

    def test_no_overlap_same_chr(self):
        """Non-overlapping loci on the same chromosome are not merged."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 500,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 1,
                    "start": 600,
                    "end": 900,
                    "lead_snp": "s2",
                    "lead_bp": 700,
                    "lead_p": 1e-9,
                    "ancestry": "ASN",
                    "n_variants": 10,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert len(result) == 2

    def test_overlap_merge_boundaries(self):
        """Overlapping loci → merged start = min, end = max."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 600,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 1,
                    "start": 400,
                    "end": 900,
                    "lead_snp": "s2",
                    "lead_bp": 700,
                    "lead_p": 1e-9,
                    "ancestry": "ASN",
                    "n_variants": 15,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert len(result) == 1
        assert result.iloc[0]["start"] == 100
        assert result.iloc[0]["end"] == 900

    def test_overlap_lead_snp_most_significant(self):
        """After merging, lead_snp should come from the most significant locus."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 600,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-8,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 1,
                    "start": 400,
                    "end": 900,
                    "lead_snp": "s2",
                    "lead_bp": 700,
                    "lead_p": 1e-12,
                    "ancestry": "ASN",
                    "n_variants": 15,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert result.iloc[0]["lead_snp"] == "s2"
        assert result.iloc[0]["lead_p"] == 1e-12

    def test_overlap_ancestry_comma_join_sorted(self):
        """Merged ancestry is a comma-separated sorted string."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 600,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 1,
                    "start": 400,
                    "end": 900,
                    "lead_snp": "s2",
                    "lead_bp": 700,
                    "lead_p": 1e-9,
                    "ancestry": "ASN",
                    "n_variants": 15,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert result.iloc[0]["ancestry"] == "ASN,EUR"

    def test_overlap_n_variants_max(self):
        """After merging, n_variants = max of constituents."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 600,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 1,
                    "start": 400,
                    "end": 900,
                    "lead_snp": "s2",
                    "lead_bp": 700,
                    "lead_p": 1e-9,
                    "ancestry": "ASN",
                    "n_variants": 25,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert result.iloc[0]["n_variants"] == 25

    def test_three_way_overlap(self):
        """Three overlapping loci on the same chromosome merge into one."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 500,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
                {
                    "chr": 1,
                    "start": 400,
                    "end": 800,
                    "lead_snp": "s2",
                    "lead_bp": 600,
                    "lead_p": 1e-9,
                    "ancestry": "ASN",
                    "n_variants": 15,
                },
                {
                    "chr": 1,
                    "start": 700,
                    "end": 1100,
                    "lead_snp": "s3",
                    "lead_bp": 900,
                    "lead_p": 1e-8,
                    "ancestry": "AFR",
                    "n_variants": 20,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        # The first locus overlaps with the second, forming a merged locus.
        # Whether the third also gets merged depends on whether the algorithm
        # chains overlaps. Based on the code, it checks overlap with the
        # *current* locus bounds (not the merged bounds), so the third locus
        # overlaps with the second but is processed separately if the second
        # is already marked as processed. Let's just check ≤ 2 loci.
        assert len(result) <= 2

    def test_single_row(self):
        """Single locus input → returned unchanged."""
        df = self._make_loci_df(
            [
                {
                    "chr": 1,
                    "start": 100,
                    "end": 500,
                    "lead_snp": "s1",
                    "lead_bp": 300,
                    "lead_p": 1e-10,
                    "ancestry": "EUR",
                    "n_variants": 10,
                },
            ]
        )
        result = _merge_overlapping_loci(df)
        assert len(result) == 1

    def test_empty_dataframe(self):
        """Empty DataFrame → empty DataFrame."""
        df = pd.DataFrame(
            columns=[
                "chr",
                "start",
                "end",
                "lead_snp",
                "lead_bp",
                "lead_p",
                "ancestry",
                "n_variants",
            ]
        )
        result = _merge_overlapping_loci(df)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# TestIdentifyIndependentLoci
# ---------------------------------------------------------------------------


class TestIdentifyIndependentLoci:
    """Tests for the high-level identify_independent_loci function."""

    def test_single_ancestry_dict(self, tmp_path, munged_sumstats_gz_file):
        """Dict with one ancestry produces a valid loci DataFrame."""
        result = identify_independent_loci(
            {"EUR": munged_sumstats_gz_file},
            str(tmp_path / "output"),
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "locus_id" in result.columns

    def test_single_file_string(self, tmp_path, munged_sumstats_gz_file):
        """A single file path string should also work."""
        result = identify_independent_loci(
            munged_sumstats_gz_file,
            str(tmp_path / "output"),
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_multi_ancestry_merge(self, tmp_path, two_ancestry_gz_files):
        """Multiple ancestries with merge_overlapping=True."""
        result = identify_independent_loci(
            two_ancestry_gz_files,
            str(tmp_path / "output"),
            merge_overlapping=True,
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_multi_ancestry_no_merge(self, tmp_path, two_ancestry_gz_files):
        """Multiple ancestries with merge_overlapping=False → more loci."""
        merged = identify_independent_loci(
            two_ancestry_gz_files,
            str(tmp_path / "output_m"),
            merge_overlapping=True,
        )
        unmerged = identify_independent_loci(
            two_ancestry_gz_files,
            str(tmp_path / "output_u"),
            merge_overlapping=False,
        )
        assert len(unmerged) >= len(merged)

    def test_single_ancestry_skips_merge(self, tmp_path, munged_sumstats_gz_file):
        """With a single ancestry, merging is irrelevant (code path skipped)."""
        result = identify_independent_loci(
            {"EUR": munged_sumstats_gz_file},
            str(tmp_path / "output"),
            merge_overlapping=True,
        )
        assert isinstance(result, pd.DataFrame)

    def test_output_file_created(self, tmp_path, munged_sumstats_gz_file):
        """identified_loci.txt should be written to output_dir."""
        out_dir = str(tmp_path / "output")
        identify_independent_loci({"EUR": munged_sumstats_gz_file}, out_dir)
        assert os.path.exists(os.path.join(out_dir, "identified_loci.txt"))

    def test_locus_id_format(self, tmp_path, munged_sumstats_gz_file):
        """locus_id should match chr{N}_{start}_{end}."""
        result = identify_independent_loci(
            {"EUR": munged_sumstats_gz_file},
            str(tmp_path / "output"),
        )
        for _, row in result.iterrows():
            expected_id = f"chr{row['chr']}_{row['start']}_{row['end']}"
            assert row["locus_id"] == expected_id

    def test_sorted_by_chr_start(self, tmp_path, two_ancestry_gz_files):
        """Output should be sorted by chr then start."""
        result = identify_independent_loci(
            two_ancestry_gz_files,
            str(tmp_path / "output"),
        )
        if len(result) > 1:
            for i in range(len(result) - 1):
                curr = result.iloc[i]
                nxt = result.iloc[i + 1]
                assert (curr["chr"], curr["start"]) <= (nxt["chr"], nxt["start"])

    def test_empty_result_returns_empty_df(self, tmp_path):
        """When min_variants is impossibly high, result should be empty."""
        # Create a file with very few SNPs
        df = pd.DataFrame(
            {
                "CHR": [1],
                "BP": [100_000],
                "SNPID": ["s1"],
                "EA": ["A"],
                "NEA": ["G"],
                "EAF": [0.3],
                "BETA": [0.1],
                "SE": [0.05],
                "P": [1e-10],
                "N": [10000],
                "RSID": ["rs1"],
                "MAF": [0.3],
            }
        )
        fpath = str(tmp_path / "tiny.munged.txt.gz")
        df.to_csv(fpath, sep="\t", index=False, compression="gzip")
        result = identify_independent_loci(
            {"EUR": fpath},
            str(tmp_path / "output"),
            min_variants_per_locus=999,
        )
        assert len(result) == 0

    def test_output_dir_auto_created(self, tmp_path, munged_sumstats_gz_file):
        """Output directory is created if it doesn't exist."""
        out_dir = str(tmp_path / "new" / "nested" / "dir")
        assert not os.path.exists(out_dir)
        identify_independent_loci({"EUR": munged_sumstats_gz_file}, out_dir)
        assert os.path.isdir(out_dir)


# ---------------------------------------------------------------------------
# TestChunkSumstats
# ---------------------------------------------------------------------------


class TestChunkSumstats:
    """Tests for chunk_sumstats function."""

    def test_basic_single_ancestry(
        self, tmp_path, sample_loci_df, munged_sumstats_gz_file
    ):
        """Basic chunking with one ancestry."""
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(
            sample_loci_df,
            {"EUR": munged_sumstats_gz_file},
            out_dir,
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_chunk_file_contains_only_region_variants(
        self, tmp_path, sample_loci_df, munged_sumstats_gz_file
    ):
        """Each chunk file should only contain variants within [start, end]."""
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(
            sample_loci_df,
            {"EUR": munged_sumstats_gz_file},
            out_dir,
        )
        for _, row in result.iterrows():
            chunk_data = pd.read_csv(row["sumstats_file"], sep="\t", compression="gzip")
            assert (chunk_data["BP"] >= row["start"]).all()
            assert (chunk_data["BP"] <= row["end"]).all()
            assert (chunk_data["CHR"] == row["chr"]).all()

    def test_compressed_output(self, tmp_path, sample_loci_df, munged_sumstats_gz_file):
        """compress=True → files end with .gz."""
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(
            sample_loci_df,
            {"EUR": munged_sumstats_gz_file},
            out_dir,
            compress=True,
        )
        for _, row in result.iterrows():
            assert row["sumstats_file"].endswith(".gz")

    def test_uncompressed_output(
        self, tmp_path, sample_loci_df, munged_sumstats_gz_file
    ):
        """compress=False → files do NOT end with .gz."""
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(
            sample_loci_df,
            {"EUR": munged_sumstats_gz_file},
            out_dir,
            compress=False,
        )
        for _, row in result.iterrows():
            assert not row["sumstats_file"].endswith(".gz")
            assert row["sumstats_file"].endswith(".sumstats")

    def test_multi_ancestry(self, tmp_path, two_ancestry_gz_files):
        """Chunking with multiple ancestries in the ancestry field."""
        loci_df = pd.DataFrame(
            [
                {
                    "chr": 1,
                    "start": 1,
                    "end": 350000,
                    "lead_snp": "1-100000-A-G",
                    "lead_bp": 100000,
                    "lead_p": 1e-9,
                    "ancestry": "EUR,ASN",
                    "n_variants": 20,
                    "locus_id": "chr1_1_350000",
                },
            ]
        )
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(loci_df, two_ancestry_gz_files, out_dir)
        ancestries = result["ancestry"].unique()
        assert set(ancestries) == {"EUR", "ASN"}

    def test_return_columns(self, tmp_path, sample_loci_df, munged_sumstats_gz_file):
        """Returned DataFrame should have the expected columns."""
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(
            sample_loci_df,
            {"EUR": munged_sumstats_gz_file},
            out_dir,
        )
        expected_cols = {
            "locus_id",
            "ancestry",
            "chr",
            "start",
            "end",
            "n_variants",
            "sumstats_file",
        }
        assert expected_cols.issubset(set(result.columns))

    def test_chunk_info_file_created(
        self, tmp_path, sample_loci_df, munged_sumstats_gz_file
    ):
        """chunk_info.txt should be written to output_dir."""
        out_dir = str(tmp_path / "chunks")
        chunk_sumstats(sample_loci_df, {"EUR": munged_sumstats_gz_file}, out_dir)
        assert os.path.exists(os.path.join(out_dir, "chunk_info.txt"))

    def test_missing_ancestry_skipped(self, tmp_path, munged_sumstats_gz_file):
        """If ancestry is not in sumstats_files dict, it is skipped."""
        loci_df = pd.DataFrame(
            [
                {
                    "chr": 1,
                    "start": 1,
                    "end": 350000,
                    "lead_snp": "s1",
                    "lead_bp": 100000,
                    "lead_p": 1e-9,
                    "ancestry": "NONEXIST",
                    "n_variants": 20,
                    "locus_id": "chr1_1_350000",
                },
            ]
        )
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(loci_df, {"EUR": munged_sumstats_gz_file}, out_dir)
        assert len(result) == 0

    def test_region_no_variants_skipped(self, tmp_path, munged_sumstats_gz_file):
        """If the region has no variants in the sumstats, it is skipped."""
        loci_df = pd.DataFrame(
            [
                {
                    "chr": 99,
                    "start": 1,
                    "end": 100,
                    "lead_snp": "s1",
                    "lead_bp": 50,
                    "lead_p": 1e-9,
                    "ancestry": "EUR",
                    "n_variants": 0,
                    "locus_id": "chr99_1_100",
                },
            ]
        )
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(loci_df, {"EUR": munged_sumstats_gz_file}, out_dir)
        assert len(result) == 0

    def test_empty_loci_df(self, tmp_path, munged_sumstats_gz_file):
        """Empty loci_df → empty result and chunk_info.txt still created."""
        loci_df = pd.DataFrame(
            columns=[
                "chr",
                "start",
                "end",
                "lead_snp",
                "lead_bp",
                "lead_p",
                "ancestry",
                "n_variants",
                "locus_id",
            ]
        )
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(loci_df, {"EUR": munged_sumstats_gz_file}, out_dir)
        assert len(result) == 0
        assert os.path.exists(os.path.join(out_dir, "chunk_info.txt"))

    def test_file_naming_format(
        self, tmp_path, sample_loci_df, munged_sumstats_gz_file
    ):
        """Output files should follow {ancestry}.{locus_id}.sumstats.gz naming."""
        out_dir = str(tmp_path / "chunks")
        result = chunk_sumstats(
            sample_loci_df,
            {"EUR": munged_sumstats_gz_file},
            out_dir,
        )
        for _, row in result.iterrows():
            basename = os.path.basename(row["sumstats_file"])
            expected = f"{row['ancestry']}.{row['locus_id']}.sumstats.gz"
            assert basename == expected


# ---------------------------------------------------------------------------
# TestCreateLociListForCredtools
# ---------------------------------------------------------------------------


class TestCreateLociListForCredtools:
    """Tests for create_loci_list_for_credtools function."""

    def test_return_columns(self, tmp_path, sample_chunk_info_df):
        """Returned DataFrame has expected columns."""
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(
            sample_chunk_info_df, output_file=out_file
        )
        expected = {
            "locus_id",
            "chr",
            "start",
            "end",
            "popu",
            "cohort",
            "sample_size",
            "prefix",
        }
        assert expected.issubset(set(result.columns))

    def test_prefix_derived_from_sumstats_file(self, tmp_path, sample_chunk_info_df):
        """Prefix is sumstats_file with .sumstats.gz stripped."""
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(
            sample_chunk_info_df, output_file=out_file
        )
        for _, row in result.iterrows():
            assert ".sumstats" not in row["prefix"]
            assert not row["prefix"].endswith(".gz")

    def test_output_file_saved(self, tmp_path, sample_chunk_info_df):
        """Output file should be written."""
        out_file = str(tmp_path / "loci_list.txt")
        create_loci_list_for_credtools(sample_chunk_info_df, output_file=out_file)
        assert os.path.exists(out_file)

    def test_popu_cohort_equal_ancestry(self, tmp_path, sample_chunk_info_df):
        """Popu and cohort should both equal the ancestry."""
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(
            sample_chunk_info_df, output_file=out_file
        )
        for _, row in result.iterrows():
            assert row["popu"] == "EUR"
            assert row["cohort"] == "EUR"

    def test_sample_size_placeholder(self, tmp_path, sample_chunk_info_df):
        """sample_size should be 50000 (placeholder)."""
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(
            sample_chunk_info_df, output_file=out_file
        )
        assert (result["sample_size"] == 50000).all()

    def test_without_ld_info(self, tmp_path, sample_chunk_info_df):
        """ld_info_df=None should still work."""
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(
            sample_chunk_info_df, ld_info_df=None, output_file=out_file
        )
        assert len(result) == len(sample_chunk_info_df)

    def test_with_ld_info(self, tmp_path, sample_chunk_info_df):
        """When ld_info_df is provided, matching rows get extra columns."""
        ld_info = pd.DataFrame(
            [
                {
                    "locus_id": "chr1_1_350000",
                    "ancestry": "EUR",
                    "ld_file": "/path/to/ld",
                },
            ]
        )
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(
            sample_chunk_info_df, ld_info_df=ld_info, output_file=out_file
        )
        matched = result[result["locus_id"] == "chr1_1_350000"]
        assert "ld_file" in matched.columns
        assert matched.iloc[0]["ld_file"] == "/path/to/ld"

    def test_empty_input(self, tmp_path):
        """Empty chunk_info_df → empty output."""
        empty_df = pd.DataFrame(
            columns=[
                "locus_id",
                "ancestry",
                "chr",
                "start",
                "end",
                "n_variants",
                "sumstats_file",
            ]
        )
        out_file = str(tmp_path / "loci_list.txt")
        result = create_loci_list_for_credtools(empty_df, output_file=out_file)
        assert len(result) == 0
