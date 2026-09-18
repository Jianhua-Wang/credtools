"""Numerical and public-interface tests for matched geometric LD aggregation."""

import json
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from credtools.cli import app
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet
from credtools.meta import (
    ensure_meta_configuration,
    meta,
    meta_all,
    meta_by_population,
    meta_configuration,
    meta_lds,
    meta_loci,
    meta_locus,
    recover_completed_locus,
)


def cohort(
    name="a", pop="EUR", size=100, snps=(1, 2), ld_snps=None, se=None, beta=None, r=None
):
    """Build a small cohort whose GWAS and LD coverage can differ."""
    snps = list(snps)
    ld_snps = snps if ld_snps is None else list(ld_snps)
    ss = pd.DataFrame(
        {
            "SNPID": [f"1-{j}-A-G" for j in snps],
            "CHR": 1,
            "BP": snps,
            "EA": "A",
            "NEA": "G",
            "BETA": beta if beta is not None else np.full(len(snps), 0.2),
            "SE": se if se is not None else np.full(len(snps), 0.1),
            "EAF": 0.3,
            "P": 0.01,
        }
    )
    ldmap = pd.DataFrame(
        {
            "SNPID": [f"1-{j}-A-G" for j in ld_snps],
            "CHR": 1,
            "BP": ld_snps,
            "A1": "A",
            "A2": "G",
            "AF2": 0.7,
        }
    )
    if r is None:
        r = np.full((len(ld_snps), len(ld_snps)), 0.5)
        np.fill_diagonal(r, 1)
    return Locus(
        pop, name, size, ss, 1, 100, LDMatrix(ldmap, np.asarray(r, dtype=np.float32))
    )


@pytest.mark.parametrize("weighting", ["ess", "se"])
def test_geometric_missingness_not_pairwise(weighting):
    """A cohort carrying only SNP 1 dilutes correlation at that end only."""
    a = cohort(size=100, snps=(1, 2), se=[0.1, 0.1])
    b = cohort("b", size=300, snps=(1, 3), se=[0.1 / np.sqrt(3), 0.1])
    result = meta_all(LocusSet([a, b]), weighting)
    np.testing.assert_allclose(
        result.ld.r,
        [[1, 0.25, 0.5 * np.sqrt(0.75)], [0.25, 1, 0], [0.5 * np.sqrt(0.75), 0, 1]],
        atol=1e-7,
    )
    assert result.sumstats.SNPID.tolist() == ["1-1-A-G", "1-2-A-G", "1-3-A-G"]
    assert np.linalg.eigvalsh(result.ld.r).min() > 0


def test_default_ess_and_se_differ_only_in_ld():
    a = cohort(size=100, se=[0.1, 0.1], r=[[1, 0.2], [0.2, 1]])
    b = cohort("b", size=300, se=[0.2, 0.2], r=[[1, 0.8], [0.8, 1]])
    inputs = LocusSet([a, b])
    ess, se = meta_all(inputs), meta_all(inputs, "se")
    assert ess.ld.r[0, 1] == pytest.approx(0.65)
    assert se.ld.r[0, 1] == pytest.approx(0.32)
    pd.testing.assert_frame_equal(ess.sumstats, se.sumstats)
    np.testing.assert_array_equal(ess.ld.r, meta_lds(inputs, "ess").r)
    np.testing.assert_allclose(ess.sumstats.SE, 1 / np.sqrt(125))


def test_se_covariance_identity_with_heterogeneous_coverage():
    a = cohort(snps=(1, 2), se=[0.1, 0.3])
    b = cohort("b", snps=(2, 3), se=[0.2, 0.4], r=[[1, -0.2], [-0.2, 1]])
    result = meta_all(LocusSet([a, b]), "se")
    total = np.array([100, 1 / 0.3**2 + 25, 1 / 0.4**2])
    transform = np.zeros((3, 4))
    transform[0, 0] = 1
    transform[1, 1:3] = np.sqrt(np.array([1 / 0.3**2, 25]) / total[1])
    transform[2, 3] = 1
    covariance = np.zeros((4, 4))
    covariance[:2, :2], covariance[2:, 2:] = a.ld.r, b.ld.r
    np.testing.assert_allclose(
        result.ld.r, transform @ covariance @ transform.T, atol=1e-7
    )


@pytest.mark.parametrize("weighting", ["ess", "se"])
def test_matching_recomputes_statistics_without_original_seed(weighting):
    a = cohort(snps=(1, 2, 4), ld_snps=(1, 3), beta=[0.2, 8, 99])
    b = cohort("b", snps=(2, 3), ld_snps=(2, 3), beta=[0.4, 0.6])
    before = a.sumstats.copy()
    result = meta_all(LocusSet([a, b]), weighting)
    assert result.sumstats.SNPID.tolist() == ["1-1-A-G", "1-2-A-G", "1-3-A-G"]
    np.testing.assert_allclose(result.sumstats.BETA, [0.2, 0.4, 0.6])
    np.testing.assert_allclose(result.sumstats.SE, 0.1)
    audit = result.meta_audit.set_index("SNPID")
    assert audit.loc["1-2-A-G", "information_fraction"] == pytest.approx(0.5)
    assert not audit.loc["1-4-A-G", "retained"]
    pd.testing.assert_frame_equal(a.sumstats, before)
    pd.testing.assert_frame_equal(a.original_sumstats, before)


@pytest.mark.parametrize("weighting", ["ess", "se"])
def test_allele_flip_and_row_order(weighting):
    a = cohort(snps=(2, 1), ld_snps=(1, 2), beta=[-0.4, 0.2], r=[[1, -0.5], [-0.5, 1]])
    a.sumstats.loc[0, ["EA", "NEA", "EAF"]] = ["G", "A", 0.8]
    a.ld.map.loc[1, ["A1", "A2", "AF2"]] = ["G", "A", 0.2]
    original_r = a.ld.r.copy()
    result = meta_all(LocusSet([a]), weighting)
    np.testing.assert_allclose(result.sumstats.BETA, [0.2, 0.4])
    np.testing.assert_allclose(result.sumstats.EAF, [0.3, 0.2])
    np.testing.assert_allclose(result.ld.map.AF2, [0.7, 0.8])
    assert result.ld.r[0, 1] == pytest.approx(0.5)
    np.testing.assert_array_equal(original_r, a.ld.r)


@pytest.mark.parametrize("weighting", ["ess", "se"])
def test_invalid_contributions_and_empty_cohort(weighting):
    a = cohort(
        snps=(1, 2, 3, 4, 5), se=[0, -1, np.inf, 0.1, np.nan], beta=[1, 1, 1, np.nan, 1]
    )
    b = cohort("b", snps=(2,), ld_snps=(2,))
    result = meta_all(LocusSet([a, b]), weighting)
    assert result.n_snps == 1
    assert result.meta_cohort_audit.n_matched.tolist() == [0, 1]
    assert result.sample_size == 200
    with pytest.raises(ValueError, match="No SNP"):
        meta_all(LocusSet([a]), weighting)
    empty = cohort("empty", snps=(10,), ld_snps=())
    assert meta_all(LocusSet([empty, b]), weighting).n_snps == 1


@pytest.mark.parametrize("where", ["sumstats", "ld"])
def test_duplicate_and_mismatched_ids_rejected(where):
    a = cohort()
    frame = a.sumstats if where == "sumstats" else a.ld.map
    frame.loc[1, "SNPID"] = frame.loc[0, "SNPID"]
    with pytest.raises(ValueError, match="duplicate"):
        meta_all(LocusSet([a]))
    frame.loc[1, "SNPID"] = "1-3-A-G"
    with pytest.raises(ValueError, match="agree"):
        meta_all(LocusSet([a]))


@pytest.mark.parametrize(
    "r",
    [
        [[1, np.nan], [np.nan, 1]],
        [[1, 0.1], [0.3, 1]],
        [[0.8, 0.2], [0.2, 1]],
        [[1, 1.2], [1.2, 1]],
    ],
)
def test_invalid_ld_rejected_without_repair(r):
    with pytest.raises(ValueError, match="LD"):
        meta_all(LocusSet([cohort(r=r)]))


def test_indefinite_ld_not_implicitly_repaired():
    r = np.array([[1, 0.9, 0.9], [0.9, 1, -0.9], [0.9, -0.9, 1]])
    result = meta_all(LocusSet([cohort(snps=(1, 2, 3), r=r)]))
    np.testing.assert_allclose(result.ld.r, r, atol=1e-7)
    assert np.linalg.eigvalsh(result.ld.r).min() < -0.7


def test_population_normalizes_locally_and_singleton_matches():
    a = cohort(size=100)
    b = cohort("b", size=300, r=[[1, 0.9], [0.9, 1]])
    c = cohort("c", pop="AFR", snps=(1, 2, 3), ld_snps=(3, 2))
    result = meta_by_population(LocusSet([a, b, c]), "se")
    assert result["EUR"].ld.r[0, 1] == pytest.approx(0.7)
    assert result["AFR"].sumstats.SNPID.tolist() == ["1-2-A-G", "1-3-A-G"]
    assert result["AFR"].sample_size == 100
    np.testing.assert_allclose(result["AFR"].sumstats.SE, 0.1)
    assert len(meta(LocusSet([a, b, c]), "meta_by_population", "ess").loci) == 2


def test_frequencies_use_present_matched_ess_not_se_weights():
    a = cohort(size=100, se=[0.1, 0.1])
    b = cohort("b", size=300, se=[0.01, 0.01])
    a.sumstats["EAF"] = [0.2, np.nan]
    b.sumstats["EAF"] = [0.6, 0.8]
    np.testing.assert_allclose(
        meta_all(LocusSet([a, b]), "se").sumstats.EAF, [0.5, 0.8]
    )


def test_configuration_rejects_cross_mode_and_legacy(tmp_path):
    ensure_meta_configuration(str(tmp_path), "meta_all", "ess")
    ensure_meta_configuration(str(tmp_path), "meta_all", "ess")
    with pytest.raises(ValueError, match="separate"):
        ensure_meta_configuration(str(tmp_path), "meta_all", "se")
    with pytest.raises(ValueError, match="separate"):
        ensure_meta_configuration(str(tmp_path), "meta_by_population", "ess")
    legacy = tmp_path / "legacy"
    legacy.mkdir()
    (legacy / "loci_info.txt").touch()
    with pytest.raises(ValueError, match="no configuration"):
        ensure_meta_configuration(str(legacy), "meta_all", "ess")


@pytest.mark.parametrize("method", [meta_all, meta_lds, meta_by_population, meta])
def test_unknown_weighting_rejected(method):
    with pytest.raises(ValueError, match="ld_weighting"):
        method(LocusSet([cohort()]), ld_weighting="typo")


@pytest.mark.parametrize("weighting", ["ess", "se"])
def test_worker_roundtrip_and_mode_aware_recovery(tmp_path, weighting):
    inputs = LocusSet([cohort(), cohort("b", size=300, se=[0.2, 0.2])])
    with (
        patch("credtools.meta.load_locus_set", return_value=inputs),
        patch("credtools.meta.compute_heterogeneity", return_value={}),
        patch("credtools.meta.heterogeneity_summary", return_value=pd.DataFrame()),
    ):
        rows, _ = meta_locus(
            ("chr1_1_100", pd.DataFrame(), str(tmp_path), "meta_all", False, weighting)
        )
    prefix = rows[0][6]
    assert f"matched_{weighting}" in prefix
    stored = np.load(prefix + ".ld.npz")["ld"]
    assert stored.dtype == np.float16
    np.testing.assert_allclose(stored, meta_all(inputs, weighting).ld.r, atol=0.0005)
    ss = pd.read_csv(prefix + ".sumstats.gz", sep="\t")
    ldmap = pd.read_csv(prefix + ".ldmap.gz", sep="\t")
    assert ss.SNPID.equals(ldmap.SNPID)
    assert len(pd.read_csv(prefix + ".meta_audit.tsv.gz", sep="\t")) == 2
    info = pd.DataFrame(
        rows,
        columns=[
            "chr",
            "start",
            "end",
            "popu",
            "sample_size",
            "cohort",
            "prefix",
            "locus_id",
        ],
    )
    assert recover_completed_locus(
        "chr1_1_100", str(tmp_path), info, meta_configuration("meta_all", weighting)
    )
    other = "se" if weighting == "ess" else "ess"
    assert (
        recover_completed_locus(
            "chr1_1_100", str(tmp_path), info, meta_configuration("meta_all", other)
        )
        is None
    )


@pytest.mark.parametrize(
    "option,expected", [([], "ess"), (["--ld-weighting", "se"], "se")]
)
def test_meta_cli_forwards_weighting(tmp_path, option, expected):
    with patch("credtools.meta.meta_loci") as run:
        result = CliRunner().invoke(app, ["meta", "input.tsv", str(tmp_path)] + option)
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["ld_weighting"] == expected


@pytest.mark.parametrize("command", ["meta", "pipeline"])
def test_cli_rejects_unknown_weighting(command):
    result = CliRunner().invoke(
        app, [command, "input.tsv", "out", "--ld-weighting", "typo"]
    )
    assert result.exit_code == 2


def test_meta_loci_forwards_weighting_to_workers(tmp_path):
    data = pd.DataFrame(
        {
            "prefix": ["x"],
            "popu": ["EUR"],
            "cohort": ["a"],
            "sample_size": [100],
            "chr": [1],
            "start": [1],
            "end": [100],
            "locus_id": ["chr1_1_100"],
        }
    )
    path = tmp_path / "input.tsv"
    data.to_csv(path, sep="\t", index=False)
    with patch("credtools.meta.Pool") as pool:
        pool.return_value.__enter__.return_value.imap_unordered.return_value = []
        meta_loci(str(path), str(tmp_path / "out"), ld_weighting="se")
        arguments = (
            pool.return_value.__enter__.return_value.imap_unordered.call_args.args[1]
        )
    assert arguments[0][-1] == "se"
    assert (
        json.loads((tmp_path / "out" / "meta_config.json").read_text())["ld_weighting"]
        == "se"
    )


@pytest.mark.parametrize(
    "option,expected", [([], "ess"), (["--ld-weighting", "se"], "se")]
)
def test_pipeline_cli_forwards_weighting(tmp_path, option, expected):
    data = pd.DataFrame(
        {
            "prefix": ["x"],
            "popu": ["EUR"],
            "cohort": ["a"],
            "sample_size": [100],
            "chr": [1],
            "start": [1],
            "end": [100],
            "locus_id": ["chr1_1_100"],
        }
    )
    path = tmp_path / "input.tsv"
    data.to_csv(path, sep="\t", index=False)
    with patch("credtools.credtools.pipeline") as run:
        result = CliRunner().invoke(
            app, ["pipeline", str(path), str(tmp_path / "out")] + option
        )
    assert result.exit_code == 0, result.output
    assert run.call_args.kwargs["ld_weighting"] == expected


@pytest.mark.parametrize("weighting", ["ess", "se"])
def test_pipeline_api_passes_weighting_only_to_meta(tmp_path, weighting):
    from credtools.credtools import pipeline

    inputs = LocusSet([cohort(), cohort("b", size=300, se=[0.2, 0.2])])
    with (
        patch("credtools.credtools.load_locus_set", return_value=inputs),
        patch("credtools.credtools.compute_heterogeneity", return_value={}),
        patch("credtools.credtools.heterogeneity_summary", return_value=pd.DataFrame()),
        patch("credtools.credtools.meta", wraps=meta) as merge,
        patch(
            "credtools.credtools.fine_map",
            side_effect=RuntimeError("stop after input validation"),
        ) as fit,
    ):
        pipeline(
            pd.DataFrame(), outdir=str(tmp_path), skip_qc=True, ld_weighting=weighting
        )
    assert merge.call_args.kwargs["ld_weighting"] == weighting
    assert fit.call_count == 1
    assert "ld_weighting" not in fit.call_args.kwargs
    assert (tmp_path / "meta_config.json").exists()
    assert f"LD Weighting: {weighting}" in (tmp_path / "run_summary.log").read_text()


@pytest.mark.parametrize("size", [0, -1, np.nan, np.inf])
def test_invalid_ess_rejected(size):
    with pytest.raises(ValueError, match="sample_size"):
        meta_all(LocusSet([cohort(size=size)]))


def test_no_meta_does_not_change_with_weighting():
    inputs = LocusSet([cohort()])
    a = meta(inputs, "no_meta", "ess")
    b = meta(inputs, "no_meta", "se")
    pd.testing.assert_frame_equal(a.loci[0].sumstats, b.loci[0].sumstats)
    np.testing.assert_array_equal(a.loci[0].ld.r, b.loci[0].ld.r)


@pytest.mark.parametrize("weighting,expected_ld", [("ess", .725), ("se", .56)])
@pytest.mark.parametrize("method", ["meta_all", "meta_by_population"])
def test_real_cli_load_merge_save_and_resume(tmp_path, weighting, expected_ld, method):
    """Exercise file loading, the real worker pool, audits, and compatible resume."""
    bps = tuple(range(1000, 21000, 1000))
    b_r = np.full((20, 20), .8, dtype=np.float32)
    np.fill_diagonal(b_r, 1)
    inputs = [cohort(snps=bps), cohort("b", size=300, snps=bps, se=[.2] * 20, r=b_r)]
    rows = []
    for locus in inputs:
        prefix = tmp_path / locus.cohort
        locus.sumstats.to_csv(str(prefix) + ".sumstats.gz", sep="\t", index=False)
        locus.ld.map.to_csv(str(prefix) + ".ldmap.gz", sep="\t", index=False)
        np.savez_compressed(str(prefix) + ".ld.npz", ld=locus.ld.r)
        rows.append({"prefix": str(prefix), "popu": locus.popu, "cohort": locus.cohort,
                     "sample_size": locus.sample_size, "chr": 1, "start": 1, "end": 30000,
                     "locus_id": "chr1_1_30000"})
    info = tmp_path / "inputs.tsv"
    pd.DataFrame(rows).to_csv(info, sep="\t", index=False)
    output = tmp_path / "output"
    command = ["meta", str(info), str(output), "--meta-method", method, "--ld-weighting", weighting]
    runner = CliRunner()
    result = runner.invoke(app, command)
    assert result.exit_code == 0, (result.output, result.exception)
    merged = pd.read_csv(output / "loci_info.txt", sep="\t")
    assert len(merged) == 1
    assert merged.columns.tolist() == ["chr", "start", "end", "popu", "sample_size", "cohort", "prefix", "locus_id"]
    assert merged["chr"].iloc[0] == 1
    assert merged.sample_size.iloc[0] == 400
    prefix = merged.prefix.iloc[0]
    matrix = np.load(prefix + ".ld.npz")["ld"]
    assert matrix.dtype == np.float16
    assert float(matrix[0, 1]) == pytest.approx(expected_ld, abs=.0005)
    summary = pd.read_csv(prefix + ".sumstats.gz", sep="\t")
    np.testing.assert_allclose(summary.SE, 1 / np.sqrt(125), rtol=1e-6)
    from pathlib import Path
    timestamp = Path(prefix + ".ld.npz").stat().st_mtime_ns
    resumed = runner.invoke(app, command + ["--skip"])
    assert resumed.exit_code == 0, (resumed.output, resumed.exception)
    assert Path(prefix + ".ld.npz").stat().st_mtime_ns == timestamp
    changed = command[:-1] + (["se"] if weighting == "ess" else ["ess"])
    rejected = runner.invoke(app, changed + ["--skip"])
    assert rejected.exit_code != 0
    assert isinstance(rejected.exception, ValueError)
    assert "separate output directory" in str(rejected.exception)
