import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from credtools.cli import CombineCred, CombinePIP, Tool, run_fine_map


class DummyPool:
    """Minimal pool stub used to validate parallel execution paths."""

    def __init__(self, processes: int) -> None:
        self.processes = processes

    def __enter__(self):
        """Return the pool instance to mimic context manager behaviour."""
        return self

    def __exit__(self, exc_type, exc, tb):
        """Ignore context manager exit arguments."""
        return False

    def imap_unordered(self, func, iterable):
        """Yield results immediately in the call order for determinism."""
        for item in iterable:
            yield func(item)


def make_fake_process_task(output_root: Path):
    def _fake_process_task(task):  # pragma: no cover - exercised via run_fine_map
        _fake_process_task.calls.append(task["locus_id"])
        locus_dir = output_root / task["locus_id"]
        locus_dir.mkdir(parents=True, exist_ok=True)
        pip_path = locus_dir / "pips.txt.gz"
        with gzip.open(pip_path, "wt") as handle:
            handle.write("SNP\tCRED\nrs1\t1\n")

        return {
            "status": "success",
            "locus_id": task["locus_id"],
            "cs_records": [{"SNP": "rs1", "CRED": 1, "locus_id": task["locus_id"]}],
        }

    _fake_process_task.calls = []
    return _fake_process_task


@pytest.mark.parametrize("processes", [1, 3])
def test_run_fine_map_parallel(tmp_path, monkeypatch, processes):
    loci_df = pd.DataFrame(
        [
            {
                "prefix": "locus1",
                "popu": "EUR",
                "cohort": "C1",
                "sample_size": 1000,
                "chr": 1,
                "start": 10,
                "end": 100,
                "locus_id": "locus1",
            },
            {
                "prefix": "locus2",
                "popu": "EUR",
                "cohort": "C2",
                "sample_size": 1200,
                "chr": 2,
                "start": 20,
                "end": 200,
                "locus_id": "locus2",
            },
        ]
    )

    inputs_path = tmp_path / "loci.tsv"
    loci_df.to_csv(inputs_path, sep="	", index=False)

    output_dir = tmp_path / "out"

    fake_task_runner = make_fake_process_task(output_dir)
    monkeypatch.setattr("credtools.cli._process_fine_map_task", fake_task_runner)
    monkeypatch.setattr("credtools.cli.Pool", DummyPool)

    run_fine_map(
        inputs=str(inputs_path),
        outdir=str(output_dir),
        tool=Tool.finemap,
        max_causal=5,
        adaptive_max_causal=False,
        set_L_by_cojo=False,
        p_cutoff=5e-8,
        collinear_cutoff=0.9,
        window_size=10_000_000,
        maf_cutoff=0.01,
        diff_freq_cutoff=0.2,
        coverage=0.95,
        timeout_minutes=30.0,
        processes=processes,
        combine_cred=CombineCred.union,
        combine_pip=CombinePIP.max,
        jaccard_threshold=0.1,
        max_iter=100,
        estimate_residual_variance=False,
        min_abs_corr=0.5,
        convergence_tol=1e-3,
        calculate_lambda_s=False,
        log_file=None,
    )

    assert set(fake_task_runner.calls) == {"locus1", "locus2"}

    for locus in ("locus1", "locus2"):
        pip_path = output_dir / locus / "pips.txt.gz"
        assert pip_path.exists()
        with gzip.open(pip_path, "rt") as handle:
            lines = [line for line in handle.read().strip().splitlines() if line]
            assert len(lines) == 2

    summary_path = output_dir / "credible_sets_summary.txt.gz"
    assert summary_path.exists()
    summary_df = pd.read_csv(summary_path, sep="	")
    assert set(summary_df["locus_id"]) == {"locus1", "locus2"}

    params_path = output_dir / "parameters.json"
    with params_path.open() as handle:
        params = json.load(handle)
    assert params["parameters"]["processes"] == processes

    summary_log = output_dir / "run_summary.log"
    assert summary_log.exists()
