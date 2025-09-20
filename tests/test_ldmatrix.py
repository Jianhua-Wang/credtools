import numpy as np

from credtools.ldmatrix import load_ld_matrix


def test_load_ld_matrix_replaces_nan_in_text(tmp_path):
    contents = (
        "\n".join(
            [
                "1",
                "0.1\t1",
                "nan\t0.2\t1",
            ]
        )
        + "\n"
    )
    ld_path = tmp_path / "ld.txt"
    ld_path.write_text(contents)

    matrix = load_ld_matrix(str(ld_path))

    assert matrix.dtype == np.float32
    assert not np.isnan(matrix).any()
    assert matrix[2, 0] == 0.0


def test_load_ld_matrix_replaces_nan_in_npz(tmp_path):
    data = np.array(
        [
            [1.0, np.nan],
            [np.nan, 1.0],
        ],
        dtype=np.float32,
    )
    ld_path = tmp_path / "ld.npz"
    np.savez(ld_path, ld=data)

    matrix = load_ld_matrix(str(ld_path))

    assert matrix.dtype == np.float32
    assert not np.isnan(matrix).any()
    assert matrix[0, 1] == 0.0
