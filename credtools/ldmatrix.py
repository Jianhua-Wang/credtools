"""Functions for reading and converting lower triangle matrices."""

import gzip
import logging
from typing import Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName
from credtools.sumstats import make_SNPID_unique, munge_bp, munge_chr

logger = logging.getLogger("LDMatrix")


class LDMatrix:
    """
    Class to store the LD matrix and the corresponding Variant IDs.

    Parameters
    ----------
    map_df : pd.DataFrame
        DataFrame containing the Variant IDs.
    r : np.ndarray
        LD matrix.

    Attributes
    ----------
    map : pd.DataFrame
        DataFrame containing the Variant IDs.
    r : np.ndarray
        LD matrix.

    Raises
    ------
    ValueError
        If the number of rows in the map file does not match the number of rows in the LD matrix.
    """

    def __init__(self, map_df: pd.DataFrame, r: np.ndarray) -> None:
        """
        Initialize the LDMatrix object.

        Parameters
        ----------
        map_df : pd.DataFrame
            DataFrame containing the Variant IDs.
        r : np.ndarray
            LD matrix.

        Raises
        ------
        ValueError
            If the number of rows in the map file does not match the number of rows in the LD matrix.
        """
        self.map = map_df
        self.r = r
        self.__check_length()

    def __repr__(self) -> str:
        """
        Return a string representation of the LDMatrix object.

        Returns
        -------
        str
            String representation showing the shapes of map and r.
        """
        return f"LDMatrix(map={self.map.shape}, r={self.r.shape})"

    def __check_length(self) -> None:
        """
        Check if the number of rows in the map file matches the number of rows in the LD matrix.

        Raises
        ------
        ValueError
            If the number of rows in the map file does not match the number of rows in the LD matrix.
        """
        if len(self.map) != len(self.r):
            raise ValueError(
                "The number of rows in the map file does not match the number of rows in the LD matrix."
            )

    def copy(self) -> "LDMatrix":
        """
        Return a copy of the LDMatrix object.

        Returns
        -------
        LDMatrix
            A copy of the LDMatrix object.
        """
        return LDMatrix(self.map.copy(), self.r.copy())


def read_lower_triangle(file_path: str, delimiter: str = "\t") -> np.ndarray:
    r"""
    Read a lower triangle matrix from a file.

    Parameters
    ----------
    file_path : str
        Path to the input text file containing the lower triangle matrix.
    delimiter : str, optional
        Delimiter used in the input file, by default "\t".

    Returns
    -------
    np.ndarray
        Lower triangle matrix.

    Raises
    ------
    ValueError
        If the input file is empty or does not contain a valid lower triangle matrix.
    FileNotFoundError
        If the specified file does not exist.

    Notes
    -----
    This function reads a lower triangular matrix where each row contains
    elements from the diagonal down to that row position.
    """
    try:
        if file_path.endswith(".gz"):
            with gzip.open(file_path, "rt") as file:
                rows = [
                    list(map(float, line.strip().split(delimiter)))
                    for line in file
                    if line.strip()
                ]
        else:
            with open(file_path, "r") as file:
                rows = [
                    list(map(float, line.strip().split(delimiter)))
                    for line in file
                    if line.strip()
                ]
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    if not rows:
        raise ValueError("The input file is empty.")

    n = len(rows)
    lower_triangle = np.zeros((n, n))

    for i, row in enumerate(rows):
        if len(row) != i + 1:
            raise ValueError(
                f"Invalid number of elements in row {i + 1}. Expected {i + 1}, got {len(row)}."
            )
        lower_triangle[i, : len(row)] = row

    return lower_triangle


def load_ld_matrix(file_path: str, delimiter: str = "\t") -> np.ndarray:
    r"""
    Convert a lower triangle matrix from a file to a symmetric square matrix.

    Parameters
    ----------
    file_path : str
        Path to the input text file containing the lower triangle matrix.
    delimiter : str, optional
        Delimiter used in the input file, by default "\t".

    Returns
    -------
    np.ndarray
        Symmetric square matrix with diagonal filled with 1.

    Raises
    ------
    ValueError
        If the input file is empty or does not contain a valid lower triangle matrix.
    FileNotFoundError
        If the specified file does not exist.

    Notes
    -----
    This function assumes that the input file contains a valid lower triangle matrix
    with each row on a new line and elements separated by the specified delimiter.
    For .npz files, it loads the first array key in the file.

    Examples
    --------
    >>> # Assuming 'lower_triangle.txt' contains:
    >>> # 1.0
    >>> # 0.1 1.0
    >>> # 0.2 0.4 1.0
    >>> # 0.3 0.5 0.6 1.0
    >>> matrix = load_ld_matrix('lower_triangle.txt')
    >>> print(matrix)
    array([[1.  , 0.1 , 0.2 , 0.3 ],
            [0.1 , 1.  , 0.4 , 0.5 ],
            [0.2 , 0.4 , 1.  , 0.6 ],
            [0.3 , 0.5 , 0.6 , 1.  ]])
    """
    if file_path.endswith(".npz"):
        ld_file_key = np.load(file_path).files[0]
        return np.load(file_path)[ld_file_key].astype(np.float32)
    lower_triangle = read_lower_triangle(file_path, delimiter)

    # Create the symmetric matrix
    symmetric_matrix = lower_triangle + lower_triangle.T

    # Fill the diagonal with 1
    np.fill_diagonal(symmetric_matrix, 1)

    # convert to float32
    symmetric_matrix = symmetric_matrix.astype(np.float32)
    return symmetric_matrix


def load_ld_map(map_path: str, delimiter: str = "\t") -> pd.DataFrame:
    r"""
    Read Variant IDs from a file.

    Parameters
    ----------
    map_path : str
        Path to the input text file containing the Variant IDs.
    delimiter : str, optional
        Delimiter used in the input file, by default "\t".

    Returns
    -------
    pd.DataFrame
        DataFrame containing the Variant IDs with columns CHR, BP, A1, A2, and SNPID.

    Raises
    ------
    ValueError
        If the input file is empty or does not contain the required columns.

    Notes
    -----
    This function assumes that the input file contains the required columns:

    - Chromosome (CHR)
    - Base pair position (BP)
    - Allele 1 (A1)
    - Allele 2 (A2)

    The function performs data cleaning including:

    - Converting chromosome and position to appropriate types
    - Validating alleles are valid DNA bases (A, C, G, T)
    - Removing variants where A1 == A2
    - Creating unique SNPID identifiers

    Examples
    --------
    >>> # Create sample map file
    >>> contents = "CHR\\tBP\\tA1\\tA2\\n1\\t1000\\tA\\tG\\n1\\t2000\\tC\\tT\\n2\\t3000\\tT\\tC"
    >>> with open('map.txt', 'w') as file:
    ...     file.write(contents)
    >>> df = load_ld_map('map.txt')
    >>> print(df)
        SNPID       CHR    BP A1 A2
    0   1-1000-A-G    1  1000  A  G
    1   1-2000-C-T    1  2000  C  T
    2   2-3000-C-T    2  3000  T  C
    """
    # TODO: use REF/ALT instead of A1/A2
    map_df = pd.read_csv(map_path, sep=delimiter)
    missing_cols = [col for col in ColName.map_cols if col not in map_df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in the input file: {missing_cols}")
    outdf = munge_chr(map_df)
    outdf = munge_bp(outdf)
    for col in [ColName.A1, ColName.A2]:
        pre_n = outdf.shape[0]
        outdf = outdf[outdf[col].notnull()]
        outdf[col] = outdf[col].astype(str).str.upper()
        outdf = outdf[outdf[col].str.match(r"^[ACGT]+$")]
        after_n = outdf.shape[0]
        logger.debug(f"Remove {pre_n - after_n} rows because of invalid {col}.")
    outdf = outdf[outdf[ColName.A1] != outdf[ColName.A2]]
    outdf = make_SNPID_unique(
        outdf, col_ea=ColName.A1, col_nea=ColName.A2, remove_duplicates=False
    )
    outdf.reset_index(drop=True, inplace=True)
    # TODO: check if allele frequency is available
    return outdf


def sort_alleles(ld: LDMatrix) -> LDMatrix:
    """
    Sort alleles in the LD map in alphabetical order. Change the sign of the LD matrix if the alleles are swapped.

    Parameters
    ----------
    ld : LDMatrix
        LDMatrix object containing the Variant IDs and the LD matrix.

    Returns
    -------
    LDMatrix
        LDMatrix object containing the Variant IDs and the LD matrix with alleles sorted.

    Notes
    -----
    This function ensures consistent allele ordering by:

    1. Sorting alleles alphabetically (A1 <= A2)
    2. Flipping the sign of LD correlations for variants where alleles were swapped
    3. Maintaining diagonal elements as 1.0

    This is important for consistent merging across different datasets.

    Examples
    --------
    >>> map_df = pd.DataFrame({
    ...     'SNPID': ['1-1000-A-G', '1-2000-C-T'],
    ...     'CHR': [1, 1],
    ...     'BP': [1000, 2000],
    ...     'A1': ['A', 'T'],
    ...     'A2': ['G', 'C']
    ... })
    >>> r_matrix = np.array([[1. , 0.1],
    ...                      [0.1, 1. ]])
    >>> ld = LDMatrix(map_df, r_matrix)
    >>> sorted_ld = sort_alleles(ld)
    >>> print(sorted_ld.map)
        SNPID       CHR    BP A1 A2
    0   1-1000-A-G    1  1000  A  G
    1   1-2000-C-T    1  2000  C  T
    >>> print(sorted_ld.r)
    array([[ 1. , -0.1],
            [-0.1,  1. ]])
    """
    ld_df = ld.r.copy()
    ld_map = ld.map.copy()
    ld_map[["sort_a1", "sort_a2"]] = np.sort(ld_map[[ColName.A1, ColName.A2]], axis=1)
    swapped_index = ld_map[ld_map[ColName.A1] != ld_map["sort_a1"]].index
    # Change the sign of the rows and columns the LD matrix if the alleles are swapped
    ld_df[swapped_index] *= -1
    ld_df[:, swapped_index] *= -1
    np.fill_diagonal(ld_df, 1)

    ld_map[ColName.A1] = ld_map["sort_a1"]
    ld_map[ColName.A2] = ld_map["sort_a2"]
    ld_map.drop(columns=["sort_a1", "sort_a2"], inplace=True)
    return LDMatrix(ld_map, ld_df)


def load_ld(
    ld_path: str, map_path: str, delimiter: str = "\t", if_sort_alleles: bool = True
) -> LDMatrix:
    r"""
    Read LD matrices and Variant IDs from files. Pair each matrix with its corresponding Variant IDs.

    Parameters
    ----------
    ld_path : str
        Path to the input text file containing the lower triangle matrix or .npz file.
    map_path : str
        Path to the input text file containing the Variant IDs.
    delimiter : str, optional
        Delimiter used in the input file, by default "\t".
    if_sort_alleles : bool, optional
        Sort alleles in the LD map in alphabetical order and change the sign of the
        LD matrix if the alleles are swapped, by default True.

    Returns
    -------
    LDMatrix
        Object containing the LD matrix and the Variant IDs.

    Raises
    ------
    ValueError
        If the number of variants in the map file does not match the number of rows in the LD matrix.

    Notes
    -----
    Future enhancements planned:

    - Support for npz files (partially implemented)
    - Support for plink bin4 format
    - Support for ldstore bcor format

    The function validates that the LD matrix and map file have consistent dimensions
    and optionally sorts alleles for consistent representation.

    Examples
    --------
    >>> ld_matrix = load_ld('data.ld', 'data.ldmap')
    >>> print(f"Loaded LD matrix with {ld_matrix.r.shape[0]} variants")
    Loaded LD matrix with 1000 variants
    """
    ld_df = load_ld_matrix(ld_path, delimiter)
    logger.info(f"Loaded LD matrix with shape {ld_df.shape} from '{ld_path}'.")
    map_df = load_ld_map(map_path, delimiter)
    logger.info(f"Loaded map file with shape {map_df.shape} from '{map_path}'.")
    if ld_df.shape[0] != map_df.shape[0]:
        raise ValueError(
            "The number of variants in the map file does not match the number of rows in the LD matrix.\n"
            f"Number of variants in the map file: {map_df.shape[0]}, number of rows in the LD matrix: {ld_df.shape[0]}"
            f"ld_path: {ld_path}, map_path: {map_path}"
        )
    ld = LDMatrix(map_df, ld_df)
    if if_sort_alleles:
        ld = sort_alleles(ld)

    return ld
