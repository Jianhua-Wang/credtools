"""Class for the input data of the fine-mapping analysis."""

import logging
import os
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix, load_ld
from credtools.sumstats import load_sumstats

logger = logging.getLogger("Locus")


class Locus:
    """
    Locus class to represent a genomic locus with associated summary statistics and linkage disequilibrium (LD) matrix.

    Parameters
    ----------
    popu : str
        Population code. e.g. "EUR". Choose from ["AFR", "AMR", "EAS", "EUR", "SAS"].
    cohort : str
        Cohort name.
    sample_size : int
        Sample size.
    sumstats : pd.DataFrame
        Summary statistics DataFrame.
    ld : Optional[LDMatrix], optional
        LD matrix, by default None.
    if_intersect : bool, optional
        Whether to intersect the LD matrix and summary statistics file, by default False.

    Attributes
    ----------
    original_sumstats : pd.DataFrame
        The original summary statistics file.
    sumstats : pd.DataFrame
        The processed summary statistics file.
    ld : LDMatrix
        The LD matrix object.
    chrom : int
        Chromosome.
    start : int
        Start position of the locus.
    end : int
        End position of the locus.
    n_snps : int
        Number of SNPs in the locus.
    prefix : str
        The prefix combining population and cohort.
    locus_id : str
        Unique identifier for the locus.
    is_matched : bool
        Whether the LD matrix and summary statistics file are matched.
    lambda_s : Optional[float]
        The estimated lambda_s parameter from estimate_s_rss function, None if not calculated.

    Notes
    -----
    If no LD matrix is provided, only ABF method can be used for fine-mapping.
    """

    def __init__(
        self,
        popu: str,
        cohort: str,
        sample_size: int,
        sumstats: pd.DataFrame,
        locus_start: int,
        locus_end: int,
        ld: Optional[LDMatrix] = None,
        if_intersect: bool = False,
    ) -> None:
        """
        Initialize the Locus object.

        Parameters
        ----------
        popu : str
            Population code. e.g. "EUR". Choose from ["AFR", "AMR", "EAS", "EUR", "SAS"].
        cohort : str
            Cohort name.
        sample_size : int
            Sample size.
        sumstats : pd.DataFrame
            Summary statistics DataFrame.
        locus_start : int
            Fixed start position for the locus.
        locus_end : int
            Fixed end position for the locus.
        ld : Optional[LDMatrix], optional
            LD matrix, by default None.
        if_intersect : bool, optional
            Whether to intersect the LD matrix and summary statistics file, by default False.

        Warnings
        --------
        If no LD matrix is provided, a warning is logged that only ABF method can be used.
        """
        self.sumstats = sumstats
        self._original_sumstats = self.sumstats.copy()
        self._popu = popu
        self._cohort = cohort
        self._sample_size = sample_size
        self._locus_start = locus_start
        self._locus_end = locus_end
        self.lambda_s = None
        if ld:
            self.ld = ld
            if if_intersect:
                inters = intersect_sumstat_ld(self)
                self.sumstats = inters.sumstats
                self.ld = inters.ld
        else:
            logger.warning("LD matrix and map file not found. Can only run ABF method.")
            self.ld = LDMatrix(pd.DataFrame(), np.array([]))

    @property
    def original_sumstats(self) -> pd.DataFrame:
        """Get the original sumstats file."""
        return self._original_sumstats

    @property
    def popu(self) -> str:
        """Get the population code."""
        return self._popu

    @property
    def cohort(self) -> str:
        """Get the cohort name."""
        return self._cohort

    @property
    def sample_size(self) -> int:
        """Get the sample size."""
        return self._sample_size

    @property
    def chrom(self) -> int:
        """Get the chromosome."""
        return self.sumstats[ColName.CHR].iloc[0]

    @property
    def start(self) -> int:
        """Get the start position."""
        return self._locus_start

    @property
    def end(self) -> int:
        """Get the end position."""
        return self._locus_end

    @property
    def n_snps(self) -> int:
        """Get the number of SNPs."""
        return len(self.sumstats)

    @property
    def prefix(self) -> str:
        """Get the prefix of the locus."""
        import hashlib

        # Check if this is a meta-analysis result (contains '+' in cohort)
        if "+" in self.cohort:
            # Create a hash of the cohort part
            cohort_hash = hashlib.md5(self.cohort.encode()).hexdigest()[:8]
            # Count number of cohorts for reference
            num_cohorts = len(self.cohort.split("+"))
            return f"{self.popu}_meta{num_cohorts}cohorts_{cohort_hash}"

        # For single cohort, use original format
        return f"{self.popu}_{self.cohort}"

    @property
    def locus_id(self) -> str:
        """Get the locus ID."""
        return f"{self.popu}_{self.cohort}_chr{self.chrom}:{self.start}-{self.end}"

    @property
    def is_matched(self) -> bool:
        """Check if the LD matrix and sumstats file are matched."""
        # check the order of SNPID in the LD matrix and the sumstats file are the exact same
        if self.ld is None:
            return False
        return self.ld.map[ColName.SNPID].equals(self.sumstats[ColName.SNPID])

    def __repr__(self) -> str:
        """
        Return a string representation of the Locus object.

        Returns
        -------
        str
            String representation of the Locus object.
        """
        return f"Locus(popu={self.popu}, cohort={self.cohort}, sample_size={self.sample_size}, chr={self.chrom}, start={self.start}, end={self.end}, sumstats={self.sumstats.shape}, ld={self.ld.r.shape})"

    def copy(self) -> "Locus":
        """
        Copy the Locus object.

        Returns
        -------
        Locus
            A copy of the Locus object.
        """
        new_locus = Locus(
            self.popu,
            self.cohort,
            self.sample_size,
            self.sumstats.copy(),
            self._locus_start,
            self._locus_end,
            self.ld.copy(),
            if_intersect=False,
        )
        new_locus.lambda_s = self.lambda_s
        return new_locus


class LocusSet:
    """
    LocusSet class to represent a set of genomic loci.

    Parameters
    ----------
    loci : List[Locus]
        List of Locus objects.

    Attributes
    ----------
    loci : List[Locus]
        List of Locus objects.
    n_loci : int
        Number of loci.
    chrom : int
        Chromosome number.
    start : int
        Start position of the locus.
    end : int
        End position of the locus.
    locus_id : str
        Unique identifier for the locus.

    Raises
    ------
    ValueError
        If the chromosomes of the loci are not the same.
    """

    def __init__(self, loci: List[Locus]) -> None:
        """
        Initialize the LocusSet object.

        Parameters
        ----------
        loci : List[Locus]
            List of Locus objects.
        """
        self.loci = loci

    @property
    def n_loci(self) -> int:
        """Get the number of loci."""
        return len(self.loci)

    @property
    def chrom(self) -> int:
        """
        Get the chromosome.

        Returns
        -------
        int
            Chromosome number.

        Raises
        ------
        ValueError
            If the chromosomes of the loci are not the same.
        """
        chrom_list = [locus.chrom for locus in self.loci]
        if len(set(chrom_list)) > 1:
            raise ValueError("The chromosomes of the loci are not the same.")
        return chrom_list[0]

    @property
    def start(self) -> int:
        """Get the start position."""
        # All loci should have the same fixed start
        starts = [locus._locus_start for locus in self.loci]
        if not all(s == starts[0] for s in starts):
            raise ValueError("All loci in LocusSet must have the same start position")
        return starts[0]

    @property
    def end(self) -> int:
        """Get the end position."""
        # All loci should have the same fixed end
        ends = [locus._locus_end for locus in self.loci]
        if not all(e == ends[0] for e in ends):
            raise ValueError("All loci in LocusSet must have the same end position")
        return ends[0]

    @property
    def locus_id(self) -> str:
        """Get the locus ID."""
        return f"{self.chrom}:{self.start}-{self.end}"

    def __repr__(self) -> str:
        """
        Return a string representation of the LocusSet object.

        Returns
        -------
        str
            String representation of the LocusSet object.
        """
        return (
            f"LocusSet(\n n_loci={len(self.loci)}, chrom={self.chrom}, start={self.start}, end={self.end}, locus_id={self.locus_id} \n"
            + "\n".join([locus.__repr__() for locus in self.loci])
            + "\n"
            + ")"
        )

    def copy(self) -> "LocusSet":
        """
        Copy the LocusSet object.

        Returns
        -------
        LocusSet
            A copy of the LocusSet object.
        """
        return LocusSet([locus.copy() for locus in self.loci])


def intersect_sumstat_ld(locus: Locus) -> Locus:
    """
    Intersect the Variant IDs in the LD matrix and the sumstats file.

    Parameters
    ----------
    locus : Locus
        Locus object containing LD matrix and summary statistics.

    Returns
    -------
    Locus
        Locus object containing the intersected LD matrix and sumstats file.

    Raises
    ------
    ValueError
        If LD matrix not found or no common Variant IDs found between the LD matrix and the sumstats file.

    Warnings
    --------
    If only a few common Variant IDs are found (≤ 10), a warning is logged.

    Notes
    -----
    This function performs the following operations:

    1. Checks if LD matrix and summary statistics are already matched
    2. Finds common SNP IDs between LD matrix and summary statistics
    3. Subsets both datasets to common variants
    4. Reorders data to maintain consistency
    5. Returns a new Locus object with intersected data
    """
    if locus.ld is None:
        raise ValueError("LD matrix not found.")
    if locus.is_matched:
        logger.info("The LD matrix and sumstats file are matched.")
        return locus
    ldmap = locus.ld.map.copy()
    r = locus.ld.r.copy()
    sumstats = locus.sumstats.copy()
    sumstats = sumstats.sort_values([ColName.CHR, ColName.BP], ignore_index=True)
    intersec_sumstats = sumstats[
        sumstats[ColName.SNPID].isin(ldmap[ColName.SNPID])
    ].copy()
    intersec_variants = intersec_sumstats[ColName.SNPID].to_numpy()
    if len(intersec_variants) == 0:
        raise ValueError(
            f"No common Variant IDs found between the LD matrix and the sumstats file for locus {locus.locus_id}."
        )
    elif len(intersec_variants) <= 10:
        logger.warning(
            f"Only a few common Variant IDs found between the LD matrix and the sumstats file(<= 10) for locus {locus.locus_id}."
        )
    ldmap["idx"] = ldmap.index
    ldmap.set_index(ColName.SNPID, inplace=True, drop=False)
    ldmap = ldmap.loc[intersec_variants].copy()
    intersec_index = ldmap["idx"].to_numpy()
    r = r[intersec_index, :][:, intersec_index]
    intersec_sumstats.reset_index(drop=True, inplace=True)
    ldmap.drop("idx", axis=1, inplace=True)
    ldmap = ldmap.reset_index(drop=True)
    intersec_ld = LDMatrix(ldmap, r)
    logger.info(
        "Intersected the Variant IDs in the LD matrix and the sumstats file. "
        f"Number of common Variant IDs: {len(intersec_index)}"
    )
    return Locus(
        locus.popu,
        locus.cohort,
        locus.sample_size,
        intersec_sumstats,
        locus._locus_start,
        locus._locus_end,
        intersec_ld,
        if_intersect=False,
    )


def check_loci_info(loci_info: pd.DataFrame) -> pd.DataFrame:
    """
    Check and validate loci information DataFrame.

    Parameters
    ----------
    loci_info : pd.DataFrame
        DataFrame containing loci information.

    Returns
    -------
    pd.DataFrame
        Validated and type-corrected loci_info DataFrame.

    Raises
    ------
    ValueError
        If required columns are missing, data types are incorrect,
        or locus_id/boundary consistency checks fail.

    Notes
    -----
    This function performs the following checks:
    1. Ensures all required columns are present
    2. Validates and converts data types
    3. Checks that loci with same locus_id have same chr, start, end
    4. Validates chromosome, start, and end values
    """
    loci_info = loci_info.copy()

    # Check for required columns
    required_cols = [
        "prefix",
        "popu",
        "cohort",
        "sample_size",
        "chr",
        "start",
        "end",
        "locus_id",
    ]
    missing_cols = [col for col in required_cols if col not in loci_info.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Type checking and conversion
    try:
        # Convert numeric columns
        loci_info["sample_size"] = loci_info["sample_size"].astype(int)
        loci_info["chr"] = loci_info["chr"].astype(int)
        loci_info["start"] = loci_info["start"].astype(int)
        loci_info["end"] = loci_info["end"].astype(int)

        # Ensure string columns are strings
        loci_info["prefix"] = loci_info["prefix"].astype(str)
        loci_info["popu"] = loci_info["popu"].astype(str)
        loci_info["cohort"] = loci_info["cohort"].astype(str)
        loci_info["locus_id"] = loci_info["locus_id"].astype(str)

    except (ValueError, TypeError) as e:
        raise ValueError(f"Failed to convert data types: {e}")

    # Validate values
    if (loci_info["sample_size"] <= 0).any():
        raise ValueError("Sample size must be positive")

    if (loci_info["chr"] <= 0).any() or (loci_info["chr"] > 25).any():
        raise ValueError("Chromosome must be between 1 and 25")

    if (loci_info["start"] <= 0).any():
        raise ValueError("Start position must be positive")

    if (loci_info["end"] <= loci_info["start"]).any():
        raise ValueError("End position must be greater than start position")

    # Check for duplicates in popu+cohort+locus_id combination
    if loci_info.duplicated(subset=["popu", "cohort", "locus_id"]).any():
        raise ValueError("Each popu+cohort+locus_id combination must be unique")

    # Check consistency: same locus_id must have same chr, start, end
    locus_boundaries = loci_info.groupby("locus_id")[["chr", "start", "end"]].nunique()
    inconsistent_loci = locus_boundaries[(locus_boundaries > 1).any(axis=1)]

    if not inconsistent_loci.empty:
        raise ValueError(
            f"Inconsistent boundaries for locus_id(s): {inconsistent_loci.index.tolist()}. "
            "Each locus_id must have consistent chr, start, end values across all rows."
        )

    return loci_info


def intersect_loci(list_loci: List[Locus]) -> List[Locus]:
    """
    Intersect the Variant IDs in the LD matrices and the sumstats files of a list of Locus objects.

    Parameters
    ----------
    list_loci : List[Locus]
        List of Locus objects.

    Returns
    -------
    List[Locus]
        List of Locus objects containing the intersected LD matrices and sumstats files.

    Raises
    ------
    NotImplementedError
        This function is not yet implemented.

    Notes
    -----
    This function is planned to intersect variant IDs across multiple loci
    to ensure consistent variant sets for multi-ancestry analysis.
    """
    raise NotImplementedError(
        "Intersect the Variant IDs in the LD matrices and the sumstats files of a list of Locus objects."
    )


def load_locus(
    prefix: str,
    popu: str,
    cohort: str,
    sample_size: int,
    locus_start: int,
    locus_end: int,
    if_intersect: bool = False,
    calculate_lambda_s: bool = False,
    **kwargs: Any,
) -> Locus:
    """
    Load the input data of the fine-mapping analysis.

    Parameters
    ----------
    prefix : str
        Prefix of the input files.
    popu : str
        Population of the input data.
    cohort : str
        Cohort of the input data.
    sample_size : int
        Sample size of the input data.
    locus_start : int
        Fixed start position for the locus.
    locus_end : int
        Fixed end position for the locus.
    if_intersect : bool, optional
        Whether to intersect the input data with the LD matrix, by default False.
    calculate_lambda_s : bool, optional
        Whether to calculate lambda_s parameter using estimate_s_rss function, by default False.
    **kwargs : Any
        Additional keyword arguments passed to loading functions.

    Returns
    -------
    Locus
        Locus object containing the input data.

    Raises
    ------
    ValueError
        If the required input files are not found.

    Notes
    -----
    The function looks for files with the following patterns:

    - Summary statistics: {prefix}.sumstat or {prefix}.sumstats.gz
    - LD matrix: {prefix}.ld or {prefix}.ld.npz
    - LD map: {prefix}.ldmap or {prefix}.ldmap.gz

    All files are required for proper functioning.

    Examples
    --------
    >>> locus = load_locus('EUR_study1', 'EUR', 'study1', 50000)
    >>> print(f"Loaded locus with {locus.n_snps} SNPs")
    Loaded locus with 10000 SNPs
    """
    if os.path.exists(f"{prefix}.sumstat"):
        sumstats_path = f"{prefix}.sumstat"
    elif os.path.exists(f"{prefix}.sumstats.gz"):
        sumstats_path = f"{prefix}.sumstats.gz"
    else:
        raise ValueError("Sumstats file not found.")

    sumstats = load_sumstats(sumstats_path, if_sort_alleles=True, **kwargs)
    if os.path.exists(f"{prefix}.ld"):
        ld_path = f"{prefix}.ld"
    elif os.path.exists(f"{prefix}.ld.npz"):
        ld_path = f"{prefix}.ld.npz"
    else:
        raise ValueError("LD matrix file not found.")
    if os.path.exists(f"{prefix}.ldmap"):
        ldmap_path = f"{prefix}.ldmap"
    elif os.path.exists(f"{prefix}.ldmap.gz"):
        ldmap_path = f"{prefix}.ldmap.gz"
    else:
        raise ValueError("LD map file not found.")
    ld = load_ld(ld_path, ldmap_path, if_sort_alleles=True, **kwargs)

    locus = Locus(
        popu,
        cohort,
        sample_size,
        sumstats,
        locus_start,
        locus_end,
        ld=ld,
        if_intersect=if_intersect,
    )

    if calculate_lambda_s:
        try:
            # Import here to avoid circular imports
            from credtools.qc import estimate_s_rss

            locus.lambda_s = estimate_s_rss(locus)
            logger.info(
                f"Calculated lambda_s for locus {locus.locus_id}: {locus.lambda_s}"
            )
        except Exception as e:
            logger.warning(
                f"Failed to calculate lambda_s for locus {locus.locus_id}: {e}"
            )
            locus.lambda_s = None

    return locus


def load_locus_set(
    locus_info: pd.DataFrame,
    if_intersect: bool = False,
    calculate_lambda_s: bool = False,
    **kwargs: Any,
) -> LocusSet:
    """
    Load the input data of the fine-mapping analysis for multiple loci.

    Parameters
    ----------
    locus_info : pd.DataFrame
        DataFrame containing the locus information with required columns:
        ['prefix', 'popu', 'cohort', 'sample_size', 'chr', 'start', 'end', 'locus_id'].
    if_intersect : bool, optional
        Whether to intersect the input data with the LD matrix, by default False.
    calculate_lambda_s : bool, optional
        Whether to calculate lambda_s parameter using estimate_s_rss function, by default False.
    **kwargs : Any
        Additional keyword arguments passed to load_locus function.

    Returns
    -------
    LocusSet
        LocusSet object containing the input data.

    Raises
    ------
    ValueError
        If required columns are missing or if the combination of popu and cohort is not unique.

    Notes
    -----
    The locus_info DataFrame must contain the following columns:

    - prefix: File prefix for each locus
    - popu: Population code
    - cohort: Cohort name
    - sample_size: Sample size for the cohort
    - chr: Chromosome number
    - start: Start position of the locus
    - end: End position of the locus
    - locus_id: Locus identifier

    All rows must have the same chr, start, end, locus_id values (representing the same locus).

    Examples
    --------
    >>> locus_info = pd.DataFrame({
    ...     'prefix': ['EUR_study1', 'ASN_study2'],
    ...     'popu': ['EUR', 'ASN'],
    ...     'cohort': ['study1', 'study2'],
    ...     'sample_size': [50000, 30000]
    ... })
    >>> locus_set = load_locus_set(locus_info)
    >>> print(f"Loaded {locus_set.n_loci} loci")
    Loaded 2 loci
    """
    # Check and validate the locus_info DataFrame
    locus_info = check_loci_info(locus_info)

    # Check that all rows have the same chr, start, end (same locus)
    if len(locus_info["chr"].unique()) > 1:
        raise ValueError("All rows must have the same chromosome")
    if len(locus_info["start"].unique()) > 1:
        raise ValueError("All rows must have the same start position")
    if len(locus_info["end"].unique()) > 1:
        raise ValueError("All rows must have the same end position")
    if len(locus_info["locus_id"].unique()) > 1:
        raise ValueError("All rows must have the same locus_id")

    # Additional check for load_locus_set: popu+cohort must be unique within this single locus
    if locus_info.duplicated(subset=["popu", "cohort"]).any():
        raise ValueError(
            "Each popu+cohort combination must be unique within a single locus"
        )

    loci = []
    for i, row in locus_info.iterrows():
        loci.append(
            load_locus(
                row["prefix"],
                row["popu"],
                row["cohort"],
                row["sample_size"],
                int(row["start"]),
                int(row["end"]),
                if_intersect,
                calculate_lambda_s,
                **kwargs,
            )
        )
    return LocusSet(loci)
