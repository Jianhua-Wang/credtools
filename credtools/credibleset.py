"""Credible Set functions."""

import json
import logging
from dataclasses import dataclass
from itertools import combinations
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

if TYPE_CHECKING:
    from credtools.ldmatrix import LDMatrix

logger = logging.getLogger("CREDTOOLS")


class CredibleSet:
    """
    Class representing credible sets from one fine-mapping tool.

    Parameters
    ----------
    tool : str
        The name of the fine-mapping tool.
    parameters : Dict[str, Any]
        Additional parameters used by the fine-mapping tool.
    coverage : float
        The coverage of the credible sets.
    n_cs : int
        The number of credible sets.
    cs_sizes : List[int]
        Sizes of each credible set.
    lead_snps : List[str]
        List of lead SNPs.
    snps : List[List[str]]
        List of SNPs for each credible set.
    pips : pd.Series
        Posterior inclusion probabilities.

    Attributes
    ----------
    tool : str
        The name of the fine-mapping tool.
    n_cs : int
        The number of credible sets.
    coverage : float
        The coverage of the credible sets.
    lead_snps : List[str]
        List of lead SNPs.
    snps : List[List[str]]
        List of SNPs for each credible set.
    cs_sizes : List[int]
        Sizes of each credible set.
    pips : pd.Series
        Posterior inclusion probabilities.
    parameters : Dict[str, Any]
        Additional parameters used by the fine-mapping tool.
    """

    def __init__(
        self,
        tool: str,
        parameters: Dict[str, Any],
        coverage: float,
        n_cs: int,
        cs_sizes: List[int],
        lead_snps: List[str],
        snps: List[List[str]],
        pips: pd.Series,
        per_locus_results: Optional[Dict[str, "CredibleSet"]] = None,
        purity: Optional[List[Optional[float]]] = None,
    ) -> None:
        """
        Initialize CredibleSet object.

        Parameters
        ----------
        tool : str
            The name of the fine-mapping tool.
        parameters : Dict[str, Any]
            Additional parameters used by the fine-mapping tool.
        coverage : float
            The coverage of the credible sets.
        n_cs : int
            The number of credible sets.
        cs_sizes : List[int]
            Sizes of each credible set.
        lead_snps : List[str]
            List of lead SNPs.
        snps : List[List[str]]
            List of SNPs for each credible set.
        pips : pd.Series
            Posterior inclusion probabilities.
        per_locus_results : Optional[Dict[str, "CredibleSet"]], optional
            Mapping of locus identifiers to their individual credible set results.
        purity : Optional[List[Optional[float]]], optional
            List of purity values for each credible set. Purity is the minimum
            absolute LD R value between all SNP pairs in a credible set.
            None if LD matrix is not available.
        """
        self._tool = tool
        self._parameters = parameters
        self._coverage = coverage
        self._n_cs = n_cs
        self._cs_sizes = cs_sizes
        self._lead_snps = lead_snps
        self._snps = snps
        self._pips = pips
        self._per_locus_results: Dict[str, "CredibleSet"] = per_locus_results or {}
        self._purity = purity
        # TODO: add results data like, if it is converged, etc.

    @property
    def tool(self) -> str:
        """Get the tool name."""
        return self._tool

    @property
    def parameters(self) -> Dict[str, Any]:
        """Get the parameters."""
        return self._parameters

    @property
    def coverage(self) -> float:
        """Get the coverage."""
        # TODO: add actual coverage, as a list of coverage for each credible set
        return self._coverage

    @property
    def n_cs(self) -> int:
        """Get the number of credible sets."""
        return self._n_cs

    @property
    def cs_sizes(self) -> List[int]:
        """Get the sizes of each credible set."""
        return self._cs_sizes

    @property
    def lead_snps(self) -> List[str]:
        """Get the lead SNPs."""
        return self._lead_snps

    @property
    def snps(self) -> List[List[str]]:
        """Get the SNPs."""
        return self._snps

    @property
    def pips(self) -> pd.Series:
        """Get the PIPs."""
        return self._pips

    @property
    def per_locus_results(self) -> Dict[str, "CredibleSet"]:
        """Get per-locus credible set results."""
        return self._per_locus_results

    @property
    def purity(self) -> Optional[List[Optional[float]]]:
        """Get the purity values for each credible set."""
        return self._purity

    def set_per_locus_results(
        self, per_locus_results: Dict[str, "CredibleSet"]
    ) -> None:
        """Attach per-locus credible set results."""
        self._per_locus_results = per_locus_results

    def __repr__(self) -> str:
        """
        Return a string representation of the CredibleSet object.

        Returns
        -------
        str
            String representation of the CredibleSet object.
        """
        return (
            f"CredibleSet(\n  tool={self.tool}, coverage={self.coverage}, n_cs={self.n_cs}, cs_sizes={self.cs_sizes}, lead_snps={self.lead_snps},"
            + f"\n  Parameters: {json.dumps(self.parameters)}\n)"
        )

    def copy(self) -> "CredibleSet":
        """
        Copy the CredibleSet object.

        Returns
        -------
        CredibleSet
            A copy of the CredibleSet object.
        """
        copied = CredibleSet(
            tool=self.tool,
            parameters=dict(self.parameters),
            coverage=self.coverage,
            n_cs=self.n_cs,
            cs_sizes=self.cs_sizes.copy(),
            lead_snps=self.lead_snps.copy(),
            snps=[list(snp) for snp in self.snps],
            pips=self.pips.copy(),
            purity=self.purity.copy() if self.purity is not None else None,
        )
        if self.per_locus_results:
            per_locus_copy = {}
            for key, value in self.per_locus_results.items():
                if value is self:
                    per_locus_copy[key] = copied
                else:
                    per_locus_copy[key] = value.copy()
            copied.set_per_locus_results(per_locus_copy)
        return copied

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for TOML storage (excluding pips).

        Returns
        -------
        Dict[str, Any]
            A dictionary representation of the CredibleSet excluding pips.
        """
        return {
            "tool": self.tool,
            "n_cs": self.n_cs,
            "coverage": self.coverage,
            "lead_snps": self.lead_snps,
            "snps": self.snps,
            "cs_sizes": self.cs_sizes,
            "parameters": self.parameters,
            "purity": self.purity,
        }

    def create_enhanced_pips_df(self, locus_set) -> pd.DataFrame:
        """
        Create DataFrame with PIPs and full sumstats information.

        Parameters
        ----------
        locus_set : LocusSet
            The locus set containing locus data.

        Returns
        -------
        pd.DataFrame
            DataFrame containing full sumstats, PIPs, R2, and credible set assignments.
        """
        from credtools.constants import ColName
        from credtools.qc import intersect_sumstat_ld

        # Collect all unique SNPIDs from PIPs
        all_snpids = self.pips.index.tolist()

        # Initialize the result DataFrame with SNPIDs
        result_df = pd.DataFrame({ColName.SNPID: all_snpids})

        # Process based on number of loci
        if locus_set.n_loci == 1:
            # Single locus case - simpler column names
            locus = locus_set.loci[0]

            # Make sure we have matched sumstats and LD
            locus_copy = locus.copy()
            locus_copy = intersect_sumstat_ld(locus_copy)

            # Merge with sumstats
            sumstats_cols = [
                ColName.SNPID,
                ColName.CHR,
                ColName.BP,
                ColName.RSID,
                ColName.EA,
                ColName.NEA,
                ColName.EAF,
                ColName.MAF,
                ColName.BETA,
                ColName.SE,
                ColName.P,
            ]

            # Get available columns from sumstats
            available_cols = [
                col for col in sumstats_cols if col in locus_copy.sumstats.columns
            ]
            result_df = result_df.merge(
                locus_copy.sumstats[available_cols], on=ColName.SNPID, how="left"
            )

            # Calculate R2 (squared correlation with lead SNP)
            if locus_copy.ld is not None and len(locus_copy.sumstats) > 0:
                # Find lead SNP (lowest p-value)
                lead_idx = locus_copy.sumstats[ColName.P].idxmin()
                # Calculate R2 for all SNPs
                r2_values = locus_copy.ld.r[lead_idx] ** 2
                # Map R2 values to SNPIDs
                snpid_to_r2 = dict(zip(locus_copy.sumstats[ColName.SNPID], r2_values))
                result_df["R2"] = result_df[ColName.SNPID].map(snpid_to_r2)
            else:
                result_df["R2"] = np.nan

        else:
            # Multiple loci case - prefixed column names
            # First, add common columns that don't need prefix
            first_locus = locus_set.loci[0]
            common_cols = [
                ColName.CHR,
                ColName.BP,
                ColName.RSID,
                ColName.EA,
                ColName.NEA,
            ]
            available_common = [
                col for col in common_cols if col in first_locus.sumstats.columns
            ]

            # Use the first locus for common columns
            if available_common:
                result_df = result_df.merge(
                    first_locus.sumstats[[ColName.SNPID] + available_common],
                    on=ColName.SNPID,
                    how="left",
                )

            # Add locus-specific columns with prefixes
            for locus in locus_set.loci:
                prefix = f"{locus.popu}_{locus.cohort}_"

                # Make sure we have matched sumstats and LD
                locus_copy = locus.copy()
                locus_copy = intersect_sumstat_ld(locus_copy)

                # Columns to add with prefix
                locus_cols = [
                    ColName.EAF,
                    ColName.MAF,
                    ColName.BETA,
                    ColName.SE,
                    ColName.P,
                ]

                for col in locus_cols:
                    if col in locus_copy.sumstats.columns:
                        col_data = locus_copy.sumstats[[ColName.SNPID, col]].copy()
                        col_data.rename(columns={col: f"{prefix}{col}"}, inplace=True)
                        result_df = result_df.merge(
                            col_data, on=ColName.SNPID, how="left"
                        )

                # Calculate R2
                if locus_copy.ld is not None and len(locus_copy.sumstats) > 0:
                    lead_idx = locus_copy.sumstats[ColName.P].idxmin()
                    r2_values = locus_copy.ld.r[lead_idx] ** 2
                    snpid_to_r2 = dict(
                        zip(locus_copy.sumstats[ColName.SNPID], r2_values)
                    )
                    result_df[f"{prefix}R2"] = result_df[ColName.SNPID].map(snpid_to_r2)
                else:
                    result_df[f"{prefix}R2"] = np.nan

                # Add per-locus PIP and CRED columns when available
                if self.per_locus_results:
                    locus_creds = self.per_locus_results.get(locus.locus_id)
                    if locus_creds is not None:
                        pip_col = f"{prefix}PIP"
                        result_df[pip_col] = (
                            result_df[ColName.SNPID]
                            .map(locus_creds.pips.to_dict())
                            .fillna(0.0)
                        )
                        cred_col = f"{prefix}CRED"
                        result_df[cred_col] = 0
                        for cs_idx, snp_list in enumerate(locus_creds.snps, 1):
                            mask = result_df[ColName.SNPID].isin(snp_list)
                            result_df.loc[mask, cred_col] = cs_idx

        # Add credible set assignments (CRED column)
        result_df["CRED"] = 0  # Default: not in any credible set
        for cs_idx, snp_list in enumerate(self.snps, 1):
            mask = result_df[ColName.SNPID].isin(snp_list)
            result_df.loc[mask, "CRED"] = cs_idx

        # Add PIP column
        result_df["PIP"] = result_df[ColName.SNPID].map(self.pips.to_dict()).fillna(0)

        # Sort by PIP descending
        result_df = result_df.sort_values("PIP", ascending=False)

        return result_df

    @classmethod
    def from_dict(cls, data: Dict[str, Any], pips: pd.Series) -> "CredibleSet":
        """
        Create CredibleSet from dictionary and pips.

        Parameters
        ----------
        data : Dict[str, Any]
            A dictionary containing the data to initialize the CredibleSet.
        pips : pd.Series
            Posterior inclusion probabilities.

        Returns
        -------
        CredibleSet
            An instance of CredibleSet initialized with the provided data and pips.
        """
        return cls(
            tool=data["tool"],
            parameters=data["parameters"],
            coverage=data["coverage"],
            n_cs=data["n_cs"],
            cs_sizes=data["cs_sizes"],
            lead_snps=data["lead_snps"],
            snps=data["snps"],
            pips=pips,
            purity=data.get("purity"),
        )


def combine_pips(pips: List[pd.Series], method: str = "max") -> pd.Series:
    """
    Combine PIPs from multiple tools.

    Parameters
    ----------
    pips : List[pd.Series]
        List of PIPs from multiple tools.
    method : str, optional
        Method to combine PIPs, by default "max".
        Options: "max", "min", "mean", "meta".
        When "meta" is selected, the method will use the formula:
        PIP_meta = 1 - prod(1 - PIP_i), where i is the index of tools,
        PIP_i = 0 when the SNP is not in the credible set of the tool.
        When "max", "min", "mean" is selected, the SNP not in the credible set
        will be excluded from the calculation.

    Returns
    -------
    pd.Series
        Combined PIPs.

    Raises
    ------
    ValueError
        If the method is not supported.
    """
    logger.info(f"Combining PIPs using method: {method}")
    pip_df = pd.DataFrame(pips).T
    pip_df = pip_df.fillna(0)
    if method == "meta":
        merged = 1 - np.prod(1 - pip_df, axis=1)
    elif method == "max":
        merged = pip_df.max(axis=1)
    elif method == "min":
        merged = pip_df.min(axis=1)
    elif method == "mean":
        merged = pip_df.mean(axis=1)
    else:
        raise ValueError(f"Method {method} is not supported.")
    return merged


def combine_creds(
    creds: List[CredibleSet],
    combine_cred: str = "union",
    combine_pip: str = "max",
    jaccard_threshold: float = 0.1,
    ld_matrices: Optional[List["LDMatrix"]] = None,
    min_purity: float = 0.0,
) -> CredibleSet:
    """
    Combine credible sets from multiple tools.

    Parameters
    ----------
    creds : List[CredibleSet]
        List of credible sets from multiple tools.
    combine_cred : str, optional
        Method to combine credible sets, by default "union".
        Options: "union", "intersection", "cluster".

        - "union": Union of all credible sets to form a merged credible set.
        - "intersection": First merge the credible sets from the same tool,
            then take the intersection of all merged credible sets.
            No credible set will be returned if no common SNPs found.
        - "cluster": Merge credible sets with Jaccard index > jaccard_threshold.
    combine_pip : str, optional
        Method to combine PIPs, by default "max".
        Options: "max", "min", "mean", "meta".

        - "meta": PIP_meta = 1 - prod(1 - PIP_i), where i is the index of tools,
            PIP_i = 0 when the SNP is not in the credible set of the tool.
        - "max": Maximum PIP value for each SNP across all tools.
        - "min": Minimum PIP value for each SNP across all tools.
        - "mean": Mean PIP value for each SNP across all tools.
    jaccard_threshold : float, optional
        Jaccard index threshold for the "cluster" method, by default 0.1.
    ld_matrices : Optional[List[LDMatrix]], optional
        List of LD matrices for purity calculation, by default None.
        If provided, purity will be calculated for merged credible sets using
        multi-ancestry approach (element-wise max across populations).
        If None, purity will not be calculated for the merged credible sets.
    min_purity : float, optional
        Minimum purity threshold for filtering credible sets, by default 0.0.
        After combining credible sets, only those with purity >= min_purity will be kept.
        Purity is the minimum absolute LD R value between all SNP pairs in a credible set.
        Set to 0.0 (default) for no filtering.

    Returns
    -------
    CredibleSet
        Combined credible set.

    Raises
    ------
    ValueError
        If the method is not supported.

    Notes
    -----
    'union' and 'intersection' methods will merge all credible sets into one.
    """
    paras = creds[0].parameters
    tool = creds[0].tool
    # filter out the creds with no credible set
    creds = [cred for cred in creds if cred.n_cs > 0]
    if len(creds) == 0:
        logger.warning("No credible sets found in the input list.")
        return CredibleSet(
            tool=tool,
            n_cs=0,
            coverage=0,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=pd.Series(),
            parameters=paras,
        )
    if len(creds) == 1:
        return creds[0]
    if combine_cred == "union":
        merged_snps_flat = []
        for cred in creds:
            snps = [i for snp in cred.snps for i in snp]
            merged_snps_flat.extend(snps)
        merged_snps = [list(set(merged_snps_flat))]
    elif combine_cred == "intersection":
        merged_snps_set = None
        for i, cred in enumerate(creds):
            snps = [item for snp in cred.snps for item in snp]
            if i == 0:
                merged_snps_set = set(snps)
            else:
                if merged_snps_set is not None:
                    merged_snps_set.intersection_update(set(snps))
        if merged_snps_set is None or len(merged_snps_set) == 0:
            logger.warning("No common SNPs found in the intersection of credible sets.")
            merged_snps = [[]]
        else:
            merged_snps = [list(merged_snps_set)]
    elif combine_cred == "cluster":
        cred_pips = []
        for cred in creds:
            cred_pip = [dict(cred.pips[cred.pips.index.isin(snp)]) for snp in cred.snps]
            cred_pips.append(cred_pip)
        merged_snps = cluster_cs(cred_pips, 1 - jaccard_threshold)
        paras["jaccard_threshold"] = jaccard_threshold
    else:
        raise ValueError(f"Method {combine_cred} is not supported.")
    merged_pips = combine_pips([cred.pips for cred in creds], combine_pip)
    paras["combine_cred"] = combine_cred
    paras["combine_pip"] = combine_pip

    # Calculate purity for merged credible sets if LD matrices are provided
    purity = None
    if ld_matrices is not None and len(ld_matrices) > 0:
        purity = []
        for snp_list in merged_snps:
            if len(snp_list) > 0:
                purity_val = calculate_cs_purity(ld_matrices, snp_list)
                purity.append(purity_val)
            else:
                purity.append(None)
        logger.info(f"Calculated purity for {len(purity)} merged credible sets: {purity}")

    merged = CredibleSet(
        tool=creds[0].tool,
        n_cs=len(merged_snps),
        coverage=creds[0].coverage,
        lead_snps=[
            str(merged_pips[merged_pips.index.isin(snp)].idxmax())
            for snp in merged_snps
        ],
        snps=merged_snps,
        cs_sizes=[len(i) for i in merged_snps],
        pips=merged_pips,
        parameters=paras,
        purity=purity,
    )

    # Apply purity filtering if requested
    if min_purity > 0:
        merged = filter_credset_by_purity(merged, min_purity=min_purity)

    return merged


def continuous_jaccard(dict1: Dict[str, float], dict2: Dict[str, float]) -> float:
    """
    Calculate modified Jaccard similarity for continuous values (PIP values).

    Formula: ∑min(xi,yi)/∑max(xi,yi) where xi, yi are PIP values or 0 if missing

    Citation: Yuan, K. et al. (2024) Nature Genetics https://doi.org/10.1038/s41588-024-01870-z.

    Parameters
    ----------
    dict1 : Dict[str, float]
        First dictionary with keys and PIP values (0-1).
    dict2 : Dict[str, float]
        Second dictionary with keys and PIP values (0-1).

    Returns
    -------
    float
        Modified Jaccard similarity index between 0 and 1.

    Raises
    ------
    ValueError
        If any values are not between 0 and 1.

    Examples
    --------
    >>> d1 = {'a': 0.8, 'b': 0.5}
    >>> d2 = {'b': 0.6, 'c': 0.3}
    >>> continuous_jaccard(d1, d2)
    0.5
    """
    # Validate input values
    for d in [dict1, dict2]:
        invalid_values = [v for v in d.values() if not (0 <= v <= 1)]
        if invalid_values:
            raise ValueError("All values must be between 0 and 1")

    # Get all keys
    all_keys = set(dict1.keys()).union(set(dict2.keys()))

    # Calculate sum of minimums and maximums
    sum_min = 0.0
    sum_max = 0.0

    for key in all_keys:
        val1 = dict1.get(key, 0.0)
        val2 = dict2.get(key, 0.0)
        sum_min += min(val1, val2)
        sum_max += max(val1, val2)

    return sum_min / sum_max if sum_max > 0 else 0.0


def create_similarity_matrix(
    dict_sets: List[List[Dict[str, float]]],
) -> Tuple[np.ndarray, List[Dict[str, float]]]:
    """
    Create a similarity matrix for all pairs of dictionaries across different sets.

    Parameters
    ----------
    dict_sets : List[List[Dict[str, float]]]
        List of m sets, where each set contains dictionaries with PIP values.

    Returns
    -------
    Tuple[np.ndarray, List[Dict[str, float]]]
        A tuple containing:
        - Similarity matrix (n_dicts x n_dicts)
        - Flattened list of dictionaries

    Examples
    --------
    >>> sets = [[{'a': 0.8, 'b': 0.5}], [{'b': 0.6, 'c': 0.3}]]
    >>> matrix, dicts = create_similarity_matrix(sets)
    """
    # Flatten all dictionaries while keeping track of their set membership
    all_dicts = []
    for dict_set in dict_sets:
        all_dicts.extend(dict_set)

    total_dicts = len(all_dicts)

    # Create similarity matrix
    similarity_matrix = np.zeros((total_dicts, total_dicts))

    # Calculate set membership ranges
    set_ranges = []
    current_idx = 0
    for dict_set in dict_sets:
        set_ranges.append((current_idx, current_idx + len(dict_set)))
        current_idx += len(dict_set)

    # Fill similarity matrix
    for i, j in combinations(range(total_dicts), 2):
        # Check if dictionaries are from the same set
        same_set = False
        for start, end in set_ranges:
            if start <= i < end and start <= j < end:
                same_set = True
                break

        if not same_set:
            similarity = continuous_jaccard(all_dicts[i], all_dicts[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity

    return similarity_matrix, all_dicts


def cluster_cs(
    dict_sets: List[List[Dict[str, float]]], threshold: float = 0.9
) -> List[List[str]]:
    """
    Cluster dictionaries from different sets based on continuous Jaccard similarity.

    Parameters
    ----------
    dict_sets : List[List[Dict[str, float]]]
        List of m sets, where each set contains dictionaries with PIP values.
    threshold : float, optional
        Clustering threshold, by default 0.9.

    Returns
    -------
    List[List[str]]
        List of merged clusters, where each cluster contains
        a list of unique SNP IDs from the dictionaries in that cluster.

    Raises
    ------
    ValueError
        If less than two sets of dictionaries are provided or if any set is empty.

    Examples
    --------
    >>> sets = [
    ...     [{'a': 0.8, 'b': 0.5}],
    ...     [{'b': 0.6, 'c': 0.3}]
    ... ]
    >>> clusters = cluster_cs(sets)
    """
    if len(dict_sets) < 2:
        raise ValueError("At least two sets of dictionaries are required")

    # Validate input
    for dict_set in dict_sets:
        if not dict_set:
            raise ValueError("Empty dictionary sets are not allowed")

    # Create similarity matrix
    similarity_matrix, all_dicts = create_similarity_matrix(dict_sets)

    # Convert similarity to distance (1 - similarity)
    distance_matrix = 1 - similarity_matrix

    # Perform hierarchical clustering
    condensed_dist = distance_matrix[np.triu_indices(len(distance_matrix), k=1)]

    if len(condensed_dist) == 0:
        logger.warning("No valid distances found for clustering")
        return [list(set(all_dicts[0].keys()))]

    linkage_matrix = linkage(condensed_dist, method="average")

    # Cut the dendrogram at the specified threshold
    clusters = fcluster(linkage_matrix, threshold, criterion="distance")

    # Group dictionaries by cluster and merge them
    cluster_groups: Dict[int, List[str]] = {}
    for idx, cluster_id in enumerate(clusters):
        if cluster_id not in cluster_groups:
            cluster_groups[cluster_id] = []

            # Merge dictionaries within cluster by merging keys (no PIP values) and removing duplicates
            current_dict = all_dicts[idx]
            cluster_groups[cluster_id].extend(current_dict.keys())

    return [
        list(set(cluster_groups[cluster_id])) for cluster_id in sorted(cluster_groups)
    ]


def calculate_cs_purity(
    ld: Union["LDMatrix", List["LDMatrix"]],
    cs_snp_ids: List[str],
) -> Optional[float]:
    """
    Calculate purity for a single credible set.

    Purity is defined as the minimum absolute LD R value between all pairs of
    SNPs in the credible set.

    For multiple LD matrices (multi-ancestry case), purity is calculated as:
    1. Extract CS submatrix from each LD matrix
    2. Take element-wise maximum of absolute values across all matrices
    3. Return the minimum value from the resulting meta-LD matrix

    This approach (similar to MultiSuSiE) ensures the credible set has high
    purity across all populations.

    Parameters
    ----------
    ld : LDMatrix or List[LDMatrix]
        LDMatrix object(s) containing both r matrix and map with SNPIDs.
        If a list is provided, meta-purity across all matrices is calculated.
    cs_snp_ids : List[str]
        List of SNPID strings in the credible set.

    Returns
    -------
    Optional[float]
        - If CS has only 1 SNP, returns 1.0
        - If CS has multiple SNPs, returns min(|R|) for all SNP pairs
        - For multiple LD matrices, returns min of element-wise max across matrices
        - If unable to calculate (e.g., SNPs not in LD matrix), returns None

    Examples
    --------
    >>> # Single LD matrix: CS with 3 SNPs having LD R values: 0.8, 0.9, 0.7
    >>> # Purity = min(|0.8|, |0.9|, |0.7|) = 0.7
    >>>
    >>> # Multiple LD matrices: same CS in EUR and AFR
    >>> # EUR: |R| values = [0.8, 0.9, 0.7]
    >>> # AFR: |R| values = [0.6, 0.85, 0.75]
    >>> # Meta |R| = max([0.8, 0.9, 0.7], [0.6, 0.85, 0.75]) = [0.8, 0.9, 0.75]
    >>> # Purity = min([0.8, 0.9, 0.75]) = 0.75
    """
    from credtools.constants import ColName

    if len(cs_snp_ids) == 1:
        return 1.0

    # Handle single LD matrix vs list of LD matrices
    if not isinstance(ld, list):
        ld_list = [ld]
    else:
        ld_list = ld

    if len(ld_list) == 0:
        return None

    # Create union of all CS SNPs across all LD matrices (MultiSuSiE approach)
    # This ensures all submatrices have the same dimensions
    union_snps = []
    for snpid in cs_snp_ids:
        # Check if SNP appears in at least one LD matrix
        for ld_matrix in ld_list:
            snpid_to_idx = {snpid: i for i, snpid in enumerate(ld_matrix.map[ColName.SNPID])}
            if snpid in snpid_to_idx:
                union_snps.append(snpid)
                break

    # Remove duplicates while preserving order
    seen = set()
    union_snps = [x for x in union_snps if not (x in seen or seen.add(x))]

    if len(union_snps) < 2:
        # Not enough SNPs found in any LD matrix
        return None

    # Create mapping from SNP ID to index in union set
    variant_to_index = {snpid: i for i, snpid in enumerate(union_snps)}
    n_union = len(union_snps)

    # Extract and expand CS submatrices from all LD matrices
    cs_submatrices = []
    for ld_matrix in ld_list:
        # Create SNPID to index mapping for this LD matrix
        snpid_to_idx = {snpid: i for i, snpid in enumerate(ld_matrix.map[ColName.SNPID])}

        # Initialize expanded LD matrix with zeros
        expand_ld = np.zeros((n_union, n_union), dtype=np.float32)

        # Find SNPs that exist in both union and this LD matrix
        present_snps = [snpid for snpid in union_snps if snpid in snpid_to_idx]

        if len(present_snps) >= 2:
            # Get indices in this LD matrix
            ld_indices = np.array([snpid_to_idx[snpid] for snpid in present_snps])
            # Get indices in union set
            union_indices = np.array([variant_to_index[snpid] for snpid in present_snps])

            # Extract submatrix for present SNPs from this LD matrix
            ld_submatrix = ld_matrix.r[np.ix_(ld_indices, ld_indices)]

            # Place LD values at correct positions in expanded matrix using meshgrid
            idx_i, idx_j = np.meshgrid(union_indices, union_indices)
            expand_ld[idx_i, idx_j] = ld_submatrix.astype(np.float32)

        # Set diagonal to 1 (for both present and missing SNPs)
        np.fill_diagonal(expand_ld, 1.0)

        cs_submatrices.append(expand_ld)

    if len(cs_submatrices) == 0:
        # No valid submatrices found
        return None

    # Calculate meta-purity across all LD matrices
    # Take element-wise maximum of absolute values (MultiSuSiE approach)
    abs_meta_R = np.abs(cs_submatrices[0])
    for submatrix in cs_submatrices[1:]:
        abs_meta_R = np.maximum(abs_meta_R, np.abs(submatrix))

    # Get upper triangle (excluding diagonal) and find minimum
    upper_tri_indices = np.triu_indices_from(abs_meta_R, k=1)
    r_values = abs_meta_R[upper_tri_indices]

    if len(r_values) == 0:
        return None

    return float(np.min(r_values))


def filter_credset_by_purity(
    credset: "CredibleSet",
    min_purity: float = 0.0,
) -> "CredibleSet":
    """
    Filter credible sets by purity threshold.

    Removes credible sets that do not meet the minimum purity requirement.
    Purity is defined as the minimum absolute LD R value between all pairs
    of SNPs in the credible set.

    Parameters
    ----------
    credset : CredibleSet
        CredibleSet object containing credible sets and their purity values.
    min_purity : float, optional
        Minimum purity threshold for filtering, by default 0.0.
        Credible sets with purity < min_purity will be removed.
        Set to 0.0 (default) for no filtering.

    Returns
    -------
    CredibleSet
        New CredibleSet object with only credible sets meeting purity threshold.
        If no credible sets pass filtering, returns empty CredibleSet (n_cs=0).

    Notes
    -----
    - If credset.purity is None or empty, no filtering is applied (returns original credset)
    - If min_purity <= 0, no filtering is applied (returns original credset)
    - Filtered credible sets maintain their original ordering
    - PIPs are preserved for all variants (not filtered)

    Examples
    --------
    >>> # Filter credible sets to keep only high-purity sets (purity >= 0.5)
    >>> filtered_cs = filter_credset_by_purity(credset, min_purity=0.5)
    >>> print(f"Original: {credset.n_cs} CS, Filtered: {filtered_cs.n_cs} CS")
    Original: 5 CS, Filtered: 3 CS

    >>> # No filtering (default)
    >>> same_cs = filter_credset_by_purity(credset, min_purity=0.0)
    >>> assert same_cs.n_cs == credset.n_cs
    """
    # No filtering if min_purity <= 0
    if min_purity <= 0:
        return credset

    # No filtering if purity values are not available
    if credset.purity is None or len(credset.purity) == 0:
        logger.warning(
            "Purity values not available for filtering. "
            "Returning original credible set without filtering."
        )
        return credset

    # No credible sets to filter
    if credset.n_cs == 0:
        return credset

    # Filter credible sets by purity threshold
    keep_indices = []
    for i, purity_val in enumerate(credset.purity):
        if purity_val is not None and purity_val >= min_purity:
            keep_indices.append(i)

    # If no credible sets pass filtering, return empty CredibleSet
    if len(keep_indices) == 0:
        logger.warning(
            f"No credible sets passed purity filtering (min_purity={min_purity}). "
            f"All {credset.n_cs} credible sets were filtered out."
        )
        return CredibleSet(
            tool=credset.tool,
            n_cs=0,
            coverage=credset.coverage,
            lead_snps=[],
            snps=[],
            cs_sizes=[],
            pips=credset.pips,
            parameters=credset.parameters,
            purity=[],
        )

    # Filter credible sets
    filtered_snps = [credset.snps[i] for i in keep_indices]
    filtered_lead_snps = [credset.lead_snps[i] for i in keep_indices]
    filtered_cs_sizes = [credset.cs_sizes[i] for i in keep_indices]
    filtered_purity = [credset.purity[i] for i in keep_indices]

    logger.info(
        f"Filtered credible sets by purity >= {min_purity}: "
        f"{credset.n_cs} → {len(keep_indices)} credible sets"
    )

    return CredibleSet(
        tool=credset.tool,
        n_cs=len(keep_indices),
        coverage=credset.coverage,
        lead_snps=filtered_lead_snps,
        snps=filtered_snps,
        cs_sizes=filtered_cs_sizes,
        pips=credset.pips,  # Keep all PIPs (not filtered)
        parameters=credset.parameters,
        purity=filtered_purity,
    )
