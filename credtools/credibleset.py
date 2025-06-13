"""Credible Set functions."""

import json
import logging
from dataclasses import dataclass
from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

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
        """
        self._tool = tool
        self._parameters = parameters
        self._coverage = coverage
        self._n_cs = n_cs
        self._cs_sizes = cs_sizes
        self._lead_snps = lead_snps
        self._snps = snps
        self._pips = pips
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
        return CredibleSet(
            tool=self.tool,
            parameters=self.parameters,
            coverage=self.coverage,
            n_cs=self.n_cs,
            cs_sizes=self.cs_sizes,
            lead_snps=self.lead_snps,
            snps=self.snps,
            pips=self.pips,
        )

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
        }

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
    )
    return merged


def continuous_jaccard(dict1: Dict[str, float], dict2: Dict[str, float]) -> float:
    """
    Calculate modified Jaccard similarity for continuous values (PIP values).

    Formula: ∑min(xi,yi)/∑max(xi,yi) where xi, yi are PIP values or 0 if missing

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
