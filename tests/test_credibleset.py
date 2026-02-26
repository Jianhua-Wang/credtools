"""Tests for credtools.credibleset module."""

import numpy as np
import pandas as pd
import pytest

from credtools.constants import ColName
from credtools.credibleset import (
    CredibleSet,
    filter_credset_by_purity,
    generate_cs_summary,
)
from credtools.ldmatrix import LDMatrix
from credtools.locus import Locus, LocusSet


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_credset(
    n_cs=2,
    tool="susie",
    purity=None,
    per_locus_results=None,
):
    """Create a CredibleSet for testing."""
    snps = [[f"snp_{i}_{j}" for j in range(3)] for i in range(n_cs)]
    lead_snps = [s[0] for s in snps]
    cs_sizes = [len(s) for s in snps]
    all_snpids = [snp for cs in snps for snp in cs]
    pips = pd.Series(
        np.linspace(0.9, 0.1, len(all_snpids)),
        index=all_snpids,
    )
    return CredibleSet(
        tool=tool,
        parameters={"max_causal": 5, "coverage": 0.95},
        coverage=0.95,
        n_cs=n_cs,
        cs_sizes=cs_sizes,
        lead_snps=lead_snps,
        snps=snps,
        pips=pips,
        purity=purity,
        per_locus_results=per_locus_results,
    )


# ---------------------------------------------------------------------------
# TestCredibleSetCopy
# ---------------------------------------------------------------------------
class TestCredibleSetCopy:
    """Tests for CredibleSet.copy method."""

    def test_basic_copy(self):
        cs = _make_credset()
        copied = cs.copy()
        assert copied.tool == cs.tool
        assert copied.n_cs == cs.n_cs
        assert copied.cs_sizes == cs.cs_sizes
        assert copied.lead_snps == cs.lead_snps
        assert copied.snps == cs.snps
        assert copied.pips.equals(cs.pips)
        # Ensure it's a deep copy
        copied._snps[0].append("extra")
        assert "extra" not in cs.snps[0]

    def test_copy_with_per_locus_results(self):
        inner = _make_credset(n_cs=1, tool="abf")
        cs = _make_credset(per_locus_results={"locus_1": inner})
        copied = cs.copy()
        assert "locus_1" in copied.per_locus_results
        assert copied.per_locus_results["locus_1"].tool == "abf"

    def test_self_referencing_copy(self):
        cs = _make_credset()
        cs.set_per_locus_results({"self": cs})
        copied = cs.copy()
        # The self-referencing copy should point to the new copied object
        assert copied.per_locus_results["self"] is copied

    def test_copy_with_purity(self):
        cs = _make_credset(purity=[0.8, 0.6])
        copied = cs.copy()
        assert copied.purity == [0.8, 0.6]
        copied.purity[0] = 0.99
        assert cs.purity[0] == 0.8  # original unchanged


# ---------------------------------------------------------------------------
# TestCredibleSetToDict
# ---------------------------------------------------------------------------
class TestCredibleSetToDict:
    """Tests for CredibleSet.to_dict method."""

    def test_to_dict_keys(self):
        cs = _make_credset(purity=[0.9, 0.7])
        d = cs.to_dict()
        expected_keys = {
            "tool",
            "n_cs",
            "coverage",
            "lead_snps",
            "snps",
            "cs_sizes",
            "parameters",
            "purity",
        }
        assert set(d.keys()) == expected_keys

    def test_roundtrip(self):
        cs = _make_credset(purity=[0.9, 0.7])
        d = cs.to_dict()
        restored = CredibleSet.from_dict(d, cs.pips)
        assert restored.tool == cs.tool
        assert restored.n_cs == cs.n_cs
        assert restored.coverage == cs.coverage
        assert restored.lead_snps == cs.lead_snps
        assert restored.snps == cs.snps
        assert restored.purity == cs.purity

    def test_to_dict_without_purity(self):
        cs = _make_credset()
        d = cs.to_dict()
        assert d["purity"] is None


# ---------------------------------------------------------------------------
# TestFilterCredsetByPurity
# ---------------------------------------------------------------------------
class TestFilterCredsetByPurity:
    """Tests for filter_credset_by_purity function."""

    def test_purity_none_returns_original(self):
        cs = _make_credset()
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == cs.n_cs  # no filtering applied

    def test_min_purity_zero_returns_original(self):
        cs = _make_credset(purity=[0.3, 0.8])
        result = filter_credset_by_purity(cs, min_purity=0.0)
        assert result.n_cs == cs.n_cs

    def test_filter_low_purity(self):
        cs = _make_credset(purity=[0.3, 0.8])
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == 1
        assert result.purity == [0.8]

    def test_all_filtered_returns_empty(self):
        cs = _make_credset(purity=[0.1, 0.2])
        result = filter_credset_by_purity(cs, min_purity=0.9)
        assert result.n_cs == 0
        assert result.snps == []
        assert result.lead_snps == []

    def test_empty_purity_list_returns_original(self):
        cs = CredibleSet(
            tool="test",
            parameters={},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
            purity=[],
        )
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == 0


# ---------------------------------------------------------------------------
# TestGenerateCsSummary
# ---------------------------------------------------------------------------
class TestGenerateCsSummary:
    """Tests for generate_cs_summary function."""

    def _make_locus_set(self):
        """Create a simple LocusSet for testing."""
        n_snps = 5
        snpids = [f"1-{1000 + i * 100}-A-G" for i in range(n_snps)]
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n_snps,
                ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                ColName.EA: ["A"] * n_snps,
                ColName.NEA: ["G"] * n_snps,
                ColName.EAF: [0.3] * n_snps,
                ColName.BETA: [0.1] * n_snps,
                ColName.SE: [0.01] * n_snps,
                ColName.P: [1e-8] * n_snps,
            }
        )
        r = np.eye(n_snps)
        ld_map = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n_snps,
                ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                ColName.A1: ["A"] * n_snps,
                ColName.A2: ["G"] * n_snps,
            }
        )
        ld = LDMatrix(ld_map, r)
        locus = Locus("EUR", "test", 10000, sumstats, 1000, 2000, ld=ld)
        return LocusSet([locus])

    def test_empty_dataframe_returns_empty(self):
        ls = self._make_locus_set()
        result = generate_cs_summary(pd.DataFrame(), "locus_1", ls)
        assert result == []

    def test_single_cs(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1, 1, 1],
                "PIP": [0.8, 0.15, 0.05],
                "SNPID": ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        assert len(result) == 1
        row = result[0]
        assert row["locus_id"] == "locus_1"
        assert row["cs_id"] == 1
        assert row["cs_size"] == 3
        assert row["lead_snp"] == "1-1000-A-G"

    def test_multiple_cs(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1, 1, 2, 2],
                "PIP": [0.8, 0.15, 0.7, 0.3],
                "SNPID": ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G", "1-1300-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        assert len(result) == 2
        assert result[0]["cs_id"] == 1
        assert result[1]["cs_id"] == 2

    def test_pip_thresholds(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1, 1, 1, 1],
                "PIP": [0.95, 0.55, 0.15, 0.05],
                "SNPID": ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G", "1-1300-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        row = result[0]
        assert row["pip_01"] == 3  # >=0.1: 0.95, 0.55, 0.15
        assert row["pip_05"] == 2  # >=0.5: 0.95, 0.55
        assert row["pip_09"] == 1  # >=0.9: 0.95

    def test_required_fields(self):
        ls = self._make_locus_set()
        causal = pd.DataFrame(
            {
                "CRED": [1],
                "PIP": [0.8],
                "SNPID": ["1-1000-A-G"],
            }
        )
        result = generate_cs_summary(causal, "locus_1", ls)
        expected_keys = {
            "locus_id",
            "cs_id",
            "lead_snp",
            "cs_size",
            "pip_01",
            "pip_05",
            "pip_09",
            "purity",
        }
        assert set(result[0].keys()) == expected_keys


# ---------------------------------------------------------------------------
# Additional imports for new tests
# ---------------------------------------------------------------------------
from credtools.credibleset import (
    calculate_cs_purity,
    cluster_cs,
    combine_creds,
    combine_pips,
    continuous_jaccard,
    create_similarity_matrix,
)


# ---------------------------------------------------------------------------
# Helper: build an LDMatrix with known r values
# ---------------------------------------------------------------------------
def _make_ld_matrix(snpids, r_matrix):
    """Create an LDMatrix from SNP IDs and a correlation matrix.

    Parameters
    ----------
    snpids : list of str
        SNP identifier strings (format: chr-bp-A-G).
    r_matrix : np.ndarray
        Square correlation matrix matching the number of SNPs.

    Returns
    -------
    LDMatrix
        LDMatrix with a proper map DataFrame and *r* array.
    """
    n = len(snpids)
    ld_map = pd.DataFrame(
        {
            ColName.SNPID: snpids,
            ColName.CHR: [1] * n,
            ColName.BP: [1000 + i * 100 for i in range(n)],
            ColName.A1: ["A"] * n,
            ColName.A2: ["G"] * n,
        }
    )
    return LDMatrix(ld_map, np.array(r_matrix, dtype=float))


# ---------------------------------------------------------------------------
# TestCombinePips
# ---------------------------------------------------------------------------
class TestCombinePips:
    """Tests for combine_pips function."""

    def test_meta_method(self):
        """Meta method: PIP_meta = 1 - prod(1 - PIP_i)."""
        pip1 = pd.Series([0.8, 0.5], index=["snp_a", "snp_b"])
        pip2 = pd.Series([0.6, 0.4], index=["snp_a", "snp_b"])
        result = combine_pips([pip1, pip2], method="meta")
        # 1 - (1-0.8)*(1-0.6) = 1 - 0.2*0.4 = 0.92
        assert np.isclose(result["snp_a"], 0.92)
        # 1 - (1-0.5)*(1-0.4) = 1 - 0.5*0.6 = 0.70
        assert np.isclose(result["snp_b"], 0.70)

    def test_meta_method_missing_snp(self):
        """Meta treats missing SNP as PIP=0."""
        pip1 = pd.Series([0.8], index=["snp_a"])
        pip2 = pd.Series([0.6], index=["snp_b"])
        result = combine_pips([pip1, pip2], method="meta")
        # snp_a: 1 - (1-0.8)*(1-0) = 0.8
        assert np.isclose(result["snp_a"], 0.8)
        # snp_b: 1 - (1-0)*(1-0.6) = 0.6
        assert np.isclose(result["snp_b"], 0.6)

    def test_min_method(self):
        """Min method picks the minimum PIP per SNP."""
        pip1 = pd.Series([0.8, 0.5], index=["snp_a", "snp_b"])
        pip2 = pd.Series([0.6, 0.9], index=["snp_a", "snp_b"])
        result = combine_pips([pip1, pip2], method="min")
        assert np.isclose(result["snp_a"], 0.6)
        assert np.isclose(result["snp_b"], 0.5)

    def test_mean_method(self):
        """Mean method averages PIP per SNP."""
        pip1 = pd.Series([0.8, 0.4], index=["snp_a", "snp_b"])
        pip2 = pd.Series([0.6, 0.2], index=["snp_a", "snp_b"])
        result = combine_pips([pip1, pip2], method="mean")
        assert np.isclose(result["snp_a"], 0.7)
        assert np.isclose(result["snp_b"], 0.3)

    def test_max_method(self):
        """Max method (default) picks the maximum PIP per SNP."""
        pip1 = pd.Series([0.3, 0.9], index=["snp_a", "snp_b"])
        pip2 = pd.Series([0.7, 0.1], index=["snp_a", "snp_b"])
        result = combine_pips([pip1, pip2], method="max")
        assert np.isclose(result["snp_a"], 0.7)
        assert np.isclose(result["snp_b"], 0.9)

    def test_invalid_method_raises(self):
        """Invalid method name raises ValueError."""
        pip1 = pd.Series([0.5], index=["snp_a"])
        with pytest.raises(ValueError, match="not supported"):
            combine_pips([pip1], method="invalid")


# ---------------------------------------------------------------------------
# TestCombineCreds
# ---------------------------------------------------------------------------
class TestCombineCreds:
    """Tests for combine_creds function."""

    def _make_simple_cred(self, snps_lists, tool="susie", pips_dict=None):
        """Build a CredibleSet with explicit SNP lists."""
        all_snpids = [snp for cs in snps_lists for snp in cs]
        if pips_dict is not None:
            pips = pd.Series(pips_dict)
        else:
            pips = pd.Series(
                np.linspace(0.9, 0.1, len(all_snpids)),
                index=all_snpids,
            )
        return CredibleSet(
            tool=tool,
            parameters={"max_causal": 5, "coverage": 0.95},
            coverage=0.95,
            n_cs=len(snps_lists),
            cs_sizes=[len(s) for s in snps_lists],
            lead_snps=[s[0] for s in snps_lists],
            snps=snps_lists,
            pips=pips,
        )

    def test_intersection_with_common_snps(self):
        """Intersection keeps only SNPs present in all input creds."""
        cred1 = self._make_simple_cred(
            [["snpA", "snpB", "snpC"]],
            pips_dict={"snpA": 0.8, "snpB": 0.5, "snpC": 0.3},
        )
        cred2 = self._make_simple_cred(
            [["snpB", "snpC", "snpD"]],
            pips_dict={"snpB": 0.7, "snpC": 0.4, "snpD": 0.2},
        )
        result = combine_creds(
            [cred1, cred2], combine_cred="intersection", combine_pip="max"
        )
        # Intersection of {A,B,C} and {B,C,D} = {B,C}
        merged_snps_flat = set(result.snps[0])
        assert merged_snps_flat == {"snpB", "snpC"}

    def test_intersection_no_common_snps_raises(self):
        """Intersection with no overlap hits idxmax on empty Series (source limitation)."""
        cred1 = self._make_simple_cred(
            [["snpA"]],
            pips_dict={"snpA": 0.9},
        )
        cred2 = self._make_simple_cred(
            [["snpB"]],
            pips_dict={"snpB": 0.8},
        )
        # The current implementation raises ValueError when trying idxmax
        # on an empty merged_snps list. This tests the actual behaviour.
        with pytest.raises(ValueError):
            combine_creds(
                [cred1, cred2], combine_cred="intersection", combine_pip="max"
            )

    def test_cluster_method(self):
        """Cluster method groups credible sets by Jaccard similarity."""
        cred1 = self._make_simple_cred(
            [["snpA", "snpB"]],
            pips_dict={"snpA": 0.8, "snpB": 0.5},
        )
        cred2 = self._make_simple_cred(
            [["snpA", "snpB"]],
            pips_dict={"snpA": 0.7, "snpB": 0.6},
        )
        result = combine_creds(
            [cred1, cred2],
            combine_cred="cluster",
            combine_pip="max",
            jaccard_threshold=0.1,
        )
        # Both share snpA, snpB so they should cluster together
        all_merged = set()
        for s in result.snps:
            all_merged.update(s)
        assert "snpA" in all_merged
        assert "snpB" in all_merged

    def test_invalid_combine_cred_raises(self):
        """Invalid combine_cred value raises ValueError."""
        cred1 = self._make_simple_cred(
            [["snpA"]],
            pips_dict={"snpA": 0.9},
        )
        cred2 = self._make_simple_cred(
            [["snpB"]],
            pips_dict={"snpB": 0.8},
        )
        with pytest.raises(ValueError, match="not supported"):
            combine_creds([cred1, cred2], combine_cred="invalid")

    def test_all_creds_n_cs_zero_returns_empty(self):
        """All input creds with n_cs=0 returns empty CredibleSet."""
        empty_cred = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5, "coverage": 0.95},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
        )
        result = combine_creds([empty_cred, empty_cred])
        assert result.n_cs == 0
        assert result.snps == []

    def test_single_cred_after_filtering_returns_that_cred(self):
        """Only one cred with n_cs>0 among inputs returns it directly."""
        good_cred = self._make_simple_cred(
            [["snpA", "snpB"]],
            pips_dict={"snpA": 0.9, "snpB": 0.5},
        )
        empty_cred = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5, "coverage": 0.95},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
        )
        result = combine_creds([good_cred, empty_cred])
        assert result is good_cred

    def test_union_method(self):
        """Union merges all SNPs from all creds into a single set."""
        cred1 = self._make_simple_cred(
            [["snpA"]],
            pips_dict={"snpA": 0.9},
        )
        cred2 = self._make_simple_cred(
            [["snpB"]],
            pips_dict={"snpB": 0.8},
        )
        result = combine_creds([cred1, cred2], combine_cred="union", combine_pip="max")
        merged = set(result.snps[0])
        assert merged == {"snpA", "snpB"}


# ---------------------------------------------------------------------------
# TestContinuousJaccard
# ---------------------------------------------------------------------------
class TestContinuousJaccard:
    """Tests for continuous_jaccard function."""

    def test_identical_dicts(self):
        """Identical dictionaries yield Jaccard = 1.0."""
        d = {"a": 0.8, "b": 0.5}
        assert np.isclose(continuous_jaccard(d, d), 1.0)

    def test_completely_different_dicts(self):
        """Completely disjoint non-zero dicts: sum_min = 0."""
        d1 = {"a": 0.8}
        d2 = {"b": 0.6}
        # min(0.8,0)=0, min(0,0.6)=0 => sum_min=0
        assert np.isclose(continuous_jaccard(d1, d2), 0.0)

    def test_partial_overlap(self):
        """Partial overlap yields expected ratio."""
        d1 = {"a": 0.8, "b": 0.5}
        d2 = {"b": 0.6, "c": 0.3}
        # keys: a, b, c
        # a: min(0.8,0)=0, max(0.8,0)=0.8
        # b: min(0.5,0.6)=0.5, max(0.5,0.6)=0.6
        # c: min(0,0.3)=0, max(0,0.3)=0.3
        # sum_min=0.5, sum_max=1.7 => ~0.2941
        expected = 0.5 / 1.7
        assert np.isclose(continuous_jaccard(d1, d2), expected)

    def test_empty_dicts(self):
        """Two empty dicts yield Jaccard = 0.0 (sum_max=0 edge case)."""
        assert np.isclose(continuous_jaccard({}, {}), 0.0)

    def test_invalid_values_above_one_raises(self):
        """Values > 1 raise ValueError."""
        with pytest.raises(ValueError, match="between 0 and 1"):
            continuous_jaccard({"a": 1.5}, {"a": 0.5})

    def test_invalid_values_below_zero_raises(self):
        """Negative values raise ValueError."""
        with pytest.raises(ValueError, match="between 0 and 1"):
            continuous_jaccard({"a": -0.1}, {"a": 0.5})


# ---------------------------------------------------------------------------
# TestCreateSimilarityMatrix
# ---------------------------------------------------------------------------
class TestCreateSimilarityMatrix:
    """Tests for create_similarity_matrix function."""

    def test_basic_two_sets(self):
        """Two sets with one dict each produce a 2x2 matrix."""
        sets = [[{"a": 0.8, "b": 0.5}], [{"b": 0.6, "c": 0.3}]]
        matrix, all_dicts = create_similarity_matrix(sets)
        assert matrix.shape == (2, 2)
        # off-diagonal should be the jaccard of the two dicts
        expected_jaccard = continuous_jaccard(sets[0][0], sets[1][0])
        assert np.isclose(matrix[0, 1], expected_jaccard)
        assert np.isclose(matrix[1, 0], expected_jaccard)

    def test_same_set_pairs_are_zero(self):
        """Pairs from the same set must have similarity = 0."""
        sets = [[{"a": 0.8}, {"b": 0.5}], [{"a": 0.6}]]
        matrix, all_dicts = create_similarity_matrix(sets)
        # indices 0 and 1 belong to the same set
        assert matrix[0, 1] == 0.0
        assert matrix[1, 0] == 0.0
        # index 0 vs 2 (different sets) should be computed
        assert matrix.shape == (3, 3)

    def test_diagonal_is_zero(self):
        """Diagonal of the similarity matrix should be zero."""
        sets = [[{"a": 0.5}], [{"a": 0.5}]]
        matrix, _ = create_similarity_matrix(sets)
        assert matrix[0, 0] == 0.0
        assert matrix[1, 1] == 0.0


# ---------------------------------------------------------------------------
# TestClusterCS
# ---------------------------------------------------------------------------
class TestClusterCS:
    """Tests for cluster_cs function."""

    def test_basic_clustering_two_sets(self):
        """Two similar dict sets cluster into at least one group."""
        sets = [
            [{"snpA": 0.8, "snpB": 0.5}],
            [{"snpA": 0.7, "snpB": 0.6}],
        ]
        result = cluster_cs(sets, threshold=0.9)
        # Should produce list of lists of SNP IDs
        assert isinstance(result, list)
        assert len(result) >= 1
        # All original SNP IDs should appear
        all_snps = set()
        for cluster in result:
            all_snps.update(cluster)
        assert "snpA" in all_snps
        assert "snpB" in all_snps

    def test_less_than_two_sets_raises(self):
        """Less than two sets raises ValueError."""
        with pytest.raises(ValueError, match="At least two sets"):
            cluster_cs([[{"a": 0.5}]])

    def test_empty_dict_set_raises(self):
        """An empty inner list raises ValueError."""
        with pytest.raises(ValueError, match="Empty dictionary sets"):
            cluster_cs([[{"a": 0.5}], []])

    def test_dissimilar_sets_form_separate_clusters(self):
        """Very different dict sets form separate clusters."""
        sets = [
            [{"snpA": 0.9}],
            [{"snpZ": 0.9}],
        ]
        result = cluster_cs(sets, threshold=0.5)
        # Completely disjoint keys => jaccard=0 => distance=1 > 0.5
        # Should form 2 clusters
        assert len(result) == 2


# ---------------------------------------------------------------------------
# TestCalculateCsPurity
# ---------------------------------------------------------------------------
class TestCalculateCsPurity:
    """Tests for calculate_cs_purity function."""

    def test_single_snp_returns_one(self):
        """Single SNP in credible set returns purity = 1.0."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r = np.eye(2)
        ld = _make_ld_matrix(snpids, r)
        result = calculate_cs_purity(ld, ["1-1000-A-G"])
        assert result == 1.0

    def test_two_snps_known_ld(self):
        """Two SNPs with known R value: purity = |R|."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r = np.array([[1.0, 0.7], [0.7, 1.0]])
        ld = _make_ld_matrix(snpids, r)
        result = calculate_cs_purity(ld, snpids)
        assert np.isclose(result, 0.7)

    def test_three_snps_min_ld(self):
        """Three SNPs: purity is the minimum off-diagonal |R|."""
        snpids = ["1-1000-A-G", "1-1100-A-G", "1-1200-A-G"]
        r = np.array(
            [
                [1.0, 0.9, 0.6],
                [0.9, 1.0, 0.8],
                [0.6, 0.8, 1.0],
            ]
        )
        ld = _make_ld_matrix(snpids, r)
        result = calculate_cs_purity(ld, snpids)
        assert np.isclose(result, 0.6)

    def test_multiple_ld_matrices_multi_ancestry(self):
        """Multi-ancestry: element-wise max of |R| across populations."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r_eur = np.array([[1.0, 0.5], [0.5, 1.0]])
        r_afr = np.array([[1.0, 0.8], [0.8, 1.0]])
        ld_eur = _make_ld_matrix(snpids, r_eur)
        ld_afr = _make_ld_matrix(snpids, r_afr)
        result = calculate_cs_purity([ld_eur, ld_afr], snpids)
        # max(|0.5|, |0.8|) = 0.8
        assert np.isclose(result, 0.8)

    def test_snps_not_in_ld_returns_none(self):
        """Return None when SNPs not found in LD matrix."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r = np.eye(2)
        ld = _make_ld_matrix(snpids, r)
        result = calculate_cs_purity(ld, ["1-9000-A-G", "1-9100-A-G"])
        assert result is None

    def test_empty_ld_list_returns_none(self):
        """Empty LD list returns None."""
        result = calculate_cs_purity([], ["1-1000-A-G", "1-1100-A-G"])
        assert result is None

    def test_negative_r_values(self):
        """Negative R values: purity uses |R|."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r = np.array([[1.0, -0.65], [-0.65, 1.0]])
        ld = _make_ld_matrix(snpids, r)
        result = calculate_cs_purity(ld, snpids)
        assert np.isclose(result, 0.65)

    def test_one_snp_missing_from_ld(self):
        """If only one CS SNP is in LD, returns None (need >= 2)."""
        snpids = ["1-1000-A-G"]
        r = np.array([[1.0]])
        ld = _make_ld_matrix(snpids, r)
        # Ask for two SNPs, only one is in the LD matrix
        result = calculate_cs_purity(ld, ["1-1000-A-G", "1-9999-A-G"])
        assert result is None


# ---------------------------------------------------------------------------
# TestCreateEnhancedPipsDf
# ---------------------------------------------------------------------------
class TestCreateEnhancedPipsDf:
    """Tests for CredibleSet.create_enhanced_pips_df method (single locus)."""

    def _make_locus_and_credset(self, with_ld=True):
        """Create a single-locus LocusSet and matching CredibleSet."""
        n_snps = 4
        snpids = [f"1-{1000 + i * 100}-A-G" for i in range(n_snps)]
        sumstats = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n_snps,
                ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                ColName.EA: ["A"] * n_snps,
                ColName.NEA: ["G"] * n_snps,
                ColName.EAF: [0.3] * n_snps,
                ColName.BETA: [0.1 * (i + 1) for i in range(n_snps)],
                ColName.SE: [0.01] * n_snps,
                ColName.P: [1e-8, 1e-6, 1e-4, 1e-2],
            }
        )
        if with_ld:
            r = np.array(
                [
                    [1.0, 0.9, 0.5, 0.2],
                    [0.9, 1.0, 0.6, 0.3],
                    [0.5, 0.6, 1.0, 0.4],
                    [0.2, 0.3, 0.4, 1.0],
                ]
            )
            ld_map = pd.DataFrame(
                {
                    ColName.SNPID: snpids,
                    ColName.CHR: [1] * n_snps,
                    ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                    ColName.A1: ["A"] * n_snps,
                    ColName.A2: ["G"] * n_snps,
                }
            )
            ld = LDMatrix(ld_map, r)
        else:
            ld = None

        locus = Locus("EUR", "test", 10000, sumstats, 1000, 2000, ld=ld)
        locus_set = LocusSet([locus])

        pips = pd.Series(
            [0.8, 0.5, 0.15, 0.05],
            index=snpids,
        )
        cs = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5, "coverage": 0.95},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[3],
            lead_snps=[snpids[0]],
            snps=[snpids[:3]],
            pips=pips,
        )
        return locus_set, cs

    def test_single_locus_with_ld_has_r2(self):
        """Single locus with LD produces R2 column with correct values."""
        locus_set, cs = self._make_locus_and_credset(with_ld=True)
        df = cs.create_enhanced_pips_df(locus_set)
        assert "R2" in df.columns
        # Lead SNP (lowest P = 1e-8 at index 0) should have R2=1.0
        lead_row = df[df[ColName.SNPID] == "1-1000-A-G"]
        assert np.isclose(lead_row["R2"].values[0], 1.0)
        # Other SNPs should have R2 values based on LD with lead
        snp2_row = df[df[ColName.SNPID] == "1-1100-A-G"]
        assert np.isclose(snp2_row["R2"].values[0], 0.81, atol=0.01)

    def test_single_locus_without_ld_hits_key_error(self):
        """Single locus without LD: Locus creates empty LDMatrix lacking SNPID column.

        The current implementation raises KeyError when intersect_sumstat_ld
        tries to access SNPID from an empty LD map DataFrame.
        """
        locus_set, cs = self._make_locus_and_credset(with_ld=False)
        with pytest.raises(KeyError):
            cs.create_enhanced_pips_df(locus_set)

    def test_pip_and_cred_columns_present(self):
        """Output DataFrame has PIP and CRED columns."""
        locus_set, cs = self._make_locus_and_credset(with_ld=True)
        df = cs.create_enhanced_pips_df(locus_set)
        assert "PIP" in df.columns
        assert "CRED" in df.columns
        # SNPs in credible set should have CRED=1
        in_cs = df[df[ColName.SNPID].isin(cs.snps[0])]
        assert (in_cs["CRED"] == 1).all()
        # SNP not in credible set should have CRED=0
        out_cs = df[~df[ColName.SNPID].isin(cs.snps[0])]
        assert (out_cs["CRED"] == 0).all()

    def test_sorted_by_pip_descending(self):
        """Output is sorted by PIP in descending order."""
        locus_set, cs = self._make_locus_and_credset(with_ld=True)
        df = cs.create_enhanced_pips_df(locus_set)
        pips = df["PIP"].values
        assert all(pips[i] >= pips[i + 1] for i in range(len(pips) - 1))


# ---------------------------------------------------------------------------
# TestCreateEnhancedPipsDfMultiLocus
# ---------------------------------------------------------------------------
class TestCreateEnhancedPipsDfMultiLocus:
    """Tests for CredibleSet.create_enhanced_pips_df with multiple loci."""

    def _make_multi_locus_setup(self, with_per_locus=False):
        """Create a two-locus LocusSet and matching CredibleSet."""
        n_snps = 3
        snpids = [f"1-{1000 + i * 100}-A-G" for i in range(n_snps)]

        def _make_sumstats():
            return pd.DataFrame(
                {
                    ColName.SNPID: snpids,
                    ColName.CHR: [1] * n_snps,
                    ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                    ColName.EA: ["A"] * n_snps,
                    ColName.NEA: ["G"] * n_snps,
                    ColName.EAF: [0.3] * n_snps,
                    ColName.BETA: [0.1, 0.2, 0.05],
                    ColName.SE: [0.01] * n_snps,
                    ColName.P: [1e-8, 1e-6, 1e-3],
                }
            )

        r = np.array(
            [
                [1.0, 0.8, 0.3],
                [0.8, 1.0, 0.4],
                [0.3, 0.4, 1.0],
            ]
        )
        ld_map = pd.DataFrame(
            {
                ColName.SNPID: snpids,
                ColName.CHR: [1] * n_snps,
                ColName.BP: [1000 + i * 100 for i in range(n_snps)],
                ColName.A1: ["A"] * n_snps,
                ColName.A2: ["G"] * n_snps,
            }
        )
        ld1 = LDMatrix(ld_map.copy(), r.copy())
        ld2 = LDMatrix(ld_map.copy(), r.copy())

        locus1 = Locus("EUR", "cohort1", 10000, _make_sumstats(), 1000, 2000, ld=ld1)
        locus2 = Locus("AFR", "cohort2", 8000, _make_sumstats(), 1000, 2000, ld=ld2)
        locus_set = LocusSet([locus1, locus2])

        pips = pd.Series([0.7, 0.5, 0.1], index=snpids)
        per_locus = None
        if with_per_locus:
            per_locus_cs = CredibleSet(
                tool="susie",
                parameters={"max_causal": 5, "coverage": 0.95},
                coverage=0.95,
                n_cs=1,
                cs_sizes=[2],
                lead_snps=[snpids[0]],
                snps=[snpids[:2]],
                pips=pd.Series([0.6, 0.4, 0.05], index=snpids),
            )
            per_locus = {
                locus1.locus_id: per_locus_cs,
                locus2.locus_id: per_locus_cs,
            }

        cs = CredibleSet(
            tool="susie",
            parameters={"max_causal": 5, "coverage": 0.95},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[2],
            lead_snps=[snpids[0]],
            snps=[snpids[:2]],
            pips=pips,
            per_locus_results=per_locus,
        )
        return locus_set, cs

    def test_multi_locus_prefixed_columns(self):
        """Multiple loci produce prefixed locus-specific columns."""
        locus_set, cs = self._make_multi_locus_setup(with_per_locus=False)
        df = cs.create_enhanced_pips_df(locus_set)
        # Should have prefixed R2 columns
        assert "EUR_cohort1_R2" in df.columns
        assert "AFR_cohort2_R2" in df.columns
        # Should have prefixed BETA columns
        assert "EUR_cohort1_BETA" in df.columns
        assert "AFR_cohort2_BETA" in df.columns
        # Should have common columns without prefix
        assert ColName.CHR in df.columns
        assert ColName.BP in df.columns

    def test_multi_locus_with_per_locus_results(self):
        """Multiple loci with per_locus_results produce per-locus PIP/CRED columns."""
        locus_set, cs = self._make_multi_locus_setup(with_per_locus=True)
        df = cs.create_enhanced_pips_df(locus_set)
        assert "EUR_cohort1_PIP" in df.columns
        assert "EUR_cohort1_CRED" in df.columns
        assert "AFR_cohort2_PIP" in df.columns
        assert "AFR_cohort2_CRED" in df.columns

    def test_multi_locus_still_has_global_pip_and_cred(self):
        """Global PIP and CRED columns present in multi-locus output."""
        locus_set, cs = self._make_multi_locus_setup(with_per_locus=False)
        df = cs.create_enhanced_pips_df(locus_set)
        assert "PIP" in df.columns
        assert "CRED" in df.columns


# ---------------------------------------------------------------------------
# TestCombineCredsWithPurity
# ---------------------------------------------------------------------------
class TestCombineCredsWithPurity:
    """Tests for combine_creds with ld_matrices and min_purity parameters."""

    def _make_cred_with_snpids(self, snpids, pips_dict):
        """Build a CredibleSet from real SNPID-style identifiers."""
        return CredibleSet(
            tool="susie",
            parameters={"max_causal": 5, "coverage": 0.95},
            coverage=0.95,
            n_cs=1,
            cs_sizes=[len(snpids)],
            lead_snps=[snpids[0]],
            snps=[snpids],
            pips=pd.Series(pips_dict),
        )

    def test_union_with_ld_matrices_calculates_purity(self):
        """Combine with ld_matrices produces purity values."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r = np.array([[1.0, 0.7], [0.7, 1.0]])
        ld = _make_ld_matrix(snpids, r)

        cred1 = self._make_cred_with_snpids(
            snpids, {"1-1000-A-G": 0.8, "1-1100-A-G": 0.5}
        )
        cred2 = self._make_cred_with_snpids(
            snpids, {"1-1000-A-G": 0.7, "1-1100-A-G": 0.6}
        )
        result = combine_creds(
            [cred1, cred2],
            combine_cred="union",
            combine_pip="max",
            ld_matrices=[ld],
        )
        assert result.purity is not None
        assert len(result.purity) == result.n_cs
        assert np.isclose(result.purity[0], 0.7)

    def test_union_with_min_purity_filters(self):
        """Combine with min_purity filters out low-purity credible sets."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        r = np.array([[1.0, 0.3], [0.3, 1.0]])
        ld = _make_ld_matrix(snpids, r)

        cred1 = self._make_cred_with_snpids(
            snpids, {"1-1000-A-G": 0.8, "1-1100-A-G": 0.5}
        )
        cred2 = self._make_cred_with_snpids(
            snpids, {"1-1000-A-G": 0.7, "1-1100-A-G": 0.6}
        )
        result = combine_creds(
            [cred1, cred2],
            combine_cred="union",
            combine_pip="max",
            ld_matrices=[ld],
            min_purity=0.5,
        )
        # Purity is 0.3, which is below 0.5, so all CS should be filtered out
        assert result.n_cs == 0

    def test_union_without_ld_matrices_no_purity(self):
        """Combine without ld_matrices produces purity=None."""
        snpids = ["1-1000-A-G", "1-1100-A-G"]
        cred1 = self._make_cred_with_snpids(
            snpids, {"1-1000-A-G": 0.8, "1-1100-A-G": 0.5}
        )
        cred2 = self._make_cred_with_snpids(
            snpids, {"1-1000-A-G": 0.7, "1-1100-A-G": 0.6}
        )
        result = combine_creds(
            [cred1, cred2],
            combine_cred="union",
            combine_pip="max",
        )
        assert result.purity is None


# ---------------------------------------------------------------------------
# TestFilterCredsetByPurityEdgeCases
# ---------------------------------------------------------------------------
class TestFilterCredsetByPurityEdgeCases:
    """Additional edge-case tests for filter_credset_by_purity."""

    def test_n_cs_zero_returns_original(self):
        """Return original when n_cs=0 and valid purity list."""
        from credtools.credibleset import filter_credset_by_purity

        cs = CredibleSet(
            tool="test",
            parameters={},
            coverage=0.95,
            n_cs=0,
            cs_sizes=[],
            lead_snps=[],
            snps=[],
            pips=pd.Series(dtype=float),
            purity=[],
        )
        result = filter_credset_by_purity(cs, min_purity=0.5)
        assert result.n_cs == 0
