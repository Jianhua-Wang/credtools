"""Wrapper for MultiSuSiE multi-ancestry fine-mapping."""

import json
import logging
from typing import List, Optional

import numpy as np
import pandas as pd

from credtools.constants import ColName, Method
from credtools.credibleset import CredibleSet
from credtools.locus import LocusSet, intersect_sumstat_ld
from credtools.wrappers.multisusie_rss import multisusie_rss

logger = logging.getLogger("MULTISUSIE")


def run_multisusie(
    locus_set: LocusSet,
    max_causal: int = 1,
    coverage: float = 0.95,
    rho: float = 0.75,
    scaled_prior_variance: float = 0.2,
    standardize: bool = False,
    pop_spec_standardization: bool = True,
    estimate_residual_variance: bool = True,
    estimate_prior_variance: bool = True,
    estimate_prior_method: str = "early_EM",
    pop_spec_effect_priors: bool = True,
    iter_before_zeroing_effects: int = 5,
    prior_tol: float = 1e-9,
    max_iter: int = 100,
    tol: float = 1e-3,
    min_abs_corr: float = 0.1,
) -> CredibleSet:
    """
    Run MultiSuSiE for multi-ancestry fine-mapping analysis.

    MultiSuSiE extends the SuSiE framework to jointly analyze summary statistics
    from multiple populations/ancestries, allowing for effect size correlation
    across populations while accounting for population-specific LD structures.
    This enables more powerful fine-mapping by leveraging shared causal variants
    across diverse populations.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object containing multiple Locus objects, each representing
        summary statistics and LD data from different populations/ancestries.
        All loci should cover the same genomic region with overlapping variants.
    max_causal : int, optional
        Maximum number of causal variants (L parameter), by default 1.
        This determines the number of single-effect components in the model
        across all populations.
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95.
        This determines the cumulative posterior probability mass
        included in each credible set.
    rho : float, optional
        Prior correlation between causal effect sizes across populations, by default 0.75.
        Higher values assume more similar effect sizes between populations,
        while lower values allow for more population-specific effects.
        Can also be provided as a K×K matrix for population-specific correlations.
    scaled_prior_variance : float, optional
        Scaled prior variance for effect sizes, by default 0.2.
        This parameter controls the expected magnitude of causal effects,
        scaled by the residual variance. Larger values allow for larger effects.
    standardize : bool, optional
        Whether to standardize genotypes to have variance 1, by default False.
        Standardization can improve numerical stability but may affect
        interpretation of effect sizes.
    pop_spec_standardization : bool, optional
        Whether to perform population-specific standardization, by default True.
        If True and standardize=True, standardizes within each population separately.
        If False, uses pooled standardization across all populations.
    estimate_residual_variance : bool, optional
        Whether to estimate residual variance from data, by default True.
        If False, residual variance is assumed to be 1 for all populations.
    estimate_prior_variance : bool, optional
        Whether to estimate prior variance adaptively, by default True.
        Adaptive estimation can improve model fit by optimizing the
        prior variance based on the data.
    estimate_prior_method : str, optional
        Method for prior variance estimation, by default "early_EM".
        Options include "early_EM", "EM", "optim", or None.
        "early_EM" provides good balance between speed and accuracy.
    pop_spec_effect_priors : bool, optional
        Whether to use population-specific effect size priors, by default True.
        Allows for different prior variances across populations,
        accommodating population-specific effect size distributions.
    iter_before_zeroing_effects : int, optional
        Number of iterations before zeroing weak effects, by default 5.
        Components with low likelihood are removed after this many iterations
        to improve computational efficiency and model parsimony.
    prior_tol : float, optional
        Tolerance for minimum prior variance, by default 1e-9.
        Components with prior variance below this threshold are excluded
        from credible set calculations.
    max_iter : int, optional
        Maximum number of iterations for the IBSS algorithm, by default 100.
        More iterations may improve convergence but increase runtime.
    tol : float, optional
        Convergence tolerance for the ELBO, by default 1e-3.
        Algorithm stops when ELBO change falls below this threshold.
    min_abs_corr : float, optional
        Minimum absolute correlation for credible set purity, by default 0.1.
        Credible sets with pairwise correlations below this threshold
        may be filtered based on purity criteria.

    Returns
    -------
    CredibleSet
        Credible set object containing:
        - Posterior inclusion probabilities aggregated across populations
        - Credible sets for each detected signal
        - Lead SNPs (highest PIP in each credible set)
        - Coverage probabilities and purity measures
        - Multi-ancestry algorithm parameters

    Notes
    -----
    MultiSuSiE implements a multi-population extension of the SuSiE model:

    y_k = Σ(l=1 to L) X_k * b_k,l + ε_k

    where:
    - y_k is the phenotype vector for population k
    - X_k is the genotype matrix for population k
    - b_k,l is the l-th single-effect vector for population k
    - ε_k is the residual error for population k
    - Effect sizes b_k,l are correlated across populations via correlation matrix ρ

    Key innovations:

    1. **Cross-population effect correlation**: Models correlation structure
       of causal effects across populations using correlation matrix ρ

    2. **Population-specific LD**: Accounts for different LD patterns
       across populations while sharing causal variant locations

    3. **Adaptive prior estimation**: Estimates population-specific or
       shared prior variances based on the data

    4. **Joint credible sets**: Constructs credible sets that aggregate
       evidence across all populations

    The algorithm workflow:
    1. Convert GWAS summary statistics to sufficient statistics for each population
    2. Initialize L single-effect regressions with correlated priors
    3. Iteratively update effects using multi-population variational Bayes
    4. Estimate residual and prior variances adaptively
    5. Monitor convergence using multi-population ELBO
    6. Construct joint credible sets across populations

    Advantages over single-population methods:
    - Increased power through multi-ancestry meta-analysis
    - Better resolution through diverse LD patterns
    - Robustness to population-specific confounding
    - Natural handling of trans-ethnic fine-mapping

    The method performs minor allele frequency (MAF) filtering automatically:
    - Variants with MAF < single_population_mac_thresh in any population are censored
    - Variants with MAF < multi_population_maf_thresh in ALL populations are removed
    - This prevents spurious associations from population-specific rare variants

    Reference:
    Zou, Y. et al. Fine-mapping from summary data with the "Sum of Single Effects" model.
    PLoS Genet. 18, e1010299 (2022).

    Examples
    --------
    >>> # Basic multi-ancestry fine-mapping
    >>> credible_set = run_multisusie(locus_set)
    >>> print(f"Found {credible_set.n_cs} credible sets")
    >>> print(f"Populations analyzed: {len(locus_set.loci)}")
    Found 1 credible sets
    Populations analyzed: 3

    >>> # Multi-ancestry with multiple causal variants and custom correlation
    >>> credible_set = run_multisusie(
    ...     locus_set,
    ...     max_causal=5,
    ...     rho=0.8,  # Higher cross-population correlation
    ...     coverage=0.99,
    ...     pop_spec_effect_priors=True
    ... )
    >>> print(f"Detected {credible_set.n_cs} independent signals")
    >>> print(f"Credible set sizes: {credible_set.cs_sizes}")
    Detected 2 independent signals
    Credible set sizes: [8, 15]

    >>> # Access aggregated posterior inclusion probabilities
    >>> top_variants = credible_set.pips.nlargest(10)
    >>> print("Top 10 variants by aggregated PIP:")
    >>> print(top_variants)
    Top 10 variants by aggregated PIP:
    rs123456    0.9234
    rs789012    0.8456
    rs345678    0.7789
    ...

    >>> # Examine cross-population evidence
    >>> for i, snps in enumerate(credible_set.snps):
    ...     lead_snp = credible_set.lead_snps[i]
    ...     pip = credible_set.pips[lead_snp]
    ...     print(f"Signal {i+1}: {lead_snp} (PIP={pip:.4f})")
    Signal 1: rs123456 (PIP=0.9234)
    Signal 2: rs789012 (PIP=0.8456)
    """
    logger.info(f"Running MultiSuSiE on {locus_set}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "rho": rho,
        "scaled_prior_variance": scaled_prior_variance,
        "standardize": standardize,
        "pop_spec_standardization": pop_spec_standardization,
        "estimate_residual_variance": estimate_residual_variance,
        "estimate_prior_variance": estimate_prior_variance,
        "estimate_prior_method": estimate_prior_method,
        "pop_spec_effect_priors": pop_spec_effect_priors,
        "iter_before_zeroing_effects": iter_before_zeroing_effects,
        "prior_tol": prior_tol,
        "max_iter": max_iter,
        "tol": tol,
        "min_abs_corr": min_abs_corr,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")

    all_variants = []
    for locus in locus_set.loci:
        locus = intersect_sumstat_ld(locus)
        all_variants.append(
            locus.sumstats[
                [ColName.SNPID, ColName.CHR, ColName.BP, ColName.EA, ColName.NEA]
            ]
        )
    all_variants = pd.concat(all_variants, axis=0)
    all_variants.drop_duplicates(subset=[ColName.SNPID], inplace=True)
    all_variants.sort_values([ColName.CHR, ColName.BP], inplace=True)
    variant_to_index = {variant: i for i, variant in enumerate(all_variants["SNPID"])}
    # TODO: make concat of all_variants as a function of locus_set
    # TODO: make intersect of loci as a function of all_variants
    # TODO: add a switch of either using all_variants or the joint_variants
    # joint_variants = []
    # for i, locus in enumerate(locus_set.loci):
    #     locus = intersect_sumstat_ld(locus)
    #     sumstat = locus.sumstats.copy()
    #     if i == 0:
    #         joint_variants = set(sumstat["SNPID"].values)
    #     else:
    #         joint_variants = joint_variants.intersection(set(sumstat["SNPID"].values))
    # joint_variants = sumstat[sumstat["SNPID"].isin(joint_variants)][["SNPID", "CHR", "BP", "EA", "NEA"]].copy()
    # joint_variants.sort_values(["CHR", "BP"], inplace=True)
    # all_variants = joint_variants.copy()
    # variant_to_index = {variant: i for i, variant in enumerate(all_variants["SNPID"])}

    z_list = []
    R_list = []
    for locus in locus_set.loci:
        sumstat = locus.sumstats.copy()
        ldmap = locus.ld.map
        ld = locus.ld.r
        sumstat[ColName.Z] = sumstat[ColName.BETA] / sumstat[ColName.SE]
        sumstat.set_index(ColName.SNPID, inplace=True)
        z = all_variants[ColName.SNPID].map(sumstat[ColName.Z]).values
        z_list.append(z)
        expand_ld = np.zeros((all_variants.shape[0], all_variants.shape[0]))
        intersec_index = ldmap[ldmap["SNPID"].isin(all_variants["SNPID"])].index
        ldmap = ldmap.loc[intersec_index]
        ld = ld[intersec_index, :][:, intersec_index]
        study_indices = np.array([variant_to_index[snp] for snp in ldmap["SNPID"]])
        idx_i, idx_j = np.meshgrid(study_indices, study_indices)
        expand_ld[idx_i, idx_j] += ld.astype(np.float32)
        np.fill_diagonal(expand_ld, 1)
        R_list.append(expand_ld)

    rho_array = np.full((len(locus_set.loci), len(locus_set.loci)), rho)
    np.fill_diagonal(rho_array, 1)
    ss_fit = multisusie_rss(
        z_list=z_list,
        R_list=R_list,
        population_sizes=[locus.sample_size for locus in locus_set.loci],
        rho=rho_array,  # type: ignore
        L=max_causal,
        coverage=coverage,
        scaled_prior_variance=scaled_prior_variance,
        max_iter=max_iter,
        tol=tol,
        pop_spec_standardization=pop_spec_standardization,
        estimate_residual_variance=estimate_residual_variance,
        estimate_prior_variance=estimate_prior_variance,
        estimate_prior_method=estimate_prior_method,
        pop_spec_effect_priors=pop_spec_effect_priors,
        iter_before_zeroing_effects=iter_before_zeroing_effects,
        prior_tol=prior_tol,
        min_abs_corr=min_abs_corr,
        float_type=np.float32,
        low_memory_mode=False,
        recover_R=False,
        single_population_mac_thresh=20,
        mac_list=None,
        multi_population_maf_thresh=0,
        maf_list=None,
    )
    pip = pd.Series(index=all_variants[ColName.SNPID].tolist(), data=ss_fit.pip)
    cs_snp = []
    for cs_snp_idx in ss_fit.sets[0]:
        if len(cs_snp_idx) > 0 and len(cs_snp_idx) < len(pip):
            snps = all_variants[ColName.SNPID].to_numpy()[cs_snp_idx]
            cs_snp.append(snps.tolist())
    cs_sizes = [len(snpids) for snpids in cs_snp]
    lead_snps = [str(pip[pip.index.isin(snpids)].idxmax()) for snpids in cs_snp]

    logger.info(f"Finished MultiSuSiE on {locus_set}")
    logger.info(f"N of credible set: {len(cs_snp)}")
    logger.info(f"Credible set size: {cs_sizes}")

    return CredibleSet(
        tool=Method.MULTISUSIE,
        n_cs=len(cs_snp),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snp,
        cs_sizes=cs_sizes,
        pips=pip,
        parameters=parameters,
    )
