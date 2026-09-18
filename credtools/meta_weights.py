"""Cohort-matched statistics and geometric ESS/SE LD aggregation.

The same eligible cohort/SNP contributions feed IVW statistics and LD. ESS is
an LD approximation; only SE weights represent the covariance of the IVW Z
scores under independent cohorts and correct within-cohort LD.
"""

import logging
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

from credtools.constants import ColName
from credtools.ldmatrix import LDMatrix
from credtools.locus import LocusSet

logger = logging.getLogger("META")


def validate_ld_weighting(ld_weighting: str) -> str:
    """Validate a public LD weighting choice before reading or writing data."""
    if ld_weighting not in ("ess", "se"):
        raise ValueError("ld_weighting must be 'ess' or 'se'")
    return str(getattr(ld_weighting, "value", ld_weighting))


def _canonical(
    frame: pd.DataFrame, first: str, second: str
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Copy a variant table and align alleles to its canonical SNP identifiers."""
    frame = frame.copy().reset_index(drop=True)
    if frame.SNPID.isna().any() or frame.SNPID.duplicated().any():
        raise ValueError("Missing or duplicate SNP identifiers in meta input")
    a = frame[first].astype(str).str.upper()
    b = frame[second].astype(str).str.upper()
    if not (a.str.fullmatch("[ACGT]+") & b.str.fullmatch("[ACGT]+") & a.ne(b)).all():
        raise ValueError("Invalid alleles in meta input")
    flip = (a > b).to_numpy()
    frame[first] = np.where(flip, b, a)
    frame[second] = np.where(flip, a, b)
    expected = (
        frame.CHR.astype(int).astype(str)
        + "-"
        + frame.BP.astype(int).astype(str)
        + "-"
        + frame[first]
        + "-"
        + frame[second]
    )
    if not frame.SNPID.equals(expected.rename("SNPID")):
        raise ValueError(
            "SNP identifiers must agree with chromosome, position and canonical alleles"
        )
    return frame, np.where(flip, -1.0, 1.0)


def matched_meta(inputs: LocusSet, ld_weighting: str = "ess"):
    """Return matched IVW statistics, geometric LD and contribution audits.

    For eligible cohort k at SNP j, w_kj is sample_size (ESS) or 1/SE_kj**2.
    a_kj = sqrt(w_kj / sum_k w_kj), and R = sum_k D(a_k) R_k D(a_k).
    Missing cohort/SNP contributions are zero; there is no pairwise
    renormalization. The output panel is the union of within-cohort matches.
    No input is mutated and no PSD repair is performed.
    """
    ld_weighting = validate_ld_weighting(ld_weighting)
    if not inputs.loci:
        raise ValueError("No input cohorts for meta-analysis")
    records = []
    cohort_rows = []
    for locus in inputs.loci:
        if not np.isfinite(locus.sample_size) or locus.sample_size <= 0:
            raise ValueError("Cohort sample_size/ESS must be finite and positive")
        ss, beta_sign = _canonical(locus.sumstats, "EA", "NEA")
        beta = pd.to_numeric(ss.BETA, errors="coerce").to_numpy(dtype=np.float64)
        se = pd.to_numeric(ss.SE, errors="coerce").to_numpy(dtype=np.float64)
        valid = np.isfinite(beta) & np.isfinite(se) & (se > 0)
        info = np.zeros(len(ss), dtype=np.float64)
        with np.errstate(
            over="ignore", divide="ignore", invalid="ignore", under="ignore"
        ):
            info[valid] = np.square(1.0 / se[valid])
        if np.any(valid & (~np.isfinite(info) | (info <= 0))):
            raise ValueError("SE information is outside the supported float64 range")
        ss["BETA"] = beta * beta_sign
        ss["SE"] = se
        ss["INFO"] = info
        eaf = pd.to_numeric(
            ss.get("EAF", pd.Series(np.nan, index=ss.index)), errors="coerce"
        )
        ss["EAF"] = np.where(beta_sign < 0, 1 - eaf, eaf)
        ss = ss.loc[valid].copy()
        ld = locus.ld
        if ld is None or ld.map.empty:
            ldmap = pd.DataFrame(columns=["SNPID", "CHR", "BP", "A1", "A2"])
            signs = np.empty(0)
        else:
            if ld.r.shape != (len(ld.map), len(ld.map)):
                raise ValueError("LD matrix must be square and match its variant map")
            ldmap, signs = _canonical(ld.map, "A1", "A2")
            if "AF2" in ldmap:
                af2 = pd.to_numeric(ldmap.AF2, errors="coerce")
                ldmap["AF2"] = np.where(signs < 0, 1 - af2, af2)
        eligible = ss.SNPID.isin(ldmap.SNPID)
        cohort_rows.append(
            {
                "popu": locus.popu,
                "cohort": locus.cohort,
                "sample_size": locus.sample_size,
                "n_sumstats": len(locus.sumstats),
                "n_valid_sumstats": len(ss),
                "n_ld_snps": len(ldmap),
                "n_matched": int(eligible.sum()),
            }
        )
        if not eligible.any():
            logger.warning(
                "No eligible GWAS/LD contribution from %s/%s", locus.popu, locus.cohort
            )
        records.append((locus, ss.set_index("SNPID", drop=False), ldmap, signs))

    template = pd.concat([r[1] for r in records], ignore_index=True)
    template = template.drop_duplicates("SNPID").sort_values(["CHR", "BP", "SNPID"])
    template = template.reset_index(drop=True)
    ids = pd.Index(template.SNPID)
    shape = (len(records), len(ids))
    information = np.zeros(shape, dtype=np.float64)
    beta = np.zeros(shape, dtype=np.float64)
    available = np.zeros(shape, dtype=bool)
    eaf = np.full(shape, np.nan, dtype=np.float64)
    af2 = np.full(shape, np.nan, dtype=np.float64)
    sizes = np.array([r[0].sample_size for r in records], dtype=np.float64)[:, None]
    for k, (_, ss, ldmap, _) in enumerate(records):
        aligned = ss.reindex(ids)
        information[k] = aligned.INFO.fillna(0).to_numpy(dtype=np.float64)
        beta[k] = aligned.BETA.fillna(0).to_numpy(dtype=np.float64)
        eaf[k] = aligned.EAF.to_numpy(dtype=np.float64)
        available[k] = ids.isin(ldmap.SNPID) & (information[k] > 0)
        if "AF2" in ldmap:
            af2[k] = (
                ldmap.set_index("SNPID").AF2.reindex(ids).to_numpy(dtype=np.float64)
            )
    matched_info = np.where(available, information, 0)
    full_total = information.sum(axis=0)
    total = matched_info.sum(axis=0)
    keep = total > 0
    if not keep.any():
        raise ValueError("No SNP has an eligible cohort GWAS/LD contribution")
    if not np.isfinite(full_total).all():
        raise ValueError("Total SE information overflow")

    out = template.loc[keep, ["SNPID", "CHR", "BP", "EA", "NEA"]].reset_index(drop=True)
    with np.errstate(over="ignore", invalid="ignore"):
        # Normalize first to avoid unnecessary overflow in information * beta.
        out["BETA"] = ((matched_info[:, keep] / total[keep]) * beta[:, keep]).sum(
            axis=0
        )
        out["SE"] = 1 / np.sqrt(total[keep])
        z = out.BETA.to_numpy() / out.SE.to_numpy()
    if not np.isfinite(z).all():
        raise ValueError("Nonfinite matched meta Z scores")
    out["P"] = np.maximum(2 * norm.sf(np.abs(z)), 1e-300)

    def average_frequency(values):
        good = available & np.isfinite(values) & (values >= 0) & (values <= 1)
        denom = (sizes * good).sum(axis=0)
        numerator = (sizes * np.where(good, values, 0)).sum(axis=0)
        result = np.full(len(ids), np.nan)
        np.divide(numerator, denom, out=result, where=denom > 0)
        return result[keep]

    out["EAF"] = average_frequency(eaf)
    out["MAF"] = np.minimum(out.EAF, 1 - out.EAF)
    ldmap_out = out[["SNPID", "CHR", "BP"]].copy()
    ldmap_out["A1"] = out.EA
    ldmap_out["A2"] = out.NEA
    if any("AF2" in r[2] for r in records):
        ldmap_out["AF2"] = average_frequency(af2)

    weights = matched_info if ld_weighting == "se" else sizes * available
    weights = weights[:, keep]
    weight_total = weights.sum(axis=0)
    if not np.isfinite(weight_total).all():
        raise ValueError("Total LD weight overflow")
    coefficients = np.sqrt(weights / weight_total).astype(np.float32)
    output_ids = pd.Index(out.SNPID)
    matrix = np.zeros((len(out), len(out)), dtype=np.float32)
    for k, (locus, _, ldmap, signs) in enumerate(records):
        selected = np.flatnonzero(weights[k] > 0)
        if not len(selected):
            continue
        source = pd.Index(ldmap.SNPID).get_indexer(output_ids[selected])
        a = coefficients[k, selected] * signs[source].astype(np.float32)
        # Bound temporary matrix allocations on large loci; input matrices are
        # never copied in full merely to match rows or flip alleles.
        for start in range(0, len(selected), 256):
            rows = slice(start, start + 256)
            raw = locus.ld.r[np.ix_(source[rows], source)].astype(np.float32)
            reverse = locus.ld.r[np.ix_(source, source[rows])].T
            if not np.isfinite(raw).all() or not np.allclose(
                raw, reverse, rtol=0, atol=2e-5
            ):
                raise ValueError(
                    "Matched LD must be finite and symmetric (no implicit repair)"
                )
            if np.any(np.abs(raw) > 1.00002):
                raise ValueError("LD correlations must be within [-1, 1]")
            raw *= a[rows, None]
            raw *= a[None, :]
            matrix[np.ix_(selected[rows], selected)] += raw
        if not np.allclose(locus.ld.r[source, source], 1, rtol=0, atol=2e-5):
            raise ValueError("Matched LD must have a unit diagonal")
    # Only remove accumulation roundoff, not negative eigenvalues.
    for start in range(0, len(out), 256):
        stop = min(start + 256, len(out))
        block = (matrix[start:stop] + matrix[:, start:stop].T) * np.float32(0.5)
        matrix[start:stop] = block
        matrix[:, start:stop] = block.T
    np.fill_diagonal(matrix, 1)

    audit = template[["SNPID"]].copy()
    audit["n_valid_cohorts"] = (information > 0).sum(axis=0)
    audit["n_matched_cohorts"] = available.sum(axis=0)
    audit["full_information"] = full_total
    audit["matched_information"] = total
    audit["information_fraction"] = total / full_total
    audit["retained"] = keep
    audit["ld_weighting"] = ld_weighting
    out = out.reindex(columns=ColName.sumstat_cols)
    return out, LDMatrix(ldmap_out, matrix), audit, pd.DataFrame(cohort_rows)
