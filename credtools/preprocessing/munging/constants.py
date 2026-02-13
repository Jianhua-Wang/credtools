"""
Constants and column definitions for munging GWAS summary statistics.

Adapted from smunger (https://github.com/Jianhua-Wang/smunger)
Original author: Jianhua Wang
License: MIT

This module re-exports shared definitions from the main credtools.constants
and provides munging-specific constants (output columns, column name aliases,
chromosome lengths, and column mapping suggestions).
"""

from credtools.constants import ColAllowNA, ColName, ColRange, ColType  # noqa: F401

# Munging-specific: output columns for credtools munge command
OUTPUT_COLS = [
    ColName.CHR,
    ColName.BP,
    ColName.SNPID,
    ColName.EA,
    ColName.NEA,
    ColName.EAF,
    ColName.BETA,
    ColName.SE,
    ColName.P,
    ColName.N,
    ColName.RSID,
]


# Common alternative column names found in GWAS files
COMMON_COLNAMES = {
    # Chromosome
    "CHROM": ColName.CHR,
    "CHR": ColName.CHR,
    "#CHROM": ColName.CHR,
    "chromosome": ColName.CHR,
    "Chromosome": ColName.CHR,
    # Base pair position
    "BP": ColName.BP,
    "POS": ColName.BP,
    "Position": ColName.BP,
    "position": ColName.BP,
    "base_pair_location": ColName.BP,
    "pos": ColName.BP,
    # SNP ID
    "SNP": ColName.SNPID,
    "SNPID": ColName.SNPID,
    "MarkerName": ColName.SNPID,
    "variant": ColName.SNPID,
    "ID": ColName.SNPID,
    # Effect allele
    "A1": ColName.EA,
    "EA": ColName.EA,
    "effect_allele": ColName.EA,
    "ALT": ColName.EA,
    "Allele1": ColName.EA,
    # Non-effect allele
    "A2": ColName.NEA,
    "NEA": ColName.NEA,
    "other_allele": ColName.NEA,
    "REF": ColName.NEA,
    "Allele2": ColName.NEA,
    # Effect allele frequency
    "EAF": ColName.EAF,
    "FRQ": ColName.EAF,
    "FREQ": ColName.EAF,
    "MAF": ColName.MAF,
    "frequency": ColName.EAF,
    "Freq1": ColName.EAF,
    # Effect size (beta)
    "BETA": ColName.BETA,
    "beta": ColName.BETA,
    "Beta": ColName.BETA,
    "effect": ColName.BETA,
    "Effect": ColName.BETA,
    # Standard error
    "SE": ColName.SE,
    "StdErr": ColName.SE,
    "stderr": ColName.SE,
    "standard_error": ColName.SE,
    # P-value
    "P": ColName.P,
    "PVAL": ColName.P,
    "P_BOLT_LMM": ColName.P,
    "pvalue": ColName.P,
    "P-value": ColName.P,
    "p_value": ColName.P,
    # Sample size
    "N": ColName.N,
    "n": ColName.N,
    "sample_size": ColName.N,
    "NMISS": ColName.N,
    # Info score
    "INFO": ColName.INFO,
    "info": ColName.INFO,
    "imputation_quality": ColName.INFO,
    # Z-score
    "Z": ColName.Z,
    "STAT": ColName.Z,
    "zscore": ColName.Z,
    "z_score": ColName.Z,
    # rsID
    "RSID": ColName.RSID,
    "rsid": ColName.RSID,
    "rs": ColName.RSID,
}


# Chromosome lengths (GRCh37/hg19)
CHROM_LENGTHS = {
    1: 249250621,
    2: 242193529,
    3: 198295559,
    4: 190214555,
    5: 181538259,
    6: 170805979,
    7: 159345973,
    8: 145138636,
    9: 138394717,
    10: 133797422,
    11: 135086622,
    12: 133275309,
    13: 114364328,
    14: 107043718,
    15: 101991189,
    16: 90338345,
    17: 83257441,
    18: 80373285,
    19: 58617616,
    20: 64444167,
    21: 46709983,
    22: 50818468,
    23: 156040895,  # X chromosome
}


def suggest_column_mapping(header_list):
    """
    Suggest column mappings for interactive configuration.

    Parameters
    ----------
    header_list : list
        List of column headers from input file.

    Returns
    -------
    dict
        Suggested mappings with confidence scores.
    """
    suggestions = {}

    for header in header_list:
        if header in COMMON_COLNAMES:
            suggestions[header] = {
                "suggested": COMMON_COLNAMES[header],
                "confidence": "high",
            }
        else:
            # Try fuzzy matching for common patterns
            header_lower = header.lower()

            if any(x in header_lower for x in ["chr", "chrom"]):
                suggestions[header] = {"suggested": ColName.CHR, "confidence": "medium"}
            elif any(x in header_lower for x in ["pos", "bp"]):
                suggestions[header] = {"suggested": ColName.BP, "confidence": "medium"}
            elif any(x in header_lower for x in ["beta", "effect"]):
                suggestions[header] = {
                    "suggested": ColName.BETA,
                    "confidence": "medium",
                }
            elif any(x in header_lower for x in ["pval", "p_val", "pvalue"]):
                suggestions[header] = {"suggested": ColName.P, "confidence": "medium"}
            else:
                suggestions[header] = {"suggested": header, "confidence": "low"}

    return suggestions
