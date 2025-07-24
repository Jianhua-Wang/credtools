"""Define constants used in the package."""

import numpy as np


class ColName:
    """Define column names."""

    # mandatory columns
    CHR = "CHR"
    BP = "BP"
    RSID = "rsID"
    EA = "EA"
    NEA = "NEA"
    P = "P"
    BETA = "BETA"
    SE = "SE"
    EAF = "EAF"

    # ld map columns
    A1 = "A1"
    A2 = "A2"

    # optional columns
    MAF = "MAF"
    N = "N"
    Z = "Z"
    INFO = "INFO"

    # columns for loci
    START = "START"
    END = "END"
    LEAD_SNP = "LEAD_SNP"
    LEAD_SNP_P = "LEAD_SNP_P"
    LEAD_SNP_BP = "LEAD_SNP_BP"

    # unique snpid, chr-bp-sorted(EA,NEA)
    SNPID = "SNPID"

    # COJO results
    COJO_P = "COJO_P"
    COJO_BETA = "COJO_BETA"
    COJO_SE = "COJO_SE"

    # posterior probability
    PIP = "PIP"
    FINEMAP = "FINEMAP"
    ABF = "ABF"
    SUSIE = "SUSIE"
    # PP_POLYFUN_FINEMAP = "PP_POLYFUN_FINEMAP"
    # PP_POLYFUN_SUSIE = "PP_POLYFUN_SUSIE"

    # ordered columns
    mandatory_cols = [CHR, BP, EA, NEA, EAF, BETA, SE, P]
    sumstat_cols = [SNPID, CHR, BP, RSID, EA, NEA, EAF, MAF, BETA, SE, P]
    loci_cols = [CHR, START, END, LEAD_SNP, LEAD_SNP_P, LEAD_SNP_BP]
    map_cols = [CHR, BP, A1, A2]


class Method:
    """Define methods."""

    FINEMAP = "FINEMAP"
    SUSIE = "SUSIE"
    ABF = "ABF"
    RSparsePro = "RSparsePro"
    SUSIEX = "SUSIEX"
    MULTISUSIE = "MULTISUSIE"


# only support autosomes
CHROM_LIST = [i for i in range(1, 24)]

chrom_len = {
    1: 249250621,
    2: 243199373,
    3: 198022430,
    4: 191154276,
    5: 180915260,
    6: 171115067,
    7: 159138663,
    8: 146364022,
    9: 141213431,
    10: 135534747,
    11: 135006516,
    12: 133851895,
    13: 115169878,
    14: 107349540,
    15: 102531392,
    16: 90354753,
    17: 81195210,
    18: 78077248,
    19: 59128983,
    20: 63025520,
    21: 48129895,
    22: 51304566,
}


class ColType:
    """Define column types."""

    CHR = int
    BP = int
    RSID = str
    EA = str
    NEA = str
    P = float
    NEGLOGP = float
    BETA = float
    OR = float
    ORSE = float
    SE = float
    EAF = float
    MAF = float
    N = int
    Z = float
    INFO = float


class ColAllowNA:
    """Define whether a column allows missing values."""

    CHR = False
    BP = False
    RSID = True
    EA = False
    NEA = False
    P = True
    BETA = True
    OR = True
    ORSE = True
    SE = True
    EAF = True
    MAF = True
    N = True
    Z = True
    INFO = True


class ColRange:
    """Define the range of values for each column."""

    INFO = (0, 1)
    CHR_MIN = 1
    CHR_MAX = 23
    BP_MIN = 0
    BP_MAX = 300000000
    P_MIN = 0
    P_MAX = 1
    SE_MIN = 0
    SE_MAX = np.inf
    EAF_MIN = 0
    EAF_MAX = 1
    MAF_MIN = 0
    MAF_MAX = 1
    ORSE_MIN = 0
    ORSE_MAX = np.inf
    NEGLOGP_MIN = 0
