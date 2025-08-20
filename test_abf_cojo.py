#!/usr/bin/env python3
"""Simple test script for ABF+COJO implementation."""

import numpy as np
import pandas as pd
from credtools.locus import Locus
from credtools.ldmatrix import LDMatrix
from credtools.wrappers.abf_cojo import run_abf_cojo
from credtools.wrappers.abf import run_abf
from credtools.constants import ColName


def create_mock_data():
    """Create mock data for testing ABF+COJO."""
    # Create mock summary statistics
    n_snps = 100
    np.random.seed(42)
    
    sumstats_data = {
        ColName.SNPID: [f"rs{i}" for i in range(n_snps)],
        ColName.CHR: [1] * n_snps,
        ColName.BP: list(range(1000000, 1000000 + n_snps * 1000, 1000)),
        ColName.RSID: [f"rs{i}" for i in range(n_snps)],
        ColName.EA: ["A"] * n_snps,
        ColName.NEA: ["G"] * n_snps,
        ColName.EAF: np.random.uniform(0.1, 0.9, n_snps),
        ColName.BETA: np.random.normal(0, 0.1, n_snps),
        ColName.SE: np.random.uniform(0.01, 0.05, n_snps),
        ColName.P: np.random.uniform(1e-8, 0.1, n_snps)
    }
    
    # Add some strong signals
    sumstats_data[ColName.P][10] = 1e-10  # Strong signal 1
    sumstats_data[ColName.BETA][10] = 0.5
    sumstats_data[ColName.SE][10] = 0.02
    
    sumstats_data[ColName.P][30] = 5e-9   # Strong signal 2 (farther away)
    sumstats_data[ColName.BETA][30] = 0.3
    sumstats_data[ColName.SE][30] = 0.02
    
    sumstats_data[ColName.P][70] = 1e-8   # Strong signal 3 (even farther)
    sumstats_data[ColName.BETA][70] = 0.25
    sumstats_data[ColName.SE][70] = 0.02
    
    sumstats = pd.DataFrame(sumstats_data)
    
    # Create mock LD matrix (diagonal with some local LD structure)
    ld_matrix = np.eye(n_snps)
    for i in range(n_snps):
        for j in range(max(0, i-3), min(n_snps, i+4)):
            if i != j:
                # Small local LD, but signals at positions 10, 30, 70 are independent
                ld_matrix[i, j] = 0.05 * np.exp(-abs(i-j)/1)
    
    # Create LD map
    ld_map_data = {
        ColName.SNPID: [f"rs{i}" for i in range(n_snps)],
        ColName.CHR: [1] * n_snps,
        ColName.BP: list(range(1000000, 1000000 + n_snps * 1000, 1000)),
        ColName.A1: ["A"] * n_snps,
        ColName.A2: ["G"] * n_snps,
        "AF2": np.random.uniform(0.1, 0.9, n_snps)  # Reference allele frequency
    }
    ld_map = pd.DataFrame(ld_map_data)
    
    ld = LDMatrix(map_df=ld_map, r=ld_matrix)
    
    # Create locus
    locus = Locus(
        popu="EUR",
        cohort="test",
        sample_size=10000,
        sumstats=sumstats,
        ld=ld
    )
    
    return locus


def test_abf_cojo():
    """Test ABF+COJO implementation."""
    print("Creating mock data...")
    locus = create_mock_data()
    
    print(f"Testing ABF+COJO on locus with {len(locus.sumstats)} SNPs")
    print(f"Sample size: {locus.sample_size}")
    
    # Test ABF+COJO
    print("\nRunning ABF+COJO...")
    try:
        abf_cojo_result = run_abf_cojo(
            locus,
            max_causal=5,
            coverage=0.95,
            var_prior=0.2,
            p_cutoff=1e-6  # More lenient for mock data
        )
        
        print(f"ABF+COJO Results:")
        print(f"  Number of credible sets: {abf_cojo_result.n_cs}")
        print(f"  Credible set sizes: {abf_cojo_result.cs_sizes}")
        print(f"  Lead SNPs: {abf_cojo_result.lead_snps}")
        print(f"  Top 5 PIPs: {abf_cojo_result.pips.nlargest(5)}")
        
    except Exception as e:
        print(f"Error in ABF+COJO: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Compare with standard ABF for single signal case
    print("\nRunning standard ABF for comparison...")
    try:
        abf_result = run_abf(
            locus,
            max_causal=1,
            coverage=0.95,
            var_prior=0.2
        )
        
        print(f"Standard ABF Results:")
        print(f"  Number of credible sets: {abf_result.n_cs}")
        print(f"  Credible set sizes: {abf_result.cs_sizes}")
        print(f"  Lead SNPs: {abf_result.lead_snps}")
        print(f"  Top 5 PIPs: {abf_result.pips.nlargest(5)}")
        
    except Exception as e:
        print(f"Error in standard ABF: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\nTest completed successfully!")
    return True


if __name__ == "__main__":
    success = test_abf_cojo()
    if success:
        print("✅ ABF+COJO implementation test passed!")
    else:
        print("❌ ABF+COJO implementation test failed!")