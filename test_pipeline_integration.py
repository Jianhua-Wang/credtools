#!/usr/bin/env python3
"""Test script for ABF+COJO integration with main pipeline."""

import numpy as np
import pandas as pd
from credtools.locus import Locus, LocusSet
from credtools.ldmatrix import LDMatrix
from credtools.credtools import fine_map
from credtools.constants import ColName


def create_single_locus_set():
    """Create a LocusSet with one locus for testing."""
    # Create mock summary statistics
    n_snps = 50
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
        ColName.P: np.random.uniform(1e-6, 0.1, n_snps)
    }
    
    # Add one strong signal
    sumstats_data[ColName.P][10] = 1e-10
    sumstats_data[ColName.BETA][10] = 0.5
    sumstats_data[ColName.SE][10] = 0.02
    
    sumstats = pd.DataFrame(sumstats_data)
    
    # Create mock LD matrix
    ld_matrix = np.eye(n_snps)
    for i in range(n_snps):
        for j in range(max(0, i-3), min(n_snps, i+4)):
            if i != j:
                ld_matrix[i, j] = 0.05 * np.exp(-abs(i-j)/1)
    
    # Create LD map
    ld_map_data = {
        ColName.SNPID: [f"rs{i}" for i in range(n_snps)],
        ColName.CHR: [1] * n_snps,
        ColName.BP: list(range(1000000, 1000000 + n_snps * 1000, 1000)),
        ColName.A1: ["A"] * n_snps,
        ColName.A2: ["G"] * n_snps,
        "AF2": np.random.uniform(0.1, 0.9, n_snps)
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
    
    # Create LocusSet
    locus_set = LocusSet([locus])
    
    return locus_set


def test_pipeline_integration():
    """Test ABF+COJO integration with the main fine_map function."""
    print("Creating single locus set...")
    locus_set = create_single_locus_set()
    
    print("Testing ABF+COJO integration with main pipeline...")
    
    # Test 1: ABF+COJO with single_input strategy
    try:
        print("\n1. Testing abf_cojo with single_input strategy:")
        result = fine_map(
            locus_set=locus_set,
            strategy="single_input",
            tool="abf_cojo",
            max_causal=5,
            coverage=0.95,
            var_prior=0.2,
            p_cutoff=1e-6
        )
        
        print(f"   ✅ Success!")
        print(f"   Tool: {result.tool}")
        print(f"   Number of credible sets: {result.n_cs}")
        print(f"   Credible set sizes: {result.cs_sizes}")
        print(f"   Top 3 PIPs: {result.pips.nlargest(3)}")
        
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test 2: Compare with standard ABF
    try:
        print("\n2. Testing standard abf for comparison:")
        abf_result = fine_map(
            locus_set=locus_set,
            strategy="single_input",
            tool="abf",
            max_causal=1,
            coverage=0.95,
            var_prior=0.2
        )
        
        print(f"   ✅ Success!")
        print(f"   Tool: {abf_result.tool}")
        print(f"   Number of credible sets: {abf_result.n_cs}")
        print(f"   Credible set sizes: {abf_result.cs_sizes}")
        print(f"   Top 3 PIPs: {abf_result.pips.nlargest(3)}")
        
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test 3: ABF+COJO with post_hoc_combine strategy (if multiple loci)
    try:
        print("\n3. Testing abf_cojo with post_hoc_combine strategy:")
        # Create a second locus by duplicating the first one with slight modifications
        locus2 = Locus(
            popu="EUR",
            cohort="test2", 
            sample_size=8000,
            sumstats=locus_set.loci[0].sumstats.copy(),
            ld=locus_set.loci[0].ld
        )
        multi_locus_set = LocusSet([locus_set.loci[0], locus2])
        
        multi_result = fine_map(
            locus_set=multi_locus_set,
            strategy="post_hoc_combine",
            tool="abf_cojo",
            max_causal=5,
            coverage=0.95,
            var_prior=0.2,
            p_cutoff=1e-6,
            combine_cred="union"
        )
        
        print(f"   ✅ Success!")
        print(f"   Tool: {multi_result.tool}")
        print(f"   Number of credible sets: {multi_result.n_cs}")
        print(f"   Credible set sizes: {multi_result.cs_sizes}")
        print(f"   Top 3 PIPs: {multi_result.pips.nlargest(3)}")
        
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n🎉 All pipeline integration tests passed!")
    return True


if __name__ == "__main__":
    success = test_pipeline_integration()
    if success:
        print("✅ ABF+COJO pipeline integration test passed!")
    else:
        print("❌ ABF+COJO pipeline integration test failed!")