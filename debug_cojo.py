#!/usr/bin/env python3
"""Debug script to understand COJO conditional analysis output."""

import numpy as np
import pandas as pd
import tempfile
from cojopy.cojopy import COJO
from credtools.constants import ColName


def test_cojo_conditional():
    """Test COJO conditional analysis to understand output format."""
    # Create simple test data
    n_snps = 20
    np.random.seed(42)
    
    # Create summary statistics in COJO format
    sumstats = pd.DataFrame({
        "SNP": [f"rs{i}" for i in range(n_snps)],
        "A1": ["A"] * n_snps,
        "A2": ["G"] * n_snps,
        "b": np.random.normal(0, 0.1, n_snps),
        "se": np.random.uniform(0.02, 0.05, n_snps),
        "p": np.random.uniform(1e-6, 0.01, n_snps),
        "freq": np.random.uniform(0.1, 0.9, n_snps),
        "N": [10000] * n_snps
    })
    
    # Add strong signals
    sumstats.loc[5, "p"] = 1e-8
    sumstats.loc[5, "b"] = 0.3
    sumstats.loc[15, "p"] = 5e-8  
    sumstats.loc[15, "b"] = 0.25
    
    # Create LD matrix
    ld_matrix = np.eye(n_snps)
    for i in range(n_snps):
        for j in range(max(0, i-2), min(n_snps, i+3)):
            if i != j:
                ld_matrix[i, j] = 0.05
    
    # Create COJO object
    c = COJO(p_cutoff=1e-6, collinear_cutoff=0.9)
    c.load_sumstats(sumstats=sumstats, ld_matrix=ld_matrix, ld_freq=None)
    
    # Run conditional selection first
    print("Running conditional selection...")
    selection_results = c.conditional_selection()
    print(f"Selection results columns: {selection_results.columns.tolist()}")
    print(f"Selection results shape: {selection_results.shape}")
    print(f"Selected SNPs: {selection_results['SNP'].tolist()}")
    print(selection_results.head())
    
    if len(selection_results) > 1:
        # Try conditional analysis
        print("\nRunning conditional analysis...")
        independent_snps = selection_results["SNP"].tolist()
        conditioning_snps = independent_snps[1:]  # Condition on all but first
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            for snp in conditioning_snps:
                f.write(f"{snp}\n")
            cond_snps_file = f.name
        
        try:
            conditional_results = c.run_conditional_analysis(cond_snps_path=cond_snps_file)
            print(f"Conditional results columns: {conditional_results.columns.tolist()}")
            print(f"Conditional results shape: {conditional_results.shape}")
            print(conditional_results.head())
        except Exception as e:
            print(f"Error in conditional analysis: {e}")
        finally:
            import os
            os.unlink(cond_snps_file)
    else:
        print("Only one signal found, cannot test conditional analysis")


if __name__ == "__main__":
    test_cojo_conditional()