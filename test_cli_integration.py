#!/usr/bin/env python3
"""Test script for ABF+COJO CLI integration."""

import os
import tempfile
import numpy as np
import pandas as pd
from credtools.cli import app
from credtools.locus import Locus
from credtools.ldmatrix import LDMatrix
from credtools.constants import ColName


def create_test_loci_file():
    """Create a test loci file for CLI testing."""
    # Create a temporary directory
    temp_dir = tempfile.mkdtemp()
    
    # Create mock data
    n_snps = 30
    np.random.seed(42)
    
    # Create summary statistics
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
        ColName.P: np.random.uniform(1e-6, 0.01, n_snps)
    }
    
    # Add strong signal
    sumstats_data[ColName.P][10] = 1e-9
    sumstats_data[ColName.BETA][10] = 0.4
    sumstats_data[ColName.SE][10] = 0.02
    
    sumstats = pd.DataFrame(sumstats_data)
    
    # Create LD matrix 
    ld_matrix = np.eye(n_snps)
    for i in range(n_snps):
        for j in range(max(0, i-2), min(n_snps, i+3)):
            if i != j:
                ld_matrix[i, j] = 0.05
    
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
    
    # Save files
    sumstats_file = f"{temp_dir}/test_locus.sumstats"
    ld_file = f"{temp_dir}/test_locus.ld.npz"
    ldmap_file = f"{temp_dir}/test_locus.ldmap"
    
    sumstats.to_csv(sumstats_file, sep="\t", index=False)
    np.savez_compressed(ld_file, ld=ld_matrix.astype(np.float16))
    ld_map.to_csv(ldmap_file, sep="\t", index=False)
    
    # Create loci list file
    loci_list_data = {
        "locus_id": ["test_locus"],
        "prefix": ["EUR_test"],
        "popu": ["EUR"],
        "cohort": ["test"],
        "sample_size": [10000],
        "sumstats_path": [sumstats_file],
        "ld_path": [ld_file],
        "ldmap_path": [ldmap_file]
    }
    loci_df = pd.DataFrame(loci_list_data)
    loci_file = f"{temp_dir}/loci_list.txt"
    loci_df.to_csv(loci_file, sep="\t", index=False)
    
    return loci_file, temp_dir


def test_cli_abf_cojo():
    """Test ABF+COJO through CLI interface."""
    print("Creating test data...")
    loci_file, temp_dir = create_test_loci_file()
    output_dir = f"{temp_dir}/output"
    
    try:
        print("Testing CLI with abf_cojo tool...")
        
        # Test finemap command with abf_cojo
        from typer.testing import CliRunner
        runner = CliRunner()
        
        result = runner.invoke(app, [
            "finemap",
            loci_file,
            output_dir,
            "--tool", "abf_cojo",
            "--max-causal", "3",
            "--p-cutoff", "1e-6"
        ])
        
        if result.exit_code == 0:
            print("✅ CLI finemap command with abf_cojo succeeded!")
            
            # Check output files
            expected_output = f"{output_dir}/test_locus"
            if os.path.exists(f"{expected_output}/pips.txt") and os.path.exists(f"{expected_output}/creds.json"):
                print("✅ Output files generated successfully!")
                
                # Read and display some results
                pips_df = pd.read_csv(f"{expected_output}/pips.txt", sep="\t", header=None)
                print(f"Generated {len(pips_df)} PIPs")
                print(f"Top 3 PIPs: {pips_df.iloc[:3, 1].values}")
                
                import json
                with open(f"{expected_output}/creds.json") as f:
                    creds = json.load(f)
                print(f"Number of credible sets: {creds['n_cs']}")
                
                return True
            else:
                print("❌ Expected output files not found")
                return False
        else:
            print(f"❌ CLI command failed with exit code {result.exit_code}")
            print(f"Output: {result.output}")
            if result.exception:
                print(f"Exception: {result.exception}")
            return False
            
    except Exception as e:
        print(f"❌ Error testing CLI: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Clean up
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    success = test_cli_abf_cojo()
    if success:
        print("🎉 CLI integration test passed!")
    else:
        print("❌ CLI integration test failed!")