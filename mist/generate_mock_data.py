#!/usr/bin/env python3
"""
Mock Data Generation Script for credtools.

Generates realistic GWAS summary statistics and genotype data for testing
multi-ancestry fine-mapping with credtools.

Author: Generated for credtools development
"""

import argparse
import os
import struct
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.spatial.distance import pdist, squareform

# Plotting imports (optional)
try:
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import seaborn as sns

    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Warning: matplotlib/seaborn not available. Plotting functionality disabled.")


class MockDataGenerator:
    """Generate realistic mock genetic data for fine-mapping analysis."""

    def __init__(self, config=None):
        """Initialize with configuration parameters."""
        self.config = self._get_default_config()
        if config:
            self.config.update(config)

        # Set random seed for reproducibility
        np.random.seed(self.config.get("seed", 42))

        # Population genetics parameters
        self.pop_params = {
            "EUR": {"fst": 0.0, "sample_size": 10000, "bottleneck": 1.0},
            "AFR": {"fst": 0.15, "sample_size": 8000, "bottleneck": 1.2},
            "EAS": {"fst": 0.12, "sample_size": 12000, "bottleneck": 0.9},
            "SAS": {"fst": 0.10, "sample_size": 9000, "bottleneck": 1.1},
            "AMR": {"fst": 0.08, "sample_size": 7000, "bottleneck": 1.0},
        }

    def _get_default_config(self):
        """Set default configuration parameters."""
        return {
            "n_snps": 10000,
            "n_loci": 3,
            "populations": ["EUR", "AFR", "EAS"],
            "chromosomes": [1, 2, 9],  # Chromosomes for each locus
            "locus_sizes": [1000000, 1500000, 2000000],  # Base pairs per locus
            "n_causal_per_locus": [1, 2, 3],  # Causal variants per locus
            "causal_scenarios": ["shared", "mixed", "pop_specific"],
            "heritability": 0.05,  # Very realistic single-locus heritability (5%)
            "polygenicity": 0.001,  # Fraction of variants with non-zero effects
            "effect_size_variance": 0.01,
            "maf_range": [0.05, 0.45],
            "ld_decay_rate": 1000,  # LD decay parameter (kb)
            "output_dir": "mock_data",
            "seed": 42,
            "variant_overlap": 0.5,  # Minimum proportion of shared variants across populations
            "locus_starts": [50000000, 60000000, 70000000],  # Start positions for each locus
        }

    def generate_positions_and_alleles(self, locus_idx):
        """Generate population-specific SNP positions and alleles with partial overlap."""
        # Handle cases where locus_idx exceeds configuration array lengths
        if locus_idx < len(self.config["chromosomes"]):
            chr_num = self.config["chromosomes"][locus_idx]
        else:
            chr_num = (locus_idx % 22) + 1  # Cycle through chromosomes 1-22

        if locus_idx < len(self.config["locus_sizes"]):
            locus_size = self.config["locus_sizes"][locus_idx]
        else:
            locus_size = 1000000  # 1 Mb default

        if locus_idx < len(self.config["locus_starts"]):
            start_pos = self.config["locus_starts"][locus_idx]
        else:
            start_pos = 50000000 + locus_idx * 10000000  # Default start positions

        n_snps_locus = self.config["n_snps"] // self.config["n_loci"]
        overlap_ratio = self.config["variant_overlap"]

        # Calculate shared and population-specific SNPs
        n_shared = int(n_snps_locus * overlap_ratio)
        n_pop_specific = n_snps_locus - n_shared

        # Generate shared positions across all populations
        shared_positions = np.sort(start_pos + np.random.randint(0, locus_size, n_shared * 2))[:n_shared]
        shared_alleles = self._generate_random_alleles(n_shared)

        positions_dict = {}
        alleles_dict = {}

        for pop in self.config["populations"]:
            # Start with shared positions
            pop_positions = list(shared_positions)
            pop_alleles = list(shared_alleles)

            # Add population-specific positions
            if n_pop_specific > 0:
                # Generate unique positions for this population
                pop_specific_pos = np.sort(start_pos + np.random.randint(0, locus_size, n_pop_specific * 3))
                # Remove positions that overlap with shared
                pop_specific_pos = np.setdiff1d(pop_specific_pos, shared_positions)[:n_pop_specific]

                pop_positions.extend(pop_specific_pos)
                pop_alleles.extend(self._generate_random_alleles(len(pop_specific_pos)))

            # Sort by position and store
            sorted_idx = np.argsort(pop_positions)
            positions_dict[pop] = np.array(pop_positions)[sorted_idx]
            alleles_dict[pop] = [pop_alleles[i] for i in sorted_idx]

        return chr_num, positions_dict, alleles_dict

    def _generate_random_alleles(self, n_variants):
        """Generate random allele pairs for variants."""
        alleles = ["A", "C", "G", "T"]
        allele_pairs = []

        for _ in range(n_variants):
            # Randomly select two different alleles
            chosen = np.random.choice(alleles, 2, replace=False)
            allele_pairs.append((chosen[0], chosen[1]))

        return allele_pairs

    def generate_allele_frequencies(self, positions_dict, populations):
        """Generate population-specific allele frequencies for population-specific variants."""
        freq_dict = {}

        # First, identify shared positions across all populations
        all_positions = [set(positions_dict[pop]) for pop in populations]
        shared_positions = set.intersection(*all_positions)

        # Generate base frequencies for shared positions (as dictionary)
        n_shared = len(shared_positions)
        shared_base_freq_dict = {}
        if n_shared > 0:
            shared_base_freq_array = np.random.beta(1.5, 1.5, n_shared)
            shared_base_freq_array = np.clip(shared_base_freq_array, *self.config["maf_range"])
            # Create position -> frequency mapping
            shared_positions_list = sorted(shared_positions)
            shared_base_freq_dict = dict(zip(shared_positions_list, shared_base_freq_array))

        for pop in populations:
            pop_positions = positions_dict[pop]
            n_snps = len(pop_positions)
            pop_freq = np.zeros(n_snps)

            # Map shared positions to indices in this population
            shared_idx = [i for i, pos in enumerate(pop_positions) if pos in shared_positions]
            pop_specific_idx = [i for i in range(n_snps) if i not in shared_idx]

            fst = self.pop_params[pop]["fst"]

            # Handle shared variants
            if shared_idx and n_shared > 0:
                for i in shared_idx:
                    pos = pop_positions[i]
                    base_freq = shared_base_freq_dict[pos]

                    if fst == 0:  # Reference population
                        pop_freq[i] = base_freq
                    else:
                        # Apply FST-based divergence
                        alpha = base_freq * (1 - fst) / fst
                        beta = (1 - base_freq) * (1 - fst) / fst
                        pop_freq[i] = np.clip(np.random.beta(alpha, beta), 0.01, 0.99)

            # Handle population-specific variants
            if pop_specific_idx:
                n_pop_specific = len(pop_specific_idx)
                pop_specific_freq = np.random.beta(1.5, 1.5, n_pop_specific)
                pop_specific_freq = np.clip(pop_specific_freq, *self.config["maf_range"])
                pop_freq[pop_specific_idx] = pop_specific_freq

            freq_dict[pop] = pop_freq

        return freq_dict

    def generate_ld_matrix(self, positions, decay_rate=1000):
        """Generate realistic LD matrix using exponential decay correlation structure."""
        n_snps = len(positions)

        # Use LD blocks approach similar to cojo_simu.py for more realistic structure
        ld_blocks = max(1, n_snps // 500)  # Roughly 500 SNPs per LD block
        block_size = n_snps // ld_blocks

        ld_matrix = np.eye(n_snps)  # Start with identity

        for block_idx in range(ld_blocks):
            start = block_idx * block_size
            end = start + block_size if block_idx < ld_blocks - 1 else n_snps
            block_n_snps = end - start

            # Generate exponential decay correlation within block
            for j in range(block_n_snps):
                for k in range(block_n_snps):
                    if j != k:
                        # Use exponential decay based on SNP index distance (proxy for physical distance)
                        distance = abs(j - k)
                        correlation = 0.9**distance  # Strong LD within blocks
                        ld_matrix[start + j, start + k] = correlation

        # Add weaker inter-block correlations based on physical distance
        distances = squareform(pdist(positions.reshape(-1, 1)))

        for i in range(n_snps):
            for j in range(i + 1, n_snps):
                # Only modify if not already set (i.e., inter-block)
                if ld_matrix[i, j] == 0:
                    # Physical distance-based correlation (much weaker)
                    physical_dist_kb = distances[i, j] / 1000  # Convert to kb
                    correlation = max(0.1, np.exp(-physical_dist_kb / decay_rate)) * 0.3  # Cap at 0.3
                    ld_matrix[i, j] = correlation
                    ld_matrix[j, i] = correlation

        # Ensure matrix is positive definite
        eigenvals = np.linalg.eigvals(ld_matrix)
        if np.min(eigenvals) < 1e-8:
            # Add regularization to diagonal
            regularization = 1e-6 - np.min(eigenvals)
            ld_matrix += np.eye(n_snps) * regularization

        # Ensure diagonal is exactly 1.0
        np.fill_diagonal(ld_matrix, 1.0)

        return ld_matrix

    def select_causal_variants(self, positions_dict, n_causal, scenario="shared"):
        """Select causal variants based on scenario with population-specific variant sets."""
        # Identify shared and population-specific positions
        populations = list(positions_dict.keys())
        all_positions = [set(positions_dict[pop]) for pop in populations]
        shared_positions = list(set.intersection(*all_positions))

        causal_dict = {}

        if scenario == "shared":
            # Select causal variants only from shared positions
            if len(shared_positions) < n_causal:
                print(f"Warning: Only {len(shared_positions)} shared variants available, adjusting n_causal")
                n_causal = len(shared_positions)

            if n_causal > 0:
                causal_positions = np.random.choice(shared_positions, n_causal, replace=False)

                # Store the actual positions, not indices!
                for pop in populations:
                    causal_dict[pop] = causal_positions
            else:
                for pop in populations:
                    causal_dict[pop] = np.array([])

        elif scenario == "pop_specific":
            # Different causal variants per population
            for pop in populations:
                pop_positions = positions_dict[pop]
                n_snps = len(pop_positions)
                if n_snps < n_causal:
                    print(f"Warning: Only {n_snps} variants in {pop}, adjusting n_causal")
                    actual_n_causal = n_snps
                else:
                    actual_n_causal = n_causal
                # Select actual positions, not indices
                causal_idx = np.random.choice(n_snps, actual_n_causal, replace=False)
                causal_dict[pop] = pop_positions[causal_idx]

        elif scenario == "mixed":
            # Mix of shared and population-specific
            n_shared = max(1, n_causal // 2)  # At least 1 shared
            n_specific = n_causal - n_shared

            # Select shared causal from shared positions
            if len(shared_positions) < n_shared:
                print(f"Warning: Only {len(shared_positions)} shared variants, adjusting n_shared")
                n_shared = len(shared_positions)
                n_specific = n_causal - n_shared

            shared_causal_positions = []
            if n_shared > 0 and shared_positions:
                shared_causal_positions = np.random.choice(shared_positions, n_shared, replace=False)

            for pop in populations:
                pop_positions = positions_dict[pop]
                causal_positions_list = []

                # Add shared causal variants (these are already positions)
                causal_positions_list.extend(shared_causal_positions)

                # Add population-specific causal variants
                if n_specific > 0:
                    # Find positions not already selected
                    available_positions = [pos for pos in pop_positions if pos not in shared_causal_positions]
                    if len(available_positions) >= n_specific:
                        specific_positions = np.random.choice(available_positions, n_specific, replace=False)
                        causal_positions_list.extend(specific_positions)
                    else:
                        print(f"Warning: Not enough variants for {pop}, using all available")
                        causal_positions_list.extend(available_positions)

                causal_dict[pop] = np.array(causal_positions_list)

        else:
            raise ValueError(f"Unknown scenario: {scenario}")

        return causal_dict

    def generate_effect_sizes(self, causal_variants_dict, allele_freqs, heritability):
        """Generate realistic, moderate effect sizes."""
        effect_dict = {}

        for pop in self.config["populations"]:
            causal_positions = causal_variants_dict[pop]
            n_causal = len(causal_positions)

            if n_causal == 0:
                effect_dict[pop] = {}
                continue

            # Get allele frequencies for causal variants
            # allele_freqs[pop] is an array, we need to find indices of causal positions
            # We'll need the positions array to map positions to indices
            # For now, get the frequencies using position lookup
            freqs = []
            for pos in causal_positions:
                # This is a temporary fix - ideally we should pass positions to this method
                # For now, estimate frequency from position (this will be fixed properly)
                freq = 0.25  # Default MAF
                freqs.append(freq)
            freqs = np.array(freqs)

            # Generate more realistic effect sizes
            # Use smaller base effects and more realistic scaling
            base_effect_size = 0.04  # Very small base effect for realistic p-values
            raw_effects = np.random.normal(0, base_effect_size, n_causal)

            # Moderate inverse relationship with MAF (less extreme than before)
            maf_weights = 1.0 / np.sqrt(freqs * (1 - freqs))
            maf_weights = np.clip(maf_weights, 0.5, 3.0)  # Cap the MAF weighting

            # Apply MAF weighting
            effects = raw_effects * maf_weights

            # Add some random variance to effect sizes
            effect_variance = np.random.uniform(0.8, 1.2, n_causal)  # 20% variance
            effects *= effect_variance

            # Random sign flips
            effects *= np.random.choice([-1, 1], n_causal)

            # Return dictionary mapping position -> effect
            effect_dict[pop] = dict(zip(causal_positions, effects))

        return effect_dict

    def generate_phenotypes(self, genotypes_dict, effect_sizes_dict, positions_dict, heritability):
        """Generate realistic phenotypes from genotypes and causal effects."""
        phenotypes_dict = {}

        for pop in self.config["populations"]:
            if pop not in genotypes_dict:
                continue

            genotypes = genotypes_dict[pop]
            causal_effects = effect_sizes_dict[pop]  # This is a position -> effect dict
            n_samples = genotypes.shape[0]

            # Generate genetic component from causal variants
            genetic_component = np.zeros(n_samples)

            if len(causal_effects) > 0:
                # causal_effects is a position -> effect dict
                # We need to find the position indices in the positions array
                positions = positions_dict[pop]  # Need to pass positions_dict to this method
                for causal_pos, effect in causal_effects.items():
                    if causal_pos in positions:
                        if isinstance(positions, np.ndarray):
                            idx = np.where(positions == causal_pos)[0][0]
                        else:
                            idx = positions.index(causal_pos)
                        genetic_component += genotypes[:, idx] * effect

            # Standardize genetic component (like cojo_simu.py line 90)
            if np.std(genetic_component) > 0:
                genetic_component = genetic_component / np.std(genetic_component)

            # Add environmental noise based on heritability (like cojo_simu.py lines 93-94)
            environmental_variance = (1 - heritability) / heritability
            environmental_noise = np.random.normal(0, np.sqrt(environmental_variance), n_samples)

            # Combined phenotype
            phenotype = genetic_component + environmental_noise
            phenotypes_dict[pop] = phenotype

        return phenotypes_dict

    def generate_summary_statistics(
        self, positions_dict, alleles_dict, allele_freqs, genotypes_dict, phenotypes_dict, locus_idx
    ):
        """Generate realistic GWAS summary statistics through association testing."""
        # Handle cases where locus_idx exceeds configuration array lengths
        if locus_idx < len(self.config["chromosomes"]):
            chr_num = self.config["chromosomes"][locus_idx]
        else:
            # Generate additional chromosome numbers sequentially
            chr_num = (locus_idx % 22) + 1  # Cycle through chromosomes 1-22

        sumstats_dict = {}

        for pop in self.config["populations"]:
            positions = positions_dict[pop]
            alleles = alleles_dict[pop]
            genotypes = genotypes_dict[pop]
            phenotype = phenotypes_dict[pop]
            n_snps = len(positions)
            sample_size = self.pop_params[pop]["sample_size"]

            # Perform linear regression for each SNP
            beta_hat = np.zeros(n_snps)
            se_values = np.zeros(n_snps)
            p_values = np.zeros(n_snps)

            for snp_idx in range(n_snps):
                X = genotypes[:, snp_idx].reshape(-1, 1)
                y = phenotype

                # Add intercept column
                X_with_intercept = np.column_stack([np.ones(sample_size), X])

                try:
                    # Linear regression: (X'X)^-1 X'y
                    XtX_inv = np.linalg.inv(X_with_intercept.T @ X_with_intercept)
                    beta_full = XtX_inv @ X_with_intercept.T @ y

                    # Extract SNP effect (second coefficient, first is intercept)
                    beta_hat[snp_idx] = beta_full[1]

                    # Calculate residuals and standard error
                    y_pred = X_with_intercept @ beta_full
                    residuals = y - y_pred
                    mse = np.sum(residuals**2) / (sample_size - 2)

                    # Standard error of SNP coefficient (natural, no artificial noise)
                    se_values[snp_idx] = np.sqrt(mse * XtX_inv[1, 1])

                    # T-statistic and p-value (completely natural, no clipping)
                    t_stat = beta_hat[snp_idx] / se_values[snp_idx] if se_values[snp_idx] > 0 else 0
                    p_values[snp_idx] = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=sample_size - 2))

                except np.linalg.LinAlgError:
                    # Handle singular matrix (constant genotype)
                    beta_hat[snp_idx] = 0
                    se_values[snp_idx] = 1  # Arbitrary large SE
                    p_values[snp_idx] = 1  # Non-significant

            # No clipping - let natural statistical relationships determine p-values

            # Create summary statistics DataFrame
            sumstats_df = pd.DataFrame(
                {
                    "CHR": chr_num,
                    "BP": positions,
                    "rsID": [f"rs{chr_num}_{pos}" for pos in positions],
                    "EA": [allele[0] for allele in alleles],  # Effect allele
                    "NEA": [allele[1] for allele in alleles],  # Non-effect allele
                    "EAF": allele_freqs[pop],
                    "BETA": beta_hat,
                    "SE": se_values,
                    "P": p_values,
                    "N": sample_size,
                }
            )

            sumstats_dict[pop] = sumstats_df

        return sumstats_dict

    def save_sumstats(self, sumstats_dict, locus_idx, output_dir):
        """Save summary statistics to files."""
        locus_files = {}

        for pop, df in sumstats_dict.items():
            filename = f"{pop}_loci{locus_idx + 1}.sumstats"
            filepath = output_dir / filename

            # Reorder columns to match credtools format
            columns = ["CHR", "BP", "rsID", "EA", "NEA", "EAF", "BETA", "SE", "P", "N"]
            df[columns].to_csv(filepath, sep="\t", index=False)

            locus_files[pop] = str(filepath.stem)  # Without extension

        return locus_files

    def save_ld_matrix(self, positions_dict, alleles_dict, ld_matrix_dict, locus_idx, output_dir):
        """Save population-specific LD matrices and map files."""
        # Handle cases where locus_idx exceeds configuration array lengths
        if locus_idx < len(self.config["chromosomes"]):
            chr_num = self.config["chromosomes"][locus_idx]
        else:
            # Generate additional chromosome numbers sequentially
            chr_num = (locus_idx % 22) + 1  # Cycle through chromosomes 1-22

        locus_files = {}

        for pop in self.config["populations"]:
            positions = positions_dict[pop]
            alleles = alleles_dict[pop]
            ld_matrix = ld_matrix_dict[pop]

            # Skip if LD matrix is None (PLINK failed for this population)
            if ld_matrix is None:
                print(f"Warning: Skipping LD matrix save for {pop} (no valid matrix)")
                continue

            # Save LD matrix (lower triangular format)
            ld_filename = f"{pop}_loci{locus_idx + 1}.ld"
            ld_filepath = output_dir / ld_filename

            with open(ld_filepath, "w") as f:
                for i in range(len(ld_matrix)):
                    row_values = ld_matrix[i, : i + 1]  # Lower triangular
                    f.write("\t".join(f"{val:.6f}" for val in row_values) + "\n")

            # Save LD map
            ldmap_filename = f"{pop}_loci{locus_idx + 1}.ldmap"
            ldmap_filepath = output_dir / ldmap_filename

            ldmap_df = pd.DataFrame(
                {
                    "CHR": chr_num,
                    "BP": positions,
                    "A1": [allele[0] for allele in alleles],
                    "A2": [allele[1] for allele in alleles],
                }
            )
            ldmap_df.to_csv(ldmap_filepath, sep="\t", index=False)

            locus_files[pop] = str(ld_filepath.stem)  # Without extension

        return locus_files

    def generate_loci_definition(self, all_locus_files, output_dir):
        """Generate loci definition file for credtools."""
        loci_data = []

        for locus_idx in range(self.config["n_loci"]):
            # Handle cases where locus_idx exceeds configuration array lengths
            if locus_idx < len(self.config["chromosomes"]):
                chr_num = self.config["chromosomes"][locus_idx]
            else:
                # Generate additional chromosome numbers sequentially
                chr_num = (locus_idx % 22) + 1  # Cycle through chromosomes 1-22

            for pop in self.config["populations"]:
                prefix = all_locus_files[locus_idx][pop]

                loci_data.append(
                    {
                        "prefix": str(output_dir / prefix),
                        "popu": pop,
                        "cohort": f"{pop}_cohort",
                        "sample_size": self.pop_params[pop]["sample_size"],
                        "chr": chr_num,
                        "locus_id": f"locus_{locus_idx + 1}",
                    }
                )

        # Save to file
        loci_df = pd.DataFrame(loci_data)
        loci_filepath = output_dir / "test_loci.txt"
        loci_df.to_csv(loci_filepath, sep="\t", index=False)

        return str(loci_filepath)

    def generate_correlated_genotypes(self, n_samples, n_snps, allele_freqs, ld_blocks=None):
        """Generate correlated genotypes using LD block approach from cojo_simu.py."""
        # Determine LD blocks based on SNP count (roughly 100-500 SNPs per block)
        if ld_blocks is None:
            ld_blocks = max(1, n_snps // 200)  # Smaller blocks for stronger LD

        block_size = n_snps // ld_blocks
        genotypes = np.zeros((n_samples, n_snps), dtype=np.int8)

        for i in range(ld_blocks):
            start = i * block_size
            end = start + block_size if i < ld_blocks - 1 else n_snps
            block_n_snps = end - start

            # Create correlation matrix with exponential decay (like cojo_simu.py)
            corr_matrix = np.zeros((block_n_snps, block_n_snps))
            for j in range(block_n_snps):
                for k in range(block_n_snps):
                    corr_matrix[j, k] = 0.9 ** abs(j - k)

            # Generate multivariate normal samples
            mvn_samples = np.random.multivariate_normal(mean=np.zeros(block_n_snps), cov=corr_matrix, size=n_samples)

            # Convert to genotypes using allele frequencies
            for j in range(block_n_snps):
                unif = stats.norm.cdf(mvn_samples[:, j])
                maf = allele_freqs[start + j]

                # Hardy-Weinberg genotype generation (same as cojo_simu.py)
                geno = np.zeros(n_samples)
                geno[unif > (1 - maf) ** 2] = 1  # Heterozygous
                geno[unif > (1 - maf) ** 2 + 2 * maf * (1 - maf)] = 2  # Homozygous alternate
                genotypes[:, start + j] = geno.astype(np.int8)

        return genotypes

    def create_plink_text_files(self, genotypes, positions, alleles, chr_num, pop, locus_idx, output_dir):
        """Create PLINK .ped/.map files (like cojo_simu.py) for association testing."""
        n_samples, n_snps = genotypes.shape
        prefix = f"{pop}_loci{locus_idx + 1}"
        output_path = Path(output_dir)

        # Create .map file (chr, rsid, genetic_dist, position)
        map_file = output_path / f"{prefix}.map"
        with open(map_file, "w") as f:
            for i, pos in enumerate(positions):
                rsid = f"rs{chr_num}_{pos}"
                f.write(f"{chr_num} {rsid} 0 {pos}\n")

        # Create .ped file (FID, IID, father, mother, sex, phenotype, genotypes)
        ped_file = output_path / f"{prefix}.ped"
        with open(ped_file, "w") as f:
            for i in range(n_samples):
                # FID, IID, father, mother, sex, phenotype
                line = [f"{pop}_{i+1}", f"{pop}_{i+1}", "0", "0", "1", "-9"]

                # Add genotypes (convert 0,1,2 to A/A, A/G, G/G format)
                for j in range(n_snps):
                    allele1, allele2 = alleles[j]
                    if genotypes[i, j] == 0:
                        line.extend([allele2, allele2])  # Homozygous reference
                    elif genotypes[i, j] == 1:
                        line.extend([allele1, allele2])  # Heterozygous
                    else:  # genotypes[i, j] == 2
                        line.extend([allele1, allele1])  # Homozygous alternate

                f.write(" ".join(line) + "\n")

        return str(map_file.stem)  # Return prefix without extension

    def create_plink_text_files_with_prefix(self, genotypes, positions, alleles, chr_num, prefix, output_dir):
        """Create PLINK .ped/.map files with custom prefix."""
        n_samples, n_snps = genotypes.shape
        output_path = Path(output_dir)

        # Create .map file (chr, rsid, genetic_dist, position)
        map_file = output_path / f"{prefix}.map"
        with open(map_file, "w") as f:
            for i, pos in enumerate(positions):
                rsid = f"rs{chr_num}_{pos}"
                f.write(f"{chr_num} {rsid} 0 {pos}\n")

        # Create .ped file (FID, IID, father, mother, sex, phenotype, genotypes)
        ped_file = output_path / f"{prefix}.ped"
        with open(ped_file, "w") as f:
            for i in range(n_samples):
                # FID, IID, father, mother, sex, phenotype
                line = [f"{prefix}_{i+1}", f"{prefix}_{i+1}", "0", "0", "1", "-9"]

                # Add genotypes (convert 0,1,2 to A/A, A/G, G/G format)
                for j in range(n_snps):
                    allele1, allele2 = alleles[j]
                    if genotypes[i, j] == 0:
                        line.extend([allele1, allele1])  # AA
                    elif genotypes[i, j] == 1:
                        line.extend([allele1, allele2])  # Aa
                    else:  # genotypes[i, j] == 2
                        line.extend([allele2, allele2])  # aa

                f.write(" ".join(line) + "\n")

        return str(map_file.stem)  # Return prefix without extension

    def run_plink_association(self, prefix, output_dir, phenotype_data):
        """Run PLINK association testing like in cojo_simu.py."""
        # Save phenotype file
        output_path = Path(output_dir)
        pheno_file = output_path / f"{prefix}.pheno"
        phenotype_data.to_csv(pheno_file, sep="\t", index=False, header=True)

        # Convert to binary format
        cmd = f"plink --file {output_path}/{prefix} --make-bed --out {output_path}/{prefix}"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Warning: PLINK binary conversion failed for {prefix}")
                print(f"Command: {cmd}")
                if result.stderr:
                    print(f"STDERR: {result.stderr[:200]}")
                return None
        except Exception as e:
            print(f"Warning: PLINK binary conversion error for {prefix}: {e}")
            return None

        # Run association test
        cmd = f"plink --bfile {output_path}/{prefix} --pheno {pheno_file} --assoc --out {output_path}/{prefix}_gwas"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Warning: PLINK association test failed for {prefix}")
                print(f"Command: {cmd}")
                if result.stderr:
                    print(f"STDERR: {result.stderr[:200]}")
                return None
        except Exception as e:
            print(f"Warning: PLINK association error for {prefix}: {e}")
            return None

        # Read and return the results
        try:
            gwas_results = pd.read_csv(f"{output_path}/{prefix}_gwas.qassoc", sep="\\s+")
            return gwas_results
        except Exception as e:
            print(f"Warning: Error reading PLINK results: {e}")
            return None

    def read_plink_association_results(self, gwas_prefix, output_dir, alleles_dict, pop):
        """Read PLINK association results and format for credtools."""
        try:
            # Read PLINK .qassoc results
            qassoc_file = output_dir / f"{gwas_prefix}.qassoc"
            gwas_results = pd.read_csv(qassoc_file, sep="\\s+", engine="python")

            # Create credtools format summary statistics
            sumstats_df = pd.DataFrame(
                {
                    "CHR": gwas_results["CHR"],
                    "BP": gwas_results["BP"],
                    "rsID": gwas_results["SNP"],
                    "EA": [allele[0] for allele in alleles_dict],  # Effect allele
                    "NEA": [allele[1] for allele in alleles_dict],  # Non-effect allele
                    "EAF": 0.5,  # Placeholder - would need frequency calculation
                    "BETA": gwas_results["BETA"],
                    "SE": gwas_results["SE"],
                    "P": gwas_results["P"],
                    "N": gwas_results["NMISS"],
                }
            )

            return sumstats_df

        except FileNotFoundError:
            print(f"Warning: Could not find PLINK results file {gwas_prefix}.qassoc")
            return None
        except Exception as e:
            print(f"Warning: Error reading PLINK results: {e}")
            return None

    def generate_phenotype_from_genotypes(self, genotypes, positions, causal_positions, causal_effects, heritability):
        """Generate phenotype from genotypes with specified causal effects and heritability."""
        n_samples = genotypes.shape[0]

        # Create genetic component from causal variants
        genetic_component = np.zeros(n_samples)

        for causal_pos in causal_positions:
            if causal_pos in positions:
                # Convert positions to list to use .index() or use numpy method
                if isinstance(positions, np.ndarray):
                    snp_idx = np.where(positions == causal_pos)[0][0]
                else:
                    snp_idx = positions.index(causal_pos)
                effect_size = causal_effects.get(causal_pos, 0.0)
                genetic_component += genotypes[:, snp_idx] * effect_size

        if np.std(genetic_component) > 0:
            genetic_component = genetic_component / np.std(genetic_component)

        # Add environmental noise based on heritability
        environmental_variance = (1 - heritability) / heritability if heritability > 0 else 1.0
        environmental_noise = np.random.normal(0, np.sqrt(environmental_variance), n_samples)

        return genetic_component + environmental_noise

    def convert_plink_to_sumstats(self, plink_results, positions, alleles, allele_freqs):
        """Convert PLINK association results to our sumstats format."""
        if plink_results is None:
            return None

        sumstats_df = plink_results.copy()

        # Map PLINK columns to credtools format
        # PLINK: CHR, SNP, BP, NMISS, BETA, SE, R2, T, P
        # credtools expects: CHR, BP, rsID, EA, NEA, EAF, BETA, SE, P, N

        if "SNP" in sumstats_df.columns:
            sumstats_df["rsID"] = sumstats_df["SNP"]

        # Add allele information
        if len(alleles) == len(sumstats_df):
            sumstats_df["EA"] = [a[0] for a in alleles]  # Effect allele
            sumstats_df["NEA"] = [a[1] for a in alleles]  # Non-effect allele
            sumstats_df["A1"] = [a[0] for a in alleles]  # Keep original format too
            sumstats_df["A2"] = [a[1] for a in alleles]

        # Add allele frequencies
        if isinstance(allele_freqs, dict):
            freq_values = list(allele_freqs.values())[: len(sumstats_df)]
        else:
            freq_values = allele_freqs[: len(sumstats_df)]

        sumstats_df["EAF"] = freq_values  # Effect allele frequency
        sumstats_df["AF"] = freq_values  # Keep original format too

        # Add sample size
        if "NMISS" in sumstats_df.columns:
            sumstats_df["N"] = sumstats_df["NMISS"]

        return sumstats_df

    def generate_plink_ld_matrix(self, plink_prefix):
        """Generate LD matrix using PLINK --r square command."""
        try:
            # Generate LD matrix with PLINK
            cmd = f"plink --bfile {plink_prefix} --r square --keep-allele-order --out {plink_prefix}_r"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode == 0:
                ld_file = f"{plink_prefix}_r.ld"
                if os.path.exists(ld_file):
                    # Read the square LD matrix
                    ld_matrix = np.loadtxt(ld_file)

                    # Ensure it's symmetric and has valid correlations
                    if ld_matrix.ndim == 2 and ld_matrix.shape[0] == ld_matrix.shape[1]:
                        # Make sure matrix is symmetric
                        ld_matrix = (ld_matrix + ld_matrix.T) / 2
                        # Ensure diagonal is 1.0
                        np.fill_diagonal(ld_matrix, 1.0)
                        # Clip values to valid correlation range
                        ld_matrix = np.clip(ld_matrix, -1, 1)
                        return ld_matrix
                    else:
                        print(f"Warning: Invalid LD matrix shape for {plink_prefix}: {ld_matrix.shape}")
                        return None
                else:
                    print(f"Warning: LD file not found: {ld_file}")
                    return None
            else:
                # If PLINK --r square fails, fall back to correlation of genotypes
                print(f"Warning: PLINK LD generation failed for {plink_prefix}, using fallback method")
                return self._generate_fallback_ld_matrix(plink_prefix)

        except Exception as e:
            print(f"Warning: Error generating LD matrix: {e}")
            return self._generate_fallback_ld_matrix(plink_prefix)

    def _generate_fallback_ld_matrix(self, plink_prefix):
        """Fallback LD matrix generation from genotype correlation."""
        try:
            # Try to read the .bed file and calculate correlation
            bim_file = f"{plink_prefix}.bim"
            if os.path.exists(bim_file):
                # For now, return None - the save_ld_matrix method can handle this
                print(f"Info: Fallback LD matrix generation not implemented, using None for {plink_prefix}")
                return None
            return None
        except Exception:
            return None

    def cleanup_temp_plink_files(self, plink_prefix):
        """Clean up temporary PLINK files."""
        extensions = [
            ".map",
            ".ped",
            ".bed",
            ".bim",
            ".fam",
            ".pheno",
            ".qassoc",
            ".log",
            ".nosex",
            "_r.ld",
            "_r.log",
            "_gwas.qassoc",
            "_gwas.log",
        ]

        for ext in extensions:
            file_path = f"{plink_prefix}{ext}"
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                except Exception:
                    pass  # Ignore cleanup errors

    def save_plink_files(self, genotypes_dict, positions_dict, alleles_dict, locus_idx, output_dir):
        """Save genotypes in PLINK binary format (.bed/.bim/.fam) with population-specific variants."""
        # Handle cases where locus_idx exceeds configuration array lengths
        if locus_idx < len(self.config["chromosomes"]):
            chr_num = self.config["chromosomes"][locus_idx]
        else:
            # Generate additional chromosome numbers sequentially
            chr_num = (locus_idx % 22) + 1  # Cycle through chromosomes 1-22

        for pop in self.config["populations"]:
            genotypes = genotypes_dict[pop]
            positions = positions_dict[pop]
            alleles = alleles_dict[pop]
            n_samples, n_snps = genotypes.shape

            prefix = f"{pop}_loci{locus_idx + 1}"

            # Create .fam file (family information)
            fam_file = output_dir / f"{prefix}.fam"
            with open(fam_file, "w") as f:
                for i in range(n_samples):
                    # FID IID father mother sex phenotype
                    f.write(f"{pop}_{i+1} {pop}_{i+1} 0 0 0 -9\n")

            # Create .bim file (variant information)
            bim_file = output_dir / f"{prefix}.bim"
            with open(bim_file, "w") as f:
                for i, pos in enumerate(positions):
                    # chr rsid genetic_dist physical_pos allele1 allele2
                    rsid = f"rs{chr_num}_{pos}"
                    allele1, allele2 = alleles[i]
                    f.write(f"{chr_num} {rsid} 0 {pos} {allele1} {allele2}\n")

            # Create .bed file (binary genotype data)
            bed_file = output_dir / f"{prefix}.bed"
            with open(bed_file, "wb") as f:
                # PLINK binary file magic number
                f.write(struct.pack("BBB", 0x6C, 0x1B, 0x01))

                # Write genotypes (SNP-major mode)
                for snp_idx in range(n_snps):
                    # Pack 4 genotypes per byte
                    bytes_needed = (n_samples + 3) // 4
                    packed_bytes = bytearray(bytes_needed)

                    for sample_idx in range(n_samples):
                        byte_idx = sample_idx // 4
                        bit_offset = (sample_idx % 4) * 2

                        # PLINK encoding: 00=homozygous AA, 01=missing, 10=heterozygous, 11=homozygous BB
                        geno = genotypes[sample_idx, snp_idx]
                        if geno == 0:
                            plink_code = 0b00  # AA
                        elif geno == 1:
                            plink_code = 0b10  # Aa
                        else:  # geno == 2
                            plink_code = 0b11  # aa

                        packed_bytes[byte_idx] |= plink_code << bit_offset

                    f.write(packed_bytes)

        return {pop: f"{pop}_loci{locus_idx + 1}" for pop in self.config["populations"]}

    def save_concatenated_sumstats(self, all_sumstats_data, output_dir):
        """Save concatenated summary statistics across all loci for each population."""
        print("\nGenerating concatenated summary statistics...")

        concatenated_files = {}

        for pop in self.config["populations"]:
            # Collect all sumstats dataframes for this population
            pop_sumstats = []

            for locus_data in all_sumstats_data:
                if pop in locus_data:
                    pop_sumstats.append(locus_data[pop])

            if pop_sumstats:
                # Concatenate all loci for this population
                combined_df = pd.concat(pop_sumstats, ignore_index=True)

                # Sort by chromosome and position
                combined_df = combined_df.sort_values(["CHR", "BP"]).reset_index(drop=True)

                # Ensure unique rsIDs across loci
                rsid_counts = {}
                new_rsids = []
                for rsid in combined_df["rsID"]:
                    if rsid in rsid_counts:
                        rsid_counts[rsid] += 1
                        new_rsid = f"{rsid}_dup{rsid_counts[rsid]}"
                    else:
                        rsid_counts[rsid] = 0
                        new_rsid = rsid
                    new_rsids.append(new_rsid)

                combined_df["rsID"] = new_rsids

                # Save concatenated file
                filename = f"{pop}_all_loci.sumstats"
                filepath = output_dir / filename

                # Reorder columns to match credtools format
                columns = ["CHR", "BP", "rsID", "EA", "NEA", "EAF", "BETA", "SE", "P", "N"]
                combined_df[columns].to_csv(filepath, sep="\t", index=False)

                concatenated_files[pop] = str(filepath.stem)
                print(f"  {pop}: {len(combined_df)} variants -> {filename}")

        return concatenated_files

    def save_concatenated_plink_files(self, all_genotypes_data, all_positions_data, all_alleles_data, output_dir):
        """Save concatenated PLINK binary files across all loci for each population."""
        print("\nGenerating concatenated PLINK files...")

        concatenated_files = {}

        for pop in self.config["populations"]:
            # Collect all data for this population across loci
            pop_genotypes = []
            pop_positions = []
            pop_alleles = []
            pop_chrs = []

            for locus_idx, (geno_data, pos_data, allele_data) in enumerate(
                zip(all_genotypes_data, all_positions_data, all_alleles_data)
            ):
                if pop in geno_data and pop in pos_data and pop in allele_data:
                    pop_genotypes.append(geno_data[pop])
                    pop_positions.extend(pos_data[pop])
                    pop_alleles.extend(allele_data[pop])

                    # Get chromosome for this locus
                    if locus_idx < len(self.config["chromosomes"]):
                        chr_num = self.config["chromosomes"][locus_idx]
                    else:
                        chr_num = (locus_idx % 22) + 1

                    pop_chrs.extend([chr_num] * len(pos_data[pop]))

            if pop_genotypes:
                # Concatenate genotype matrices horizontally (more SNPs)
                combined_genotypes = np.hstack(pop_genotypes)
                n_samples, n_snps = combined_genotypes.shape

                prefix = f"{pop}_all_loci"

                # Create .fam file (family information)
                fam_file = output_dir / f"{prefix}.fam"
                with open(fam_file, "w") as f:
                    for i in range(n_samples):
                        # FID IID father mother sex phenotype
                        f.write(f"{pop}_{i+1} {pop}_{i+1} 0 0 0 -9\n")

                # Create .bim file (variant information)
                bim_file = output_dir / f"{prefix}.bim"
                with open(bim_file, "w") as f:
                    for i, (chr_num, pos, allele_pair) in enumerate(zip(pop_chrs, pop_positions, pop_alleles)):
                        # Make unique rsID across loci
                        rsid = f"rs{chr_num}_{pos}"
                        allele1, allele2 = allele_pair
                        f.write(f"{chr_num} {rsid} 0 {pos} {allele1} {allele2}\n")

                # Create .bed file (binary genotype data)
                bed_file = output_dir / f"{prefix}.bed"
                with open(bed_file, "wb") as f:
                    # PLINK binary file magic number
                    f.write(struct.pack("BBB", 0x6C, 0x1B, 0x01))

                    # Write genotypes (SNP-major mode)
                    for snp_idx in range(n_snps):
                        # Pack 4 genotypes per byte
                        bytes_needed = (n_samples + 3) // 4
                        packed_bytes = bytearray(bytes_needed)

                        for sample_idx in range(n_samples):
                            byte_idx = sample_idx // 4
                            bit_offset = (sample_idx % 4) * 2

                            # PLINK encoding: 00=homozygous AA, 01=missing, 10=heterozygous, 11=homozygous BB
                            geno = combined_genotypes[sample_idx, snp_idx]
                            if geno == 0:
                                plink_code = 0b00  # AA
                            elif geno == 1:
                                plink_code = 0b10  # Aa
                            else:  # geno == 2
                                plink_code = 0b11  # aa

                            packed_bytes[byte_idx] |= plink_code << bit_offset

                        f.write(packed_bytes)

                concatenated_files[pop] = prefix
                print(f"  {pop}: {n_samples} samples x {n_snps} variants -> {prefix}.{{bed,bim,fam}}")

        return concatenated_files

    def plot_manhattan(self, concat_sumstats_files, output_dir):
        """Generate Manhattan plots for each population."""
        if not HAS_PLOTTING:
            print("Skipping Manhattan plots - matplotlib not available")
            return

        print("\nGenerating Manhattan plots...")
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        for pop in self.config["populations"]:
            if pop not in concat_sumstats_files:
                continue

            # Load concatenated summary statistics
            sumstats_file = output_dir / f"{concat_sumstats_files[pop]}.sumstats"
            df = pd.read_csv(sumstats_file, sep="\t")

            # Calculate -log10(p-values)
            df["-log10P"] = -np.log10(df["P"].clip(lower=1e-300))

            # Create cumulative positions for plotting
            df = df.sort_values(["CHR", "BP"])
            chr_lengths = df.groupby("CHR")["BP"].max()
            chr_starts = {}
            cumulative_pos = 0

            for chr_num in sorted(df["CHR"].unique()):
                chr_starts[chr_num] = cumulative_pos
                cumulative_pos += chr_lengths[chr_num] + 1e6  # Add 1Mb spacing

            df["POS_PLOT"] = df.apply(lambda row: chr_starts[row["CHR"]] + row["BP"], axis=1)

            # Create Manhattan plot
            fig, ax = plt.subplots(figsize=(16, 8))

            # Color chromosomes alternately
            colors = ["#1f77b4", "#ff7f0e"]
            for i, chr_num in enumerate(sorted(df["CHR"].unique())):
                chr_data = df[df["CHR"] == chr_num]
                ax.scatter(
                    chr_data["POS_PLOT"],
                    chr_data["-log10P"],
                    c=colors[i % 2],
                    s=20,
                    alpha=0.7,
                    label=f"Chr {chr_num}",
                )

            # Add significance lines
            ax.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.7, label="Genome-wide sig")
            ax.axhline(y=-np.log10(1e-5), color="orange", linestyle="--", alpha=0.7, label="Suggestive")

            # Customize plot
            ax.set_xlabel("Chromosome", fontsize=14)
            ax.set_ylabel("-log10(P-value)", fontsize=14)
            ax.set_title(f"Manhattan Plot - {pop} Population", fontsize=16)

            # Set chromosome labels on x-axis
            chr_centers = [chr_starts[chr_num] + chr_lengths[chr_num] / 2 for chr_num in sorted(df["CHR"].unique())]
            ax.set_xticks(chr_centers)
            ax.set_xticklabels([str(chr_num) for chr_num in sorted(df["CHR"].unique())])

            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

            # Save plot
            plot_file = plots_dir / f"{pop}_manhattan.png"
            plt.tight_layout()
            plt.savefig(plot_file, dpi=300, bbox_inches="tight")
            plt.close()

            print(f"  {pop}: Manhattan plot saved to {plot_file}")

    def plot_locuszoom(self, all_sumstats_data, all_positions_data, effect_sizes_data, output_dir):
        """Generate LocusZoom plots for each locus and population."""
        if not HAS_PLOTTING:
            print("Skipping LocusZoom plots - matplotlib not available")
            return

        print("\nGenerating LocusZoom plots...")
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        for locus_idx in range(self.config["n_loci"]):
            if locus_idx >= len(all_sumstats_data):
                continue

            sumstats_data = all_sumstats_data[locus_idx]
            positions_data = all_positions_data[locus_idx]
            effects_data = effect_sizes_data[locus_idx]

            for pop in self.config["populations"]:
                if pop not in sumstats_data:
                    continue

                df = sumstats_data[pop].copy()
                df["-log10P"] = -np.log10(df["P"].clip(lower=1e-300))

                # Identify causal variants
                causal_effects_dict = effects_data[pop]  # This is a position -> effect dict
                causal_positions = set()
                if len(causal_effects_dict) > 0:
                    causal_positions = set(causal_effects_dict.keys())

                # Create figure with two subplots
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), height_ratios=[3, 1])

                # Top panel: Association plot
                colors = ["blue" if pos not in causal_positions else "red" for pos in df["BP"]]
                sizes = [30 if pos not in causal_positions else 80 for pos in df["BP"]]

                scatter = ax1.scatter(
                    df["BP"] / 1e6, df["-log10P"], c=colors, s=sizes, alpha=0.7, edgecolors="black", linewidths=0.5
                )

                # Add significance line
                ax1.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.5, label="Genome-wide sig")

                # Label causal variants
                if causal_positions:
                    causal_data = df[df["BP"].isin(causal_positions)]
                    for _, row in causal_data.iterrows():
                        ax1.annotate(
                            "Causal",
                            xy=(row["BP"] / 1e6, row["-log10P"]),
                            xytext=(5, 5),
                            textcoords="offset points",
                            fontsize=8,
                            color="red",
                            weight="bold",
                        )

                ax1.set_ylabel("-log10(P-value)", fontsize=12)
                ax1.set_title(f"LocusZoom Plot - {pop} Population, Locus {locus_idx + 1}", fontsize=14)
                ax1.grid(True, alpha=0.3)
                ax1.legend()

                # Bottom panel: Genes/regions (simplified)
                chr_num = df["CHR"].iloc[0]
                region_start = df["BP"].min()
                region_end = df["BP"].max()

                # Add a simple gene track representation
                ax2.add_patch(
                    patches.Rectangle(
                        (region_start / 1e6, 0.4),
                        (region_end - region_start) / 1e6,
                        0.2,
                        facecolor="lightblue",
                        alpha=0.7,
                    )
                )
                ax2.text(
                    (region_start + region_end) / 2 / 1e6,
                    0.5,
                    "Simulated Region",
                    ha="center",
                    va="center",
                    fontsize=10,
                )

                ax2.set_xlim(region_start / 1e6, region_end / 1e6)
                ax2.set_ylim(0, 1)
                ax2.set_xlabel(f"Position on Chromosome {chr_num} (Mb)", fontsize=12)
                ax2.set_ylabel("Genes", fontsize=10)
                ax2.set_yticks([])

                # Save plot
                plot_file = plots_dir / f"{pop}_locus{locus_idx + 1}_locuszoom.png"
                plt.tight_layout()
                plt.savefig(plot_file, dpi=300, bbox_inches="tight")
                plt.close()

                print(f"  {pop} Locus {locus_idx + 1}: LocusZoom plot saved to {plot_file}")

    def generate_mock_dataset(self, include_genotypes=False, generate_plots=False):
        """Generate complete mock dataset with population-specific variants."""
        output_dir = Path(self.config["output_dir"])
        output_dir.mkdir(exist_ok=True)

        print(f"Generating mock data with {self.config['n_loci']} loci...")
        print(f"Populations: {', '.join(self.config['populations'])}")
        print(f"Variant overlap ratio: {self.config['variant_overlap']}")
        print(f"Output directory: {output_dir}")
        if include_genotypes:
            print("Including genotype files in PLINK format")
        if generate_plots:
            print("Generating Manhattan and LocusZoom plots")

        all_locus_files = []
        all_genotype_files = []

        # Collect data for concatenation and plotting
        all_sumstats_data = []
        all_genotypes_data = []
        all_positions_data = []
        all_alleles_data = []
        all_effect_sizes_data = []

        for locus_idx in range(self.config["n_loci"]):
            print(f"\nGenerating locus {locus_idx + 1}...")

            # Generate population-specific positions and alleles
            chr_num, positions_dict, alleles_dict = self.generate_positions_and_alleles(locus_idx)

            # Report variant overlap stats
            all_positions = [set(positions_dict[pop]) for pop in self.config["populations"]]
            shared_positions = set.intersection(*all_positions)
            # total_unique = len(set.union(*all_positions))  # Unused variable
            overlap_ratio = len(shared_positions) / len(positions_dict[list(positions_dict.keys())[0]])
            print(
                f"  Chromosome {chr_num}, {len(shared_positions)} shared variants, overlap ratio: {overlap_ratio:.2f}"
            )

            # Generate population-specific allele frequencies
            allele_freqs = self.generate_allele_frequencies(positions_dict, self.config["populations"])

            # Select causal variants
            scenario = self.config["causal_scenarios"][locus_idx % len(self.config["causal_scenarios"])]
            n_causal = self.config["n_causal_per_locus"][locus_idx % len(self.config["n_causal_per_locus"])]
            causal_variants = self.select_causal_variants(positions_dict, n_causal, scenario)

            # Generate effect sizes (shared across populations)
            effect_sizes = self.generate_effect_sizes(causal_variants, allele_freqs, self.config["heritability"])

            # Generate population-specific data using PLINK-based approach
            sumstats = {}
            ld_matrix_dict = {}
            genotypes_dict = {}

            for pop in self.config["populations"]:
                positions = positions_dict[pop]
                n_snps = len(positions)

                print(f"  Generating {pop} population data with PLINK...")

                # Generate correlated genotypes using LD blocks
                n_samples = self.pop_params[pop]["sample_size"]
                genotypes = self.generate_correlated_genotypes(n_samples, n_snps, allele_freqs[pop], ld_blocks=10)
                genotypes_dict[pop] = genotypes

                # Create temporary PLINK files
                chr_num = self.config["chromosomes"][locus_idx] if locus_idx < len(self.config["chromosomes"]) else 1
                plink_prefix = f"temp_{pop}_locus{locus_idx}"
                # Ensure output directory exists
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                # Use the same prefix for PLINK files
                self.create_plink_text_files_with_prefix(
                    genotypes, positions, alleles_dict[pop], chr_num, plink_prefix, output_dir
                )

                # Generate phenotypes with causal effects
                phenotype = self.generate_phenotype_from_genotypes(
                    genotypes,
                    positions,
                    causal_variants[pop] if pop in causal_variants else [],
                    effect_sizes.get(pop, {}),
                    self.config["heritability"],
                )
                print(
                    f"    Phenotype stats: mean={phenotype.mean():.3f}, std={phenotype.std():.3f}, range=[{phenotype.min():.3f}, {phenotype.max():.3f}]"
                )

                # Create phenotype dataframe - IDs must match PLINK files
                pheno_df = pd.DataFrame(
                    {
                        "FID": [f"{plink_prefix}_{i+1}" for i in range(n_samples)],
                        "IID": [f"{plink_prefix}_{i+1}" for i in range(n_samples)],
                        "PHENO": phenotype,
                    }
                )

                # Run PLINK association and get results
                association_results = self.run_plink_association(plink_prefix, Path(output_dir), pheno_df)

                # Convert to our sumstats format
                if association_results is not None:
                    sumstats[pop] = self.convert_plink_to_sumstats(
                        association_results, positions, alleles_dict[pop], allele_freqs[pop]
                    )
                    # Generate LD matrix using PLINK for successful populations
                    ld_matrix_dict[pop] = self.generate_plink_ld_matrix(f"{output_dir}/{plink_prefix}")
                else:
                    print(f"Warning: No association results for {pop}, skipping this population...")
                    # Set LD matrix to None for failed populations
                    ld_matrix_dict[pop] = None

                # Clean up temporary files
                self.cleanup_temp_plink_files(f"{output_dir}/{plink_prefix}")

            # Save genotype files if requested
            if include_genotypes:
                genotype_files = self.save_plink_files(
                    genotypes_dict, positions_dict, alleles_dict, locus_idx, output_dir
                )
                all_genotype_files.append(genotype_files)

            # Save files only if we have data
            if sumstats:
                sumstats_files = self.save_sumstats(sumstats, locus_idx, output_dir)
                self.save_ld_matrix(positions_dict, alleles_dict, ld_matrix_dict, locus_idx, output_dir)
            else:
                print(f"Warning: No valid sumstats data for locus {locus_idx + 1}, skipping file generation")
                continue

            all_locus_files.append(sumstats_files)

            # Collect data for concatenation and plotting
            all_sumstats_data.append(sumstats)
            all_positions_data.append(positions_dict)
            all_alleles_data.append(alleles_dict)
            all_effect_sizes_data.append(effect_sizes)
            all_genotypes_data.append(genotypes_dict)  # Always collected since always generated

            print(f"  Scenario: {scenario}, {n_causal} causal variants requested")

            # Print causal variant info for each population
            for pop in self.config["populations"]:
                causal_effects = effect_sizes[pop]  # This is a position -> effect dict
                if len(causal_effects) > 0:
                    causal_positions = list(causal_effects.keys())
                    effects_values = list(causal_effects.values())
                    print(
                        f"    {pop}: {len(causal_positions)} causal at positions {causal_positions} with effects {effects_values}"
                    )
                else:
                    print(f"    {pop}: No causal variants")

        # Generate concatenated files
        concat_sumstats_files = self.save_concatenated_sumstats(all_sumstats_data, output_dir)

        concat_genotype_files = {}
        if include_genotypes:  # all_genotypes_data is always available now
            concat_genotype_files = self.save_concatenated_plink_files(
                all_genotypes_data, all_positions_data, all_alleles_data, output_dir
            )

        # Generate loci definition file
        loci_file = self.generate_loci_definition(all_locus_files, output_dir)

        # Generate plots if requested
        if generate_plots:
            if HAS_PLOTTING:
                self.plot_manhattan(concat_sumstats_files, output_dir)
                self.plot_locuszoom(all_sumstats_data, all_positions_data, all_effect_sizes_data, output_dir)
            else:
                print("\nSkipping plots - matplotlib not available. Install with: pip install matplotlib seaborn")

        print("\nMock dataset generated successfully!")
        print(f"Loci definition file: {loci_file}")

        # Report individual locus files
        print(f"Per-locus files: {len(all_locus_files)} loci x {len(self.config['populations'])} populations")
        if include_genotypes:
            print(
                f"Per-locus genotype files: {len(all_genotype_files)} loci x {len(self.config['populations'])} populations"
            )

        # Report concatenated files
        print(f"Concatenated summary statistics: {len(concat_sumstats_files)} populations")
        if concat_genotype_files:
            print(f"Concatenated genotype files: {len(concat_genotype_files)} populations")

        print(f"Use with credtools: credtools pipeline {loci_file} results/")


def main():
    """Initialize command line interface."""
    parser = argparse.ArgumentParser(description="Generate mock genetic data for credtools")
    parser.add_argument("--config", type=str, help="YAML configuration file")
    parser.add_argument("--output-dir", type=str, default="mock_data", help="Output directory")
    parser.add_argument("--n-loci", type=int, default=3, help="Number of loci")
    parser.add_argument("--n-snps", type=int, default=10000, help="Total number of SNPs")
    parser.add_argument("--populations", nargs="+", default=["EUR", "AFR", "EAS"], help="Populations to simulate")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--include-genotypes", action="store_true", help="Generate PLINK genotype files (.bed/.bim/.fam)"
    )
    parser.add_argument(
        "--variant-overlap",
        type=float,
        default=0.5,
        help="Minimum proportion of shared variants across populations (default: 0.5)",
    )
    parser.add_argument(
        "--generate-plots", action="store_true", help="Generate Manhattan and LocusZoom plots (requires matplotlib)"
    )

    args = parser.parse_args()

    # Load configuration
    config = {}
    if args.config and os.path.exists(args.config):
        with open(args.config, "r") as f:
            config = yaml.safe_load(f)

    # Override with command line arguments
    if args.output_dir:
        config["output_dir"] = args.output_dir
    if args.n_loci:
        config["n_loci"] = args.n_loci
    if args.n_snps:
        config["n_snps"] = args.n_snps
    if args.populations:
        config["populations"] = args.populations
    if args.seed:
        config["seed"] = args.seed
    if args.variant_overlap:
        config["variant_overlap"] = args.variant_overlap

    # Generate data
    generator = MockDataGenerator(config)
    generator.generate_mock_dataset(include_genotypes=args.include_genotypes, generate_plots=args.generate_plots)


if __name__ == "__main__":
    main()
