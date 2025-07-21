#!/usr/bin/env python3
#
# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV
# Integrase Subtype C and CRF02_AG"
#
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at
# the University of the Western Cape.
#
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties
#
#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import os
import argparse
import datetime
from pathlib import Path
import subprocess
import stat
import re
import pandas as pd
import matplotlib.pyplot as plt
import logging
import sys

class MaestroRunner:

    def __init__(self, base_dir, config_file, pdb_file, output_dir):
        self.base_dir = Path(base_dir).resolve()
        self.config_file = Path(config_file).resolve()
        self.pdb_file = Path(pdb_file).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.maestro_exe = self.find_and_prepare_maestro()

    def find_and_prepare_maestro(self):
        """Locate maestro executable in amia/MAESTRO_linux_x64, chmod u+x"""
        maestro_dir = self.base_dir / "amia" / "MAESTRO_linux_x64"
        if not maestro_dir.exists():
            raise FileNotFoundError(f"Directory not found: {maestro_dir}")

        for file in maestro_dir.iterdir():
            if file.is_file() and "maestro" in file.name.lower():
                st = file.stat()
                os.chmod(file, st.st_mode | stat.S_IXUSR)
                print(f"✅ Maestro executable prepared: {file}")
                return file.resolve()

        raise FileNotFoundError(f"Maestro executable not found in: {maestro_dir}")

    def parse_mutations(self, mutation_lines):
        """
        Parse your special CSV mutation format,
        skipping header, flatten mutations,
        and validate format WTres + position + Mutres (e.g., I84M)
        """
        mutations = []

        if not mutation_lines:
            return mutations

        # Skip header line (assumed first line)
        data_lines = mutation_lines[1:]

        for line in data_lines:
            cells = line.strip().split(',')
            for cell in cells:
                cell = cell.strip()
                if cell:
                    # Validate format like I84M: WT residue uppercase letter, digits, Mut residue uppercase letter
                    if re.match(r'^[A-Z]\d+[A-Z]$', cell):
                        mutations.append(cell)
                    else:
                        print(f"⚠️ Skipping invalid mutation format: {cell}")

        return mutations

    def run_maestro_for_mutations(self, mutations):
        for mut in mutations:
            wt_res = mut[0]                             # e.g., I
            pos = ''.join(filter(str.isdigit, mut))    # e.g., 84
            mutant_res = mut[-1]                       # e.g., M
            evalmut_str = f'{wt_res}{pos}{{{mutant_res}}}'

            # Construct result file path under output_dir
            resultfile = self.output_dir / f"{wt_res}{pos}{mutant_res}_result.txt"

            cmd = [
                str(self.maestro_exe),
                str(self.config_file),
                str(self.pdb_file),
                f'--evalmut={evalmut_str}',
                '--bu',
                f'--resultfile={resultfile}'  # Use = syntax here
            ]

            print(f"🚀 Running: {' '.join(cmd)}")
            subprocess.run(cmd, cwd=self.output_dir, check=True)
            print(f"✅ Result saved to: {resultfile}")

    def parse_maestro_results(self):
    
    #Parse newly created result files in output_dir and extract structured MAESTRO results,
    #explicitly mapping each entry to its mutation (variant) in the 'structure' column.
    #Returns:
    #    results_list: List of dictionaries with fields:
    #        structure, seqlength, pH, mutation, score, delta_score, ddG, ddG_confidence
    
        results_list = []
        headers = ['structure', 'seqlength', 'pH', 'mutation', 'score', 'delta_score', 'ddG', 'ddG_confidence']

        for file in self.output_dir.iterdir():
            if file.is_file() and file.name.endswith("_result.txt"):
                # Extract variant name from filename
                variant_name = file.stem.replace("_result", "")
                print(f"📂 Parsing: {file} as variant {variant_name}")

                with open(file) as f:
                    lines = [line.strip() for line in f if line.strip()]
                    for line in lines:
                        if line.startswith("#") or line.startswith("structure"):
                            continue
                        parts = re.split(r'\s+', line)
                        if len(parts) != len(headers):
                            print(f"⚠️ Unexpected line format in {file.name}: {line}")
                            continue

                        record = dict(zip(headers, parts))

                        # Overwrite 'structure' with the variant name for clarity
                        record['structure'] = variant_name

                        results_list.append(record)

        print(f"✅ Parsed {len(results_list)} result entries from {self.output_dir}")
        return results_list

def plot_variant_ddG(csv_path, sort_by_ddG=True, figsize=(14, 6), save_path=None):
    """
    Parses a CSV with mutation and ddG data, filters out wildtype, 
    and plots ddG values for variants.

    Parameters:
        csv_path (str): Path to the input CSV file.
        sort_by_ddG (bool): Whether to sort variants by ddG before plotting.
        figsize (tuple): Size of the plot (width, height).
        save_path (str or None): If provided, saves the plot to this path.
    """
    # Load CSV
    df = pd.read_csv(csv_path)

    # Filter out wildtype entries
    variants = df[df['ddG'] != 0.00000].copy()

    # Optional sorting
    if sort_by_ddG:
        variants.sort_values(by='ddG', ascending=False, inplace=True)

    # Plot
    plt.figure(figsize=figsize)
    plt.bar(variants['mutation'], variants['ddG'], color='steelblue', edgecolor='black')
    plt.xticks(rotation=90, ha='right')
    plt.axhline(0, color='gray', linestyle='--')
    plt.ylabel('ddG')
    plt.xlabel('Mutation')
    plt.title('ddG Values for Protein Variants Compared to Wildtype')
    plt.tight_layout()

    # Save or show
    if save_path:
        plt.savefig(save_path,  dpi=1200)
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()


def infer_base_dir(path_with_amia):
    """Extract base directory path before 'amia' folder"""
    parts = Path(path_with_amia).parts
    if 'amia' not in parts:
        raise ValueError("Path does not contain 'amia' directory.")
    idx = parts.index('amia')
    base_parts = parts[:idx]
    return Path(*base_parts)

def main():
    parser = argparse.ArgumentParser(description="Run Maestro with automatic base_dir and config detection")
    parser.add_argument("--pdb_file", required=True, help="PDB file for mutation")
    parser.add_argument("--output_dir", required=True, help="Directory to store output")
    parser.add_argument("--mutations", required=True, help="CSV file with mutation lines, comma-separated")
    args = parser.parse_args()

    pdb_path = Path(args.pdb_file).resolve()
    output_dir = Path(args.output_dir).resolve()
    mutations_file = Path(args.mutations).resolve()

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = output_dir / "maestro_run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

    class LoggerWriter:
        def __init__(self, level):
            self.level = level
        def write(self, message):
            if message.strip():  # Avoid blank lines
                self.level(message.strip())
        def flush(self): pass

    sys.stdout = LoggerWriter(logging.info)
    sys.stderr = LoggerWriter(logging.error)

    print(f"📄 Logging output to: {log_file}")

    # Infer base_dir from pdb_file path or fallback to mutations file path
    try:
        base_dir = infer_base_dir(pdb_path)
    except ValueError:
        try:
            base_dir = infer_base_dir(mutations_file)
        except ValueError:
            print("❌ Could not infer base_dir because 'amia' not found in pdb_file or mutations file path.")
            exit(1)

    config_file = base_dir / "amia" / "MAESTRO_linux_x64" / "config.xml"
    if not config_file.exists():
        print(f"❌ Config file not found: {config_file}")
        exit(1)

    print("🕒 Start time:", datetime.datetime.now())
    print(f"ℹ️ Base directory inferred as: {base_dir}")
    print(f"ℹ️ Using config file: {config_file}")

    with open(mutations_file) as f:
        mutation_lines = f.readlines()

    runner = MaestroRunner(base_dir, config_file, pdb_path, output_dir)
    mutations = runner.parse_mutations(mutation_lines)
    print(f"🔍 Parsed mutations: {mutations}")

    runner.run_maestro_for_mutations(mutations)

    # Parse and save results
    results = runner.parse_maestro_results()

    if results:
        df = pd.DataFrame(results)
        csv_path = output_dir / "maestro_results_summary.csv"
        df.to_csv(csv_path, index=False)
        print(f"💾 Results summary saved to: {csv_path}")

        # Generate high-resolution ddG plot
        plot_path = output_dir / "ddG_plot.png"
        plot_variant_ddG(csv_path, save_path=plot_path)
    else:
        print("⚠️ No results parsed, CSV not created.")

    print("✅ Finished at:", datetime.datetime.now())

if __name__ == "__main__":
    main()
