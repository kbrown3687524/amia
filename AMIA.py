#!/usr/bin/env python3
#
# Project Title: "The development of an automated computational workflow to prioritize potential resistance variants identified in HIV
# Integrase Subtype C"
#
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at
# the University of the Western Cape.
#
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties
#
#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import warnings,  os, logging, ntpath, argparse, sys
import subprocess
import click
from pathlib import Path

@click.command()
@click.option('--pdb_file', required=True, type=click.Path(exists=True), help="Path to input PDB file")
@click.option('--output_dir', required=True, type=click.Path(), help="Output directory for variant systems")
@click.option('--mutations', required=True, type=click.Path(exists=True), help="CSV file with variant residues")
@click.option('--mode', type=click.Choice(['single', 'multiple']), default='single', help="Mutation mode")
@click.option('--smiles', default="", show_default=True,
              help="Ligand SMILES string (default is Dolutegravir)")
@click.option('--center', nargs=3, type=float, required=True,
              help="Docking box center coordinates (format: x y z)")
def run_pipeline(pdb_file, output_dir, mutations, mode, smiles, center):
    """
    Main AMIA CLI pipeline: mutation intro → docking → contact analysis → stability.
    """

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Correct script execution order:
    scripts = [
        os.path.join(script_dir, "mutintro.py"),
        os.path.join(script_dir, "contacts.py"),
        os.path.join(script_dir, "autodock.py"),
        os.path.join(script_dir, "maestroana.py"),
    ]

    for i, script in enumerate(scripts, start=1):
        script_name = os.path.basename(script)
        click.echo(f"\n🔹 Running Script {i}: {script_name}")

        # Base command
        cmd = ["python", script, "--pdb_file", pdb_file, "--output_dir", output_dir]

        if "mutintro.py" in script:
            cmd.extend(["--mutations", mutations, "--mode", mode])

        elif "autodock.py" in script:
            cmd.extend(["--smiles", smiles])
            cmd.extend(["--center"] + [str(c) for c in center])

        elif "maestroana.py" in script:
            cmd.extend(["--mutations", mutations])

        # Run the script
        result = subprocess.run(cmd)

        if result.returncode != 0:
            click.echo(f"❌ Script failed: {script_name}")
            exit(1)

    click.echo("\n✅ AMIA pipeline complete.")

if __name__ == "__main__":
    run_pipeline()

