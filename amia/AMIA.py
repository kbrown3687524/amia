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

def run_pipeline(pdb_file, output_dir, mutations, mode):
    """Run the variant modeling pipeline."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    script_dir = os.path.dirname(os.path.realpath(__file__))
    scripts = [
        os.path.join(script_dir, "mutintro.py"),
        os.path.join(script_dir, "contacts.py"),
        os.path.join(script_dir, "foldxana.py")
    ]

    for i, script in enumerate(scripts, start=1):
        click.echo(f"\n🔹 Running Script {i}: {script}")

        # Base command
        cmd = ["python", script, "--pdb_file", pdb_file, "--output_dir", output_dir]

        # Only add mutations and mode if it's mutintro.py
        if "mutintro.py" in script:
            cmd.extend(["--mutations", mutations, "--mode", mode])

        result = subprocess.run(cmd)

        if result.returncode != 0:
            click.echo(f"❌ Script failed: {script}")
            exit(1)

    click.echo("\n✅ Pipeline complete.")

if __name__ == "__main__":
    run_pipeline()

