#!/usr/bin/env python3
#
# Project Title: "The development of an automated computational workflow to prioritize potential resistance variants identified in HIV
# Integrase Subtype C"
#
# This script is developed for the fulfillment of a Master's degree at the South African National Bioinformatics Institute,
# the University of the Western Cape.
#
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme.
# Licensing and usage are governed by these organizations.
#
# Author: Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
# Author: Ruben Cloete (Supervisor) - Lecturer at SANBI (ruben@sanbi.ac.za)

import os
import subprocess
import click
import yaml
from pathlib import Path

@click.command()
@click.option('--config', type=click.Path(exists=True), required=True, help="Path to YAML configuration file")
def run_pipeline(config):
    """
    Main AMIA CLI pipeline: mutation intro → contact analysis → stability → (optional docking) → (optional PASSer).
    Configuration is loaded from a YAML file.
    """
    with open(config, 'r') as f:
        cfg = yaml.safe_load(f)

    pdb_file = cfg['pdb_file']
    output_dir = cfg['output_dir']
    mutations = cfg['mutations']
    mode = cfg.get('mode', 'single')
    smiles = cfg.get('smiles', '')
    compound_name = cfg.get('compound_name', 'Ligand')
    center = cfg['center']

    run_docking = cfg.get('run_docking', False)
    run_passer = cfg.get('run_passer', False)

    passer_dir = cfg.get('passer_dir', output_dir)
    passer_txt = cfg.get('passer_txt', 'passer_all_results.txt')
    passer_html = cfg.get('passer_html', 'passer_summary.html')
    passer_file = cfg.get('passer_file', '')

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    script_dir = os.path.dirname(os.path.realpath(__file__))

    scripts = [
        #os.path.join(script_dir, "mutintro.py"),
        #os.path.join(script_dir, "contacts.py"),
        #os.path.join(script_dir, "maestroana.py"),
    ]

    for i, script in enumerate(scripts, start=1):
        script_name = os.path.basename(script)
        click.echo(f"\n🔹 Running Script {i}: {script_name}")

        cmd = ["python", script, "--pdb_file", pdb_file, "--output_dir", output_dir]

        if "mutintro.py" in script:
            cmd.extend(["--mutations", mutations, "--mode", mode])
        elif "maestroana.py" in script:
            cmd.extend(["--mutations", mutations])

        result = subprocess.run(cmd)
        if result.returncode != 0:
            click.echo(f"❌ Script failed: {script_name}")
            exit(1)

    if run_passer:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        passer_script = os.path.join(script_dir, "autoallo.py")
        passer_dir = passer_dir or output_dir

        click.echo(f"\n🔹 Running PASSer on directory: {passer_dir}")
        passer_cmd = [
            "python", passer_script,
            "--pdb_file", pdb_file,
            "--pdb_dir", passer_dir,
            "--output_dir", output_dir,
            "--text_output", os.path.join(output_dir, passer_txt),
            "--html_output", os.path.join(output_dir, passer_html)
        ]

        
        result = subprocess.run(passer_cmd)
        if result.returncode != 0:
            click.echo("❌ PASSer script failed.")
            exit(1)

        click.echo("✅ PASSer API submission complete.")

        
    if run_docking:
        autodock_script = os.path.join(script_dir, "autodock.py")
        click.echo(f"\n🔹 Running AutoDock: {autodock_script}")
        dock_cmd = [
            "python", autodock_script,
            "--pdb_file", pdb_file,
            "--output_dir", output_dir,
            "--smiles", smiles,
            "--compound_name", compound_name,
            "--center"
        ] + [str(c) for c in center]

        result = subprocess.run(dock_cmd)
        if result.returncode != 0:
            click.echo("❌ AutoDock script failed.")
            exit(1)
        click.echo("✅ Docking complete.")

    

        if passer_file:
            passer_cmd.extend(["--pdb_file", passer_file])

        result = subprocess.run(passer_cmd)
        if result.returncode != 0:
            click.echo("❌ PASSer script failed.")
            exit(1)

        click.echo("✅ PASSer API submission complete.")

    click.echo("\n✅ AMIA pipeline complete.")

if __name__ == "__main__":
    run_pipeline()
