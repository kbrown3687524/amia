#!/usr/bin/env python3

import os
import subprocess
import click
import yaml
from pathlib import Path

def checkpoint_exists(checkpoint_dir, name):
    return (checkpoint_dir / f"{name}.done").exists()

def mark_checkpoint(checkpoint_dir, name):
    (checkpoint_dir / f"{name}.done").touch()

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
    output_dir = Path(cfg['output_dir'])
    mutations = cfg['mutations']
    mode = cfg.get('mode', 'single')
    smiles = cfg.get('smiles', '')
    compound_name = cfg.get('compound_name', 'Ligand')
    center = cfg['center']

    run_docking = cfg.get('run_docking', False)
    run_passer = cfg.get('run_passer', False)

    passer_dir = cfg.get('passer_dir', str(output_dir))
    passer_txt = cfg.get('passer_txt', 'passer_all_results.txt')
    passer_html = cfg.get('passer_html', 'passer_summary.html')
    passer_file = cfg.get('passer_file', '')

    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = output_dir / ".checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)

    script_dir = Path(__file__).resolve().parent

    scripts = [
        ("mutintro", script_dir / "mutintro.py"),
        ("contacts", script_dir / "contacts.py"),
        ("maestroana", script_dir / "maestroana.py"),
    ]

    for step_name, script in scripts:
        if checkpoint_exists(checkpoint_dir, step_name):
            click.echo(f"✅ Skipping {step_name} (checkpoint exists)")
            continue

        click.echo(f"\n🔹 Running Script: {script.name}")

        cmd = ["python", str(script), "--pdb_file", pdb_file, "--output_dir", str(output_dir)]

        if step_name == "mutintro":
            cmd.extend(["--mutations", mutations, "--mode", mode])
        elif step_name == "maestroana":
            cmd.extend(["--mutations", mutations])

        result = subprocess.run(cmd)
        if result.returncode != 0:
            click.echo(f"❌ Script failed: {script.name}")
            exit(1)

        mark_checkpoint(checkpoint_dir, step_name)
        click.echo(f"✅ Checkpoint created for {step_name}")

    if run_passer:
        step_name = "passer"
        if checkpoint_exists(checkpoint_dir, step_name):
            click.echo(f"✅ Skipping PASSer (checkpoint exists)")
        else:
            passer_script = script_dir / "autoallo.py"
            click.echo(f"\n🔹 Running PASSer on directory: {passer_dir}")

            passer_cmd = [
                "python", str(passer_script),
                "--pdb_file", pdb_file,
                "--pdb_dir", passer_dir,
                "--output_dir", str(output_dir),
                "--text_output", str(output_dir / passer_txt),
                "--html_output", str(output_dir / passer_html)
            ]

            if passer_file:
                passer_cmd.extend(["--pdb_file", passer_file])

            result = subprocess.run(passer_cmd)
            if result.returncode != 0:
                click.echo("❌ PASSer script failed.")
                exit(1)

            mark_checkpoint(checkpoint_dir, step_name)
            click.echo("✅ PASSer API submission complete.")

    if run_docking:
        step_name = "docking"
        if checkpoint_exists(checkpoint_dir, step_name):
            click.echo(f"✅ Skipping docking (checkpoint exists)")
        else:
            autodock_script = script_dir / "autodock.py"
            click.echo(f"\n🔹 Running AutoDock: {autodock_script.name}")

            dock_cmd = [
                "python", str(autodock_script),
                "--pdb_file", pdb_file,
                "--output_dir", str(output_dir),
                "--smiles", smiles,
                "--compound_name", compound_name,
                "--center"
            ] + [str(c) for c in center]

            result = subprocess.run(dock_cmd)
            if result.returncode != 0:
                click.echo("❌ AutoDock script failed.")
                exit(1)

            mark_checkpoint(checkpoint_dir, step_name)
            click.echo("✅ Docking complete.")

    click.echo("\n✅ AMIA pipeline complete.")

if __name__ == "__main__":
    run_pipeline()
