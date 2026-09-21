#!/usr/bin/env python3

import os
import subprocess
import click
import yaml
from pathlib import Path

def get_last_completed(checkpoint_file):
    if checkpoint_file.exists():
        return checkpoint_file.read_text().strip()
    return None

def update_last_completed(checkpoint_file, step_name):
    checkpoint_file.write_text(step_name)

@click.command()
@click.option('--config', type=click.Path(exists=True), required=True, help="Path to YAML configuration file")
@click.option('--force', is_flag=True, default=False, help="Force rerun all steps, ignoring checkpoints")
def run_pipeline(config, force):
    """
    AMIA pipeline with smart checkpointing and config-based step control.
    Resumes from the last successful step, unless --force is used.
    """
    # Load config
    with open(config, 'r') as f:
        cfg = yaml.safe_load(f)

    # Required input
    pdb_file = cfg['pdb_file']
    output_dir = Path(cfg['output_dir'])
    mutations = cfg['mutations']
    mode = cfg.get('mode', 'single')

    # Optional docking/passser inputs
    run_docking = cfg.get('run_docking', False)
    run_passer = cfg.get('run_passer', False)
    run_trajstat = cfg.get('run_trajstat', False)

    smiles = cfg.get('smiles', '')
    compound_name = cfg.get('compound_name', 'Ligand')
    center = cfg.get('center', [])

    passer_dir = cfg.get('passer_dir', str(output_dir))
    passer_txt = cfg.get('passer_txt', 'passer_all_results.txt')
    passer_html = cfg.get('passer_html', 'passer_summary.html')
    passer_file = cfg.get('passer_file', '')
    trajstat_systems = cfg.get('trajstat_systems', str(output_dir))
    trajstat_start_fr = cfg.get('trajstat_start_fr', 0)

    # Setup paths
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = output_dir / ".checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    last_step_file = checkpoint_dir / "last_completed.txt"

    last_completed = None if force else get_last_completed(last_step_file)

    script_dir = Path(__file__).resolve().parent

    # Ordered pipeline steps, mutintro/contacts/maestroana always first
    pipeline_steps = [
        ("mutintro", "mutintro.py"),
        ("contacts", "contacts.py"),
    ]

    # Add optional steps based on config
    run_maestroana = cfg.get('run_maestroana', False)
    if run_maestroana:
        pipeline_steps.append(("maestroana", "maestroana.py"))

    if run_passer:
        pipeline_steps.append(("passer", "autoallo.py"))
    if run_docking:
        pipeline_steps.append(("docking", "autodock.py"))
    if run_trajstat:
        pipeline_steps.append(("trajstat", "trajstat.py"))

    step_names = [name for name, _ in pipeline_steps]

    start_index = 0
    if last_completed and last_completed in step_names:
        start_index = step_names.index(last_completed) + 1

    if force:
        click.echo("⚠️ Force flag enabled: ignoring checkpoints and rerunning all steps.")

    # Main pipeline execution
    for step_index in range(start_index, len(pipeline_steps)):
        step_name, script_file = pipeline_steps[step_index]
        full_script_path = script_dir / script_file

        click.echo(f"\n🔹 Running step {step_index + 1}/{len(pipeline_steps)}: {step_name} ({script_file})")

        # Build command
        cmd = ["python", str(full_script_path), "--output_dir", str(output_dir)]

        if step_name != "trajstat":
            cmd += ["--pdb_file", pdb_file]

        if step_name == "mutintro":
            cmd += ["--mutations", mutations, "--mode", mode]
        elif step_name == "maestroana":
            cmd += ["--mutations", mutations, "--mode", mode]
        elif step_name == "passer":
            cmd += [
                "--pdb_dir", passer_dir,
                "--text_output", str(output_dir / passer_txt),
                "--html_output", str(output_dir / passer_html)
            ]
            if passer_file:
                cmd += ["--pdb_file", passer_file]
        elif step_name == "docking":
            cmd += [
                "--smiles", smiles,
                "--compound_name", compound_name,
                "--center"
            ] + [str(c) for c in center]
        elif step_name == "trajstat":
            cmd += [
                "--systems", str(trajstat_systems),
                "--start_fr", str(trajstat_start_fr),
            ]

        # Execute step
        result = subprocess.run(cmd)
        if result.returncode != 0:
            click.echo(f"❌ Step failed: {step_name}")
            exit(1)

        update_last_completed(last_step_file, step_name)
        click.echo(f"✅ Step complete: {step_name} (checkpoint updated)")

    click.echo("\n✅ AMIA pipeline fully complete.")

if __name__ == "__main__":
    run_pipeline()
