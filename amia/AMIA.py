#!/usr/bin/env python3

import os
import subprocess
import sys
import click
import yaml
from pathlib import Path
from amia.report import generate_report

def get_last_completed(checkpoint_file):
    if checkpoint_file.exists():
        return checkpoint_file.read_text().strip()
    return None

def update_last_completed(checkpoint_file, step_name):
    checkpoint_file.write_text(step_name)


def _config_path(value, config_dir):
    path = Path(value).expanduser()
    return path if path.is_absolute() else config_dir / path

@click.command()
@click.option('--config', type=click.Path(exists=True), required=True, help="Path to YAML configuration file")
@click.option('--force', is_flag=True, default=False, help="Force rerun all steps, ignoring checkpoints")
def run_pipeline(config, force):
    """
    AMIA pipeline with smart checkpointing and config-based step control.
    Resumes from the last successful step, unless --force is used.
    """
    # Load config
    config_path = Path(config).expanduser().resolve()
    config_dir = config_path.parent
    with config_path.open('r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)

    # Required input
    pdb_file = str(_config_path(cfg['pdb_file'], config_dir))
    output_dir = _config_path(cfg['output_dir'], config_dir)
    mutations = str(_config_path(cfg['mutations'], config_dir))
    mode = cfg.get('mode', 'single')

    # Optional docking/passser inputs
    run_docking = cfg.get('run_docking', False)
    run_passer = cfg.get('run_passer', False)

    smiles = cfg.get('smiles', '')
    compound_name = cfg.get('compound_name', 'Ligand')
    center = cfg.get('center', [])

    passer_dir = str(_config_path(cfg.get('passer_dir', str(output_dir)), config_dir))
    passer_txt = cfg.get('passer_txt', 'passer_all_results.txt')
    passer_html = cfg.get('passer_html', 'passer_summary.html')
    passer_file = (str(_config_path(cfg['passer_file'], config_dir))
                   if cfg.get('passer_file') else '')

    # Setup paths
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = output_dir / ".checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    last_step_file = checkpoint_dir / "last_completed.txt"

    last_completed = None if force else get_last_completed(last_step_file)
    completed_steps = []
    if last_completed:
        completed_steps.append(last_completed)

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
        cmd = [sys.executable, str(full_script_path), "--pdb_file", pdb_file, "--output_dir", str(output_dir)]

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

        # Execute step
        try:
            result = subprocess.run(cmd, check=False)
        except OSError as error:
            generate_report(output_dir, config_path, "failed", completed_steps, str(error))
            raise click.ClickException(f"Unable to start step '{step_name}': {error}") from error
        if result.returncode != 0:
            click.echo(f"❌ Step failed: {step_name}")
            generate_report(
                output_dir,
                config_path,
                "failed",
                completed_steps,
                f"Step '{step_name}' returned exit code {result.returncode}.",
            )
            raise click.exceptions.Exit(result.returncode)

        update_last_completed(last_step_file, step_name)
        completed_steps.append(step_name)
        click.echo(f"✅ Step complete: {step_name} (checkpoint updated)")

    report_path = generate_report(output_dir, config_path, "complete", completed_steps)
    click.echo(f"📄 Run report written to: {report_path}")
    click.echo("\n✅ AMIA pipeline fully complete.")

if __name__ == "__main__":
    run_pipeline()
