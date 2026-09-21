#!/usr/bin/env python3

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
import csv
import pandas as pd


class MaestroRunner:

    def __init__(self, base_dir, config_file, pdb_file, output_dir):
        self.base_dir = Path(base_dir).resolve()
        self.config_file = Path(config_file).resolve()
        self.pdb_file = Path(pdb_file).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.maestro_exe = self.find_and_prepare_maestro()

    def find_and_prepare_maestro(self):
        configured_exe = os.environ.get("AMIA_MAESTRO")
        if configured_exe:
            maestro_exe = Path(configured_exe).expanduser().resolve()
            if maestro_exe.is_file():
                return maestro_exe
            raise FileNotFoundError(f"AMIA_MAESTRO does not point to a file: {maestro_exe}")

        platform_names = {
            "linux": "MAESTRO_linux_x64",
            "darwin": "MAESTRO_macos_x64",
            "win32": "MAESTRO_windows_x64",
        }
        platform_dir = platform_names.get(sys.platform)
        roots = [self.base_dir / "amia", self.base_dir]
        search_dirs = [root / platform_dir for root in roots if platform_dir]
        search_dirs += [root / "MAESTRO" for root in roots]
        search_dirs += [root / "MAESTRO_linux_x64" for root in roots]
        for maestro_dir in search_dirs:
            if not maestro_dir.is_dir():
                continue
            for file in maestro_dir.iterdir():
                if file.is_file() and file.name.lower().startswith("maestro"):
                    if os.name != "nt":
                        os.chmod(file, file.stat().st_mode | stat.S_IXUSR)
                    print(f"✅ Maestro executable prepared: {file}")
                    return file.resolve()
        raise FileNotFoundError(
            "Maestro executable not found. Install a native build and set AMIA_MAESTRO "
            "to its full path. The repository bundle is Linux-only."
        )

    def parse_mutations(self, mutation_lines):
        mutations = []
        if not mutation_lines:
            return mutations
        data_lines = mutation_lines[1:]
        for line in data_lines:
            cells = line.strip().split(',')
            for cell in cells:
                cell = cell.strip()
                if cell:
                    if re.match(r'^[A-Z]\d+[A-Z]$', cell):
                        mutations.append(cell)
                    else:
                        print(f"⚠️ Skipping invalid mutation format: {cell}")
        return mutations

    def run_maestro_for_mutations(self, mutations):
        for mut in mutations:
            wt_res = mut[0]
            pos = ''.join(filter(str.isdigit, mut))
            mutant_res = mut[-1]
            evalmut_str = f'{wt_res}{pos}{{{mutant_res}}}'
            resultfile = self.output_dir / f"{wt_res}{pos}{mutant_res}_result.txt"
            cmd = [
                str(self.maestro_exe),
                str(self.config_file),
                str(self.pdb_file),
                f'--evalmut={evalmut_str}',
                '--bu',
                f'--resultfile={resultfile}'
            ]
            print(f"🚀 Running: {' '.join(cmd)}")
            subprocess.run(cmd, cwd=self.output_dir, check=True)
            print(f"✅ Result saved to: {resultfile}")

    def run_maestro_combined_mutations(self, mutations, label="combined"):
        if not mutations:
            print("⚠️ No mutations provided for combined evaluation.")
            return
        evalmut_parts = []
        for mut in mutations:
            wt_res = mut[0]
            pos = ''.join(filter(str.isdigit, mut))
            mutant_res = mut[-1]
            evalmut_parts.append(f'{wt_res}{pos}{{{mutant_res}}}')
        combined_evalmut = ','.join(evalmut_parts)
        resultfile = self.output_dir / f"{label}_result.txt"
        cmd = [
            str(self.maestro_exe),
            str(self.config_file),
            str(self.pdb_file),
            f'--evalmut={combined_evalmut}',
            '--bu',
            f'--resultfile={resultfile}'
        ]
        print(f"🚀 Running combined mutation set: {' '.join(cmd)}")
        subprocess.run(cmd, cwd=self.output_dir, check=True)
        print(f"✅ Combined result saved to: {resultfile}")

    def parse_maestro_results(self):
        results_list = []
        headers = ['structure', 'seqlength', 'pH', 'mutation', 'score', 'delta_score', 'ddG', 'ddG_confidence']
        for file in self.output_dir.iterdir():
            if file.is_file() and file.name.endswith("_result.txt"):
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
                        record['structure'] = variant_name
                        results_list.append(record)
        print(f"✅ Parsed {len(results_list)} result entries from {self.output_dir}")
        return results_list

    def parse_combined_maestro_results(self):
        resultfile = self.output_dir / "combined_result.txt"
        if not resultfile.exists():
            print(f"⚠️ Combined result file not found: {resultfile}")
            return []
        headers = ['structure', 'seqlength', 'pH', 'mutation', 'score', 'delta_score', 'ddG', 'ddG_confidence']
        results_list = []
        print(f"📂 Parsing combined results from {resultfile}")
        with open(resultfile) as f:
            lines = [line.strip() for line in f if line.strip()]
            for line in lines:
                if line.startswith("#") or line.startswith("structure"):
                    continue
                parts = re.split(r'\s+', line)
                if len(parts) != len(headers):
                    print(f"⚠️ Unexpected line format in {resultfile.name}: {line}")
                    continue
                record = dict(zip(headers, parts))
                record['structure'] = "combined_variant"
                results_list.append(record)
        print(f"✅ Parsed {len(results_list)} entries from combined results")
        return results_list


def plot_variant_ddG(csv_path, sort_by_ddG=True, figsize=(14, 6), save_path=None):
    df = pd.read_csv(csv_path)
    variants = df[df['ddG'] != 0.00000].copy()
    if sort_by_ddG:
        variants.sort_values(by='ddG', ascending=False, inplace=True)
    plt.figure(figsize=figsize)
    plt.bar(variants['mutation'], variants['ddG'], color='steelblue', edgecolor='black')
    plt.xticks(rotation=90, ha='right')
    plt.axhline(0, color='gray', linestyle='--')
    plt.ylabel('ddG')
    plt.xlabel('Mutation')
    plt.title('ddG Values for Protein Variants Compared to Wildtype')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=1200)
        print(f"📊 Plot saved to: {save_path}")
    else:
        plt.show()


def plot_combined_ddG(csv_path, figsize=(10, 6), save_path=None):
    df = pd.read_csv(csv_path)

    # Filter out zero ddG entries
    variants = df[df['ddG'].astype(float) != 0.00000].copy()

    # Convert ddG to float if not already
    variants['ddG'] = variants['ddG'].astype(float)

    # Sort by ddG (optional)
    variants.sort_values(by='ddG', ascending=False, inplace=True)

    plt.figure(figsize=figsize)
    plt.bar(variants['structure'], variants['ddG'], color='darkorange', edgecolor='black')
    plt.xticks(rotation=45, ha='right')
    plt.axhline(0, color='gray', linestyle='--')
    plt.ylabel('ddG')
    plt.xlabel('System/Group')
    plt.title('ddG Values for Combined Mutation Sets by System')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=1200)
        print(f"📊 Combined ddG plot saved to: {save_path}")
    else:
        plt.show()



def infer_base_dir(path_with_amia):
    parts = Path(path_with_amia).parts
    if 'amia' not in parts:
        raise ValueError("Path does not contain 'amia' directory.")
    return Path(*parts[:parts.index('amia')])


def main():
    parser = argparse.ArgumentParser(description="Run Maestro with automatic base_dir and config detection")
    parser.add_argument("--pdb_file", required=True, help="PDB file for mutation")
    parser.add_argument("--output_dir", required=True, help="Directory to store output")
    parser.add_argument("--mutations", required=True, help="CSV file with mutation lines, comma-separated")
    parser.add_argument("--mode", choices=["single", "multiple"], default="single", help="Mutation mode: single or multiple (simultaneous mutations)")
    args = parser.parse_args()

    pdb_path = Path(args.pdb_file).resolve()
    output_dir = Path(args.output_dir).resolve()
    mutations_file = Path(args.mutations).resolve()
    mode = args.mode.lower()

    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "maestro_run.log"
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )

    print(f"📄 Logging output to: {log_file}")

    try:
        base_dir = infer_base_dir(pdb_path)
    except ValueError:
        try:
            base_dir = infer_base_dir(mutations_file)
        except ValueError:
            print("❌ Could not infer base_dir because 'amia' not found in pdb_file or mutations file path.")
            exit(1)

    platform_names = {
        "linux": "MAESTRO_linux_x64",
        "darwin": "MAESTRO_macos_x64",
        "win32": "MAESTRO_windows_x64",
    }
    maestro_roots = [base_dir / "amia", base_dir]
    maestro_dirs = [root / platform_names[sys.platform] for root in maestro_roots
                    if sys.platform in platform_names]
    maestro_dirs += [root / "MAESTRO" for root in maestro_roots]
    maestro_dirs += [root / "MAESTRO_linux_x64" for root in maestro_roots]
    config_file = next((directory / "config.xml" for directory in maestro_dirs
                        if (directory / "config.xml").exists()), None)
    if config_file is None:
        print("❌ Maestro config.xml not found. Install a native Maestro bundle.")
        exit(1)

    print("🕒 Start time:", datetime.datetime.now())
    print(f"ℹ️ Base directory inferred as: {base_dir}")
    print(f"ℹ️ Using config file: {config_file}")
    print(f"🧬 Mode: {mode}")

    # Load mutations BEFORE runner
    with open(mutations_file) as f:
        mutation_lines = f.readlines()

    # Initialize runner AFTER paths are resolved and mutation lines are loaded
    runner = MaestroRunner(base_dir, config_file, pdb_path, output_dir)
    mutations = runner.parse_mutations(mutation_lines)
    print(f"🔍 Parsed mutations: {mutations}")

    if mode == "multiple":
        df = pd.read_csv(mutations_file)
        for col in df.columns:
            mutations = df[col].dropna().astype(str).str.strip()
            mutations = [m for m in mutations if m]
            if not mutations:
                print(f"⚠️ Skipping column {col}: no valid mutations.")
                continue
            group_id = col.strip().replace(" ", "_").replace(":", "")
            print(f"🔄 Running group {group_id} with mutations: {mutations}")
            runner.run_maestro_combined_mutations(mutations, label=group_id)

        # Parse all *_result.txt files (group files included)
        results = runner.parse_maestro_results()

        csv_path = output_dir / "maestro_multiple_results_summary.csv"
    else:
        runner.run_maestro_for_mutations(mutations)
        results = runner.parse_maestro_results()
        csv_path = output_dir / "maestro_results_summary.csv"

    if results:
        df = pd.DataFrame(results)
        df.to_csv(csv_path, index=False)
        print(f"💾 Results summary saved to: {csv_path}")
        plot_path = output_dir / ("multiple_ddG_plot.png" if mode == "multiple" else "ddG_plot.png")
        if mode == "multiple":
            plot_combined_ddG(csv_path, save_path=plot_path)
        else:
            plot_variant_ddG(csv_path, save_path=plot_path)
    else:
        print("⚠️ No results parsed, CSV not created.")

    print("✅ Finished at:", datetime.datetime.now())


if __name__ == "__main__":
    main()
