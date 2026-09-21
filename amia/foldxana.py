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
# Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
# Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import os, logging, argparse, sys
import shutil
import subprocess
import pandas as pd
import numpy as np
import datetime
from pathlib import Path

class FoldXAna:

    def foldx_stability(self, output_dir, pdb_file):
        output_dir = Path(output_dir).resolve()
        pdb_file = Path(pdb_file).resolve()
        configured_exe = os.environ.get("AMIA_FOLDX")
        foldx_exe = Path(configured_exe).expanduser() if configured_exe else None
        if not foldx_exe:
            foldx_dir = Path(__file__).resolve().parent.parent / "foldx"
            candidates = sorted(file for file in foldx_dir.glob("foldx*") if file.is_file())
            foldx_exe = candidates[0] if candidates else None
        if not foldx_exe:
            foldx_exe = Path(shutil.which("foldx") or "")
        if not foldx_exe or not foldx_exe.is_file():
            print("FoldX executable not found. Set AMIA_FOLDX to its full path.")
            return

        print(f"✅ Using FoldX executable at: {foldx_exe}")

        # Run WT stability
        command = [str(foldx_exe), "--command=Stability", f"--pdb={pdb_file.name}",
               f"--output-dir={output_dir}"]
        subprocess.run(command, cwd=pdb_file.parent, check=True)

        # Run variant stability
        for var_file in output_dir.glob('*_auto.pdb'):
            subprocess.run([str(foldx_exe), "--command=Stability", f"--pdb={var_file.name}"],
                           cwd=output_dir, check=True)


    def stability_changes(self, output_dir, pdb_file):
        output_dir = Path(output_dir).resolve()
        pdb_file = Path(pdb_file).resolve()
        var_stability = []
        wt_stability = []
        stability_diff = []
        files = []
        for file in output_dir.glob('*.fxout'):
            if pdb_file.stem in file.name:
                wt_stability.append(float(file.read_text().split('\t')[1]))
            else:
                files.append(file.name.split('_auto')[0])
                var_stability.append(float(file.read_text().split('\t')[1]))
        for i in var_stability:
            for j in wt_stability:
                diff = j - i
                stability_diff.append(diff)
        wt_stability_df = pd.DataFrame(wt_stability)
        var_stability_df = pd.DataFrame(var_stability)
        stability_diff_df = pd.DataFrame(stability_diff)
        files_df = pd.DataFrame(files)
        hortizontal_concat = pd.concat([wt_stability_df, files_df, var_stability_df, stability_diff_df],ignore_index=True, axis=1)
        hortizontal_concat.columns =["WT System Stability", "Variant System", "Variant System Stability", "System Stability Difference"]
        df1 = hortizontal_concat.replace(np.nan, '', regex=True)
        result = df1.to_html(index=False, border=2)
        report_path = output_dir / "stability_index.html"
        text_file = report_path.open("w", encoding="utf-8")
        text_file.write(result)
        text_file.write('\n<style>' +
                        '\n' + 'table {text-align: center;}' +
                        '\n' + 'table thead th {text-align: center;}' +
                        '\n' + 'table, th, td {' +
                        '\n' + '  border: 1px solid black;' +
                        '\n' + '  border-collapse: collapse;' +
                        '\n' + '}' +
                        '\n' + 'th, td {' +
                        '\n' + '  border-style: solid;' +
                        '\n' + '}' +
                        '\n' + '</style>')
        text_file.close()
        with report_path.open('r', encoding='utf-8') as file:
            data = file.readlines()
        data[2] = '    <tr style="text-align: center; background: #1abc9c;">\n'
        with report_path.open('w', encoding='utf-8') as file:
            file.writelines(data)
        file.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze FoldX stability changes for mutations")
    parser.add_argument("--pdb_file", required=True, help="Path to the PDB file that the mutations will be introduced into")
    parser.add_argument("--output_dir", required=True, help="Directory to store stability analysis results")
    args = parser.parse_args()

    print("🕒 Start time:", datetime.datetime.now())

    analyzer = FoldXAna()

    # Call methods with correct order: (pdb_file, output_dir)
    analyzer.foldx_stability(args.output_dir, args.pdb_file)
    analyzer.stability_changes(args.output_dir, args.pdb_file)

    print("✅ Finished at:", datetime.datetime.now())

if __name__ == "__main__":
    main()
