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
import pandas as pd
import numpy as np
import datetime
from pathlib import Path

class FoldXAna:

    from pathlib import Path
import os

class FoldXAna:

    def foldx_stability(self, output_dir, pdb_file):
        # Get parent directory of this script (amia/)
        script_dir = Path(__file__).resolve().parent

        # FoldX lives in amia_project/foldx/
        foldx_dir = script_dir.parent / "foldx"

        # Try to find the foldx binary
        foldx_exe = None
        for file in foldx_dir.iterdir():
            if file.is_file() and "foldx" in file.name.lower():
                foldx_exe = file.resolve()
                break

        if not foldx_exe:
            print(f"❌ FoldX executable not found in: {foldx_dir}")
            return

        print(f"✅ Using FoldX executable at: {foldx_exe}")

        # Run WT stability
        os.chdir(Path(pdb_file).parent)
        os.system(f"{foldx_exe} --command=Stability --pdb={Path(pdb_file).name} --output-dir={output_dir}")

        # Run variant stability
        os.chdir(output_dir)
        for var_file in os.listdir(output_dir):
            if var_file.endswith('_auto.pdb'):
                os.system(f"{foldx_exe} --command=Stability --pdb={var_file}")


    def stability_changes(self, output_dir, pdb_file):
        os.chdir(output_dir)
        var_stability = []
        wt_stability = []
        stability_diff = []
        files = []
        os.chdir(output_dir)
        for file in os.listdir(output_dir):
            if file.endswith('.fxout'):
                if str(os.path.basename(pdb_file)).split('.')[0] in str(file):
                    wt_stability.append(float(open(file).read().split('\t')[1]))
                elif str(os.path.basename(pdb_file)).split('.')[0] not in str(file):
                    files.append(str(file).split('_auto')[0])
                    var_stability.append(float(open(file).read().split('\t')[1]))
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
        text_file = open("stability_index.html", "w")
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
        with open("stability_index.html", 'r', encoding='utf-8') as file:
            data = file.readlines()
        data[2] = '    <tr style="text-align: center; background: #1abc9c;">\n'
        with open("stability_index.html", 'w', encoding='utf-8') as file:
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
