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
# Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
# Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import os, logging, argparse, sys
import pandas as pd
import numpy as np
import datetime


class FoldXAna:
    def foldx_stability(self, output_dir, pdb_file):
        print(os.getcwd())
        foldx_exe = ''
        for file in os.listdir(os.getcwd()):
            if 'foldx' == str(file):
                os.chdir(file)
                for file2 in os.listdir():
                    if 'foldx' in str(file2):
                        foldx_exe = '{0}{1}{2}'.format(str(os.getcwd()), str('/'), str(file2))
        print(foldx_exe)
        os.chdir(os.path.dirname(pdb_file))
        os.system(str(foldx_exe) + ' --command=Stability  --pdb=' + str(os.path.basename(pdb_file)) + ' --output-dir=' + str(output_dir))
        os.chdir(output_dir)
        for var_file in os.listdir(output_dir):
            if var_file.endswith('_auto.pdb'):
                os.system(str(foldx_exe) + ' --command=Stability  --pdb=' + str(var_file))

    def stability_changes(self, output_dir, pdb_file):
        os.chdir(output_dir)
        var_stability = []
        wt_stability = []
        stability_diff = []
        os.chdir(output_dir)
        for file in os.listdir(output_dir):
            if file.endswith('.fxout'):
                if str(os.path.basename(pdb_file)).split('.')[0] in str(file):
                    wt_stability.append(float(open(file).read().split('\t')[1]))
                elif str(os.path.basename(pdb_file)).split('.')[0] not in str(file):
                    var_stability.append(float(open(file).read().split('\t')[1]))
        for i in var_stability:
            for j in wt_stability:
                diff = j - i
                stability_diff.append(diff)
        wt_stability_df = pd.DataFrame(wt_stability)
        var_stability_df = pd.DataFrame(var_stability)
        stability_diff_df = pd.DataFrame(stability_diff)
        hortizontal_concat = pd.concat([wt_stability_df, var_stability_df, stability_diff_df],ignore_index=True, axis=1)
        hortizontal_concat.columns =["WT System Stability", "Variant System Stability", "System Stability Difference"]
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


def folds(*argv):
    print(datetime.datetime.now())
    p = FoldXAna()
    p.foldx_stability(argv[0], argv[1])
    p.stability_changes(argv[0], argv[1])
    print(datetime.datetime.now())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file", help="Path to the PDB file that the mutations will be introduced into")
    parser.add_argument("--output_dir", help="Path to the directory that the variant PDB systems will be stored in")
    args = parser.parse_args()
    pdb_file = str(args.pdb_file)
    output_dir = str(args.output_dir)
    folds(output_dir, pdb_file)
