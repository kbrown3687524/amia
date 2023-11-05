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


if __name__ == '__main__':
    base_path = os.getcwd()
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file", required=True,
                        help="Path to the PDB file that the mutations will be introduced into")
    parser.add_argument("--output_dir", required=True,
                        help="Path to the directory that the variant PDB systems will be stored in")
    parser.add_argument("--mutations", required=True,
                        help="Path to the csv file containing the variant residues that need to be analysed")
    parser.add_argument('--mode', choices=['single', 'multiple'], required=False, default='single')
    args = parser.parse_args()
    pdb_file = str(args.pdb_file)
    output_dir = str(args.output_dir)
    mutant_data = str(args.mutations)
    mode = str(args.mode)
    os.system('python mutintro.py --mode ' + str(mode) + ' --pdb_file ' + str(pdb_file) + ' --output_dir ' +
              str(output_dir) + ' --mutations ' + str(mutant_data))
    os.system('python contacts.py --pdb_file ' + str(pdb_file) + ' --output_dir ' + str(output_dir))
    os.system('python foldxana.py --pdb_file ' + str(pdb_file) + ' --output_dir ' + str(output_dir))
