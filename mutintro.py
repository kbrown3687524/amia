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


from pymol import cmd
import pymol
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import warnings,  os, logging, ntpath, argparse, sys
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
import pandas as pd
import datetime

class MutationIntro:
    def __int__(self):
        path = os.getcwd()
        return path

    def mutant_processing(self, mutant_list):
        """Mutation table provided as a standard .csv file with each system
        as a column heading with the mutations present listed beneath accordingly.
        The script then adds each of the systems grouped with the subsequent mutation(s)
        into a dictionary for individual or simultaneuous mutant introduction"""
        single_list = {}
        multi_list = {}
        self.mutant_list = mutant_list
        list_path = os.path.dirname(self.mutant_list)
        if str(self.mutant_list).endswith('.csv'):
            mutant_df = pd.read_csv(self.mutant_list)
            for col in mutant_df.columns:
                df2 = pd.DataFrame(mutant_df[str(col)])
                system = str(col)
                for row in df2.iterrows():
                    if str(row[1][0]) != 'nan':
                        mutation = str(row[1][0])
                        single_list.setdefault(system, []).append(mutation)
                        multi_list.setdefault(system, []).append(mutation)
        return single_list, multi_list
        print(single_list, multi_list)

    def individual_introduction(self, pdb_file, output_dir, mutant_data):
        """The mutations are individually separated and introduced within a single structure afterwhich,
        the structure is saved to the output directory under the mutation that  was introduced.
        This is achieved by calling and implementing the mutagenesis wizard commands which are able to select
        a target residue at a specific chain and position and alter it to the variant residue by automatically selecting
        the rotamer with the least steric clashes"""

        three_letter = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN',
                        'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER',
                        'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        if os.getcwd() != os.path.dirname(pdb_file):
            os.chdir(os.path.dirname(pdb_file))
        structure_id = os.path.basename(pdb_file)
        file_name = os.path.basename(pdb_file)
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(structure_id, file_name)
        ppbuilder = PPBuilder()
        mutation_subset = self.mutant_processing(mutant_data)[0]
        pdb_path = os.path.dirname(pdb_file)
        for key in mutation_subset:
            os.chdir(pdb_path)
            for mutant in mutation_subset[str(key)]:
                initial_residue = mutant[0]
                mutated_residue = mutant[len(mutant) - 1]
                residue_pos = mutant[1:len(mutant) - 1]
                cmd.reinitialize()
                if os.getcwd() == os.path.dirname(pdb_file):
                    cmd.load(os.path.basename(pdb_file))
                elif os.getcwd() != os.path.dirname(pdb_file):
                    os.chdir(os.path.dirname(pdb_file))
                    cmd.load(os.path.basename(pdb_file))
                for i in cmd.get_chains(ntpath.basename(pdb_file).split('.')[0]):
                    cmd.select('var_' + mutant, 'resn ' + str(three_letter[initial_residue]) + ' and resi '
                               + residue_pos)
                    cmd.wizard('mutagenesis')
                    cmd.refresh_wizard()
                    cmd.get_wizard().do_select("/" + 'var_' + mutant + "//" + i + "/" + residue_pos + "/")
                    cmd.get_wizard().set_mode(three_letter[mutated_residue])
                    cmd.get_wizard().apply()
                cmd.set_wizard()
                if os.getcwd() == os.path.dirname(output_dir):
                    os.chdir(os.path.basename(output_dir))
                elif os.getcwd() != os.path.dirname(output_dir):
                    os.chdir(output_dir)
                cmd.save('{0}_{1}.pdb'.format(str(mutant), str('auto')))
        cmd.set_wizard("done")


    def simultaneous_introduction(self, pdb_file, output_dir, mutant_data):
        """The mutations are grouped according to their system and introduced within a single structure afterwhich,
        the structure is saved to the output directory under the mutations that  were introduced.
        This is achieved by calling and implementing the mutagenesis wizard commands which are able to select
        target residues at a specific chain and position and alter them to the respective variant residues
        by automatically selecting the rotamer with the least steric clashes"""
        three_letter = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN',
                           'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER',
                           'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        if os.getcwd() != os.path.dirname(pdb_file):
            os.chdir(os.path.dirname(pdb_file))
        structure_id = os.path.basename(pdb_file)
        file_name = os.path.basename(pdb_file)
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(structure_id, file_name)
        ppbuilder = PPBuilder()
        mutation_subset = self.mutant_processing(mutant_data)[1]
        for key in mutation_subset:
            mutant_list = ''
            cmd.reinitialize()
            if os.getcwd() == os.path.dirname(pdb_file):
                cmd.load(os.path.basename(pdb_file))
            elif os.getcwd() != os.path.dirname(pdb_file):
                os.chdir(os.path.dirname(pdb_file))
                cmd.load(os.path.basename(pdb_file))
            for mutant in mutation_subset[str(key)]:
                mutant_list = mutant_list + str(mutant) + '_'
                initial_residue = mutant[0]
                mutated_residue = mutant[len(mutant) - 1]
                residue_pos = mutant[1:len(mutant) - 1]
                for i in cmd.get_chains(ntpath.basename(pdb_file).split('.')[0]):
                   cmd.select('var_' + mutant, 'resn ' + str(three_letter[initial_residue]) + ' and resi '
                               + residue_pos)
                   cmd.wizard('mutagenesis')
                   cmd.refresh_wizard()
                   cmd.get_wizard().do_select("/" + 'var_' + mutant + "//" + i + "/" + residue_pos + "/")
                   cmd.get_wizard().set_mode(three_letter[mutated_residue])
                   cmd.get_wizard().apply()
            cmd.set_wizard()
            if os.getcwd() == os.path.dirname(output_dir):
               os.chdir(os.path.basename(output_dir))
            elif os.getcwd() != os.path.dirname(output_dir):
                os.chdir(output_dir)
            cmd.save('{0}_{1}.pdb'.format(str(mutant_list[:-1]), str('auto')))
        cmd.set_wizard("done")

    def foldx_emin(self, foldx_exe ,output_dir):
        """Once mutations have been successfully introduced, it is necessary to optimize the structures as the change
        of residues may result in increased sterich classhes that are not necessarily present within the biological
        structure. Here FoldX software is implemented through the command line to optimize each of the resultant variant
        structures. Please ensure that the appropriate FoldX for the OS has been downlaoded and unzipped within amia
        folder"""
        os.chdir(output_dir)
        for var_sys in os.listdir(output_dir):
            if var_sys.endswith('_auto.pdb'):
                os.system(str(foldx_exe) + ' --command=Optimize --pdb=' + str(var_sys) + ' --output-file=' + str(var_sys))

def main(*argv):
    print(datetime.datetime.now())
    p = MutationIntro()
    for file in os.listdir(os.getcwd()):
        if 'foldx' == str(file):
            os.chdir(file)
            for file2 in os.listdir():
                if 'foldx' in str(file2):
                    foldx_exe = '{0}/{1}'.format(str(os.getcwd()), str(file2))
    if argv[0] == 'multiple':
        p.simultaneous_introduction(argv[1], argv[2], argv[3])
    if argv[0] == 'single':
        p.individual_introduction(argv[1], argv[2], argv[3])
    #p.foldx_emin(foldx_exe, argv[2])
    print(datetime.datetime.now())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file", required=True, help="Path to the PDB file that the mutations will be introduced into")
    parser.add_argument("--output_dir", required=True, help="Path to the directory that the variant PDB systems will be stored in")
    parser.add_argument("--mutations", required=True, help="Path to the csv file containing the variant residues that need to be analysed")
    parser.add_argument('--mode', choices=['single', 'multiple'], required=False, default='single')
    args = parser.parse_args()
    pdb_file = str(args.pdb_file)
    output_dir = str(args.output_dir)
    mutant_data = str(args.mutations)
    mode = str(args.mode)
    main(mode, pdb_file, output_dir, mutant_data)

