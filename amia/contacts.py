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
#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

from pymol import cmd
import pymol
import warnings, os, logging, ntpath, argparse, sys
warnings.filterwarnings("ignore")
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
import pandas as pd
from amia.get_raw_distances import *
import numpy as np
import datetime

class ContactAnalysis:

    def pre_processing(self, var_sys):
        """The mutations present within the structure are extracted from the file name and returned for analysis"""
        mutations = str(var_sys).split('_')[0:-1]
        return mutations

    def macromolecules(self, var_sys):
        """All molecules not identified as a part of the standard 20 amino acids are extrapolated from the PDB file
        for contact analysis. This is achieved by iterating over every 'residue' within the Biopython structure
        and comparing it to a dictionary of standard residues. Once all non-protein macromolecules have are identified
        they are returned for further processing and analysis"""
        three_to_one ={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                            'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                            'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        mol_contacts = {}
        non_prot_molecules = []
        parser = PDBParser(PERMISSIVE=1)
        structure_id = os.path.basename(var_sys)
        file_name = os.path.basename(var_sys)
        structure = parser.get_structure(structure_id, file_name)
        ppbuilder = PPBuilder()
        for model in structure:
            residues = model.get_residues()
            for residue in residues:
                residue = residue.get_resname()
                if residue not in three_to_one and residue not in non_prot_molecules:
                    non_prot_molecules.append(residue)
        for mol in non_prot_molecules:
            cmd.select(str(mol), 'resn ' + str(mol))
            cmd.select(str(mol) + '_a', 'resn ' + str(mol) + ' around 6 and not org and not resn A+C+G+T')
            cmd.distance('macromol_cont_' + str(mol) + str(mol) + '_a', str(mol), str(mol) + '_a', '3.6', mode='2')
            D = get_raw_distances('macromol_cont_' + str(mol) + str(mol) + '_a')
            if len(D) > 0:
                mol_contacts[mol] = (len(D))
        return mol_contacts

    def contact_calculator(self, res, res_pos):
        """Comparing differences between the number of contacts in the WT and variant structures at the point of
        residue introduction provides insight into how the folding of the secondary structure may be impacted.
        Here a residue is provided and selected as an object within PyMOL along with the surrounding residues.
        The script calculates the distances between the residue of interest and each of the surrounding residues
        and returns it as a dataframe with the chain, residue in contact with target residue and the number of contacts
        present between them"""
        contacts = {}
        three_letter = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN',
                        'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER',
                        'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        for i in cmd.get_chains():
            cmd.select('init_' + str(res) + '_' + str(i), 'resn ' + three_letter[res] +
                       ' and resi ' + str(res_pos) + ' and chain ' + str(i) + ' and not resn A+C+G+T')
            cmd.select('init_' + str(res) + '_' + str(i) + '_a', 'resn ' + three_letter[res] +
                       ' and resi ' + str(res_pos) + ' and chain ' + str(i) + ' around 6 and not resn A+C+G+T')
            seq = cmd.get_fastastr('init_' + str(res) + '_' + str(i) + '_a')
            seq3 = []
            if len(seq.split(str('\n'))) > 1:
                seq2 = seq.split(str('\n'))[1]
                for res_a in seq2:
                    if res_a not in seq3:
                        seq3.append(res_a)
            for res_a2 in seq3:
                cmd.select('res_a_' + str(res_a2), 'resn ' + three_letter[res_a2] + ' and chain ' + str(i))
                cmd.distance('dist_' + str(i) + '_' + str(res_a2), 'init_' + str(res) + '_' + str(i),
                             'res_a_' + str(res_a2), '3.6', mode='2')
                D = get_raw_distances('dist_' + str(i) + '_' + str(res_a2))
                if len(D) > 0:
                    contacts['{0},{1}'.format(str(i), str(res_a2))] = (len(D))
        if not bool(contacts) == False:
            return contacts

    def res_contacts(self, output_dir, pdb_file):
        """Each of the variant files generated previously are iterated over and the mutations within the files are
        extracted and seperated based on the WT residue, residue position and variant residue.
        The WT structure provided previously is then loaded for each variant and the WT residue and residue position are
        passed for contact determination. After the contacts for the WT are calculated the environment is completely
        emptied and the same process is carried out for each variant system and its respective variant residues.
        Once the WT and variant contacts have been calculated the script stores them as a dataframe which is
        exported with an HTML extension for comparative analysis"""
        variants_dataframe = pd.DataFrame()
        os.chdir(output_dir)
        for file in os.listdir():
            if file.endswith('.pdb'):
                for var in self.pre_processing(file):
                    var_dict = {'Variant': [str(var)]}
                    var_df = pd.DataFrame(var_dict)
                    edit_contact = {}
                    edit_contact2 = {}
                    initial_residue = var[0]
                    mutated_residue = var[len(var) - 1]
                    residue_pos = var[1:len(var) - 1]
                    cmd.reinitialize()
                    cmd.load(pdb_file)
                    wt_contacts = self.contact_calculator(initial_residue, residue_pos)
                    cmd.reinitialize()
                    cmd.load(file)
                    var_contacts = self.contact_calculator(mutated_residue, residue_pos)
                    if var_contacts != None and wt_contacts != None:
                        missing_wt = []
                        missing_var = []
                        for key in wt_contacts:
                            if key not in var_contacts:
                                missing_wt.append(key)
                        for key2 in var_contacts:
                            if key2 not in wt_contacts:
                                missing_var.append(key2)
                        for val in missing_wt:
                            var_contacts[val] = 0
                        for val2 in missing_var:
                            wt_contacts[val2] = 0
                        res1 = {key: val for key, val in
                                sorted(wt_contacts.items(), key=lambda ele: ele[0], reverse=False)}
                        res2 = {key: val for key, val in
                                sorted(var_contacts.items(), key=lambda ele: ele[0], reverse=False)}
                        for index in res1:
                            edit_contact.setdefault('WT Chain', []).append(index.split(',')[0])
                            edit_contact.setdefault('WT Surrounding Residue(s)', []).append(index.split(',')[1])
                            edit_contact.setdefault('WT No. of Contacts', []).append(res1[index])
                        for index2 in res2:
                            edit_contact2.setdefault('Variant Chain', []).append(index2.split(',')[0])
                            edit_contact2.setdefault('Variant Surrounding Residue(s)', []).append(index2.split(',')[1])
                            edit_contact2.setdefault('Variant No. of Contacts', []).append(res2[index2])
                        wt_contacts_df = pd.DataFrame.from_dict(edit_contact)
                        var_contacts_df = pd.DataFrame.from_dict(edit_contact2)
                        horizontal_concat = pd.concat([wt_contacts_df, var_contacts_df], axis=1)
                        horizontal_concat['Contact Difference'] = horizontal_concat.apply(
                            lambda x: x['Variant No. of Contacts'] - x['WT No. of Contacts'], axis=1)
                    variant_df = pd.concat([var_df, horizontal_concat], axis=1)
                    df1 = variant_df.replace(np.nan, '', regex=True)
                    variants_dataframe = pd.concat([variants_dataframe, df1], axis=0)
            result = variants_dataframe.to_html(index=False, border=2)
            text_file = open("mutation_index.html", "w")
            text_file.write(result)
        with open("mutation_index.html", 'r', encoding='utf-8') as file:
            data = file.readlines()
        data[2] = '    <tr style="text-align: center; background: #1abc9c;">\n'
        with open("mutation_index.html", 'w', encoding='utf-8') as file:
            file.writelines(data)
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
        file.close()

    def macromol_contacts(self, output_dir, pdb_file):
        """Changes in residues within the protein sequence may impact the contacts with other biomolecules within
        the system. Providing a comparison between the other macromolecules and the protein structure within in the
        WT and variant systems provides further insight to the potential interaction changes that may contribute to
        treatment failures"""
        variants_dataframe = pd.DataFrame()
        for file in os.listdir(output_dir):
            if file.endswith('.pdb'):
                var_dict = {'Variant System': [str(file[:-9])]}
                var_df = pd.DataFrame(var_dict)
                edit_contact = {}
                edit_contact2 = {}
                cmd.reinitialize()
                if os.getcwd() == os.path.dirname(pdb_file):
                    cmd.load(os.path.basename(pdb_file))
                elif os.getcwd() != os.path.dirname(pdb_file):
                    os.chdir(os.path.dirname(pdb_file))
                    cmd.load(os.path.basename(pdb_file))
                wt_mol_contacts = self.macromolecules(pdb_file)
                os.chdir(output_dir)
                cmd.reinitialize()
                cmd.load(file)
                var_mol_contacts = self.macromolecules(file)
                res1 = {key: val for key, val in
                        sorted(wt_mol_contacts.items(), key=lambda ele: ele[0], reverse=False)}
                res2 = {key: val for key, val in
                        sorted(var_mol_contacts.items(), key=lambda ele: ele[0], reverse=False)}
                for index in res1:
                    edit_contact.setdefault('Macromolecules', []).append(index)
                    edit_contact.setdefault('WT Surrounding Residues Contact(s)', []).append(res1[index])
                for index2 in res2:
                    edit_contact2.setdefault('Macromolecules', []).append(index2)
                    edit_contact2.setdefault('Variant Surrounding Residues Contact(s)', []).append(res2[index2])
                wt_contacts_df = pd.DataFrame.from_dict(edit_contact)
                var_contacts_df = pd.DataFrame.from_dict(edit_contact2)
                horizontal_concat = pd.merge(wt_contacts_df, var_contacts_df, on = "Macromolecules", how = "inner")
                horizontal_concat['Contact Difference'] = horizontal_concat.apply(lambda x: x['Variant Surrounding Residues Contact(s)'] - x['WT Surrounding Residues Contact(s)'], axis=1)
                variant_df = pd.concat([var_df, horizontal_concat], axis=1)
                df1 = variant_df.replace(np.nan, '', regex=True)
                variants_dataframe = pd.concat([variants_dataframe, df1], axis=0)
        result = variants_dataframe.to_html(index=False, border=2)
        text_file = open("macromolecules_index.html", "w")
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
        with open("macromolecules_index.html", 'r', encoding='utf-8') as file:
            data = file.readlines()
        data[2] = '    <tr style="text-align: center; background: #1abc9c;">\n'
        with open("macromolecules_index.html", 'w', encoding='utf-8') as file:
            file.writelines(data)
        file.close()
        return print('script complete')

def conts(*argv):
    print(datetime.datetime.now())
    print('conts')
    p = ContactAnalysis()
    #pymol.finish_launching(['pymol'])
    p.res_contacts(argv[0], argv[1])
    p.macromol_contacts(argv[0], argv[1])
    print(datetime.datetime.now())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file", help="Path to the PDB file that the mutations will be introduced into")
    parser.add_argument("--output_dir", help="Path to the directory that the variant PDB systems will be stored in")
    args = parser.parse_args()
    pdb_file = str(args.pdb_file)
    output_dir = str(args.output_dir)
    conts(output_dir, pdb_file)


