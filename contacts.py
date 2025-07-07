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
            cmd.select(str(mol) + '_a', 'resn ' + str(mol) + ' around 5 and not org and not resn A+C+G+T')
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
            cmd.select('init_' + str(res) + '_' + str(i), 'resn ' + three_letter[res] + ' and resi ' + str(res_pos) + ' and chain ' + str(i) + ' and not resn DA+DC+DG+DT+A+C+G+T')
            cmd.select('init_' + str(res) + '_' + str(i) + '_a', 'resn ' + three_letter[res] + ' and resi ' + str(res_pos) + ' and chain ' + str(i) + ' around 6 and not resn DA+DC+DG+DT+A+C+G+T')
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
        """Generates a comparative table of WT and Variant residue contacts per mutation site,
        exporting CSV and HTML with table borders, and mutation names without file suffix."""
        
        import os
        import pandas as pd
        
        data = []
        os.chdir(output_dir)

        for file in os.listdir():
            if file.endswith('.pdb'):
                mutation_name = file.replace('_auto.pdb', '')  # Clean mutation name

                for var in self.pre_processing(file):
                    initial_residue = var[0]
                    mutated_residue = var[-1]
                    residue_pos = var[1:-1]

                    # WT contact calculation
                    cmd.reinitialize()
                    cmd.load(pdb_file)
                    wt_contacts = self.contact_calculator(initial_residue, residue_pos) or {}

                    # Variant contact calculation
                    cmd.reinitialize()
                    cmd.load(file)
                    var_contacts = self.contact_calculator(mutated_residue, residue_pos) or {}

                    # Always include the mutated residue position
                    expected_key = None
                    for key in set(wt_contacts.keys()).union(var_contacts.keys()):
                        expected_key = key
                        break

                    # If no contacts were found, try to reconstruct the key using default info
                    if not expected_key:
                        # Use PyMOL to extract the chain/residue directly
                        cmd.reinitialize()
                        cmd.load(file)
                        resi_query = f"resi {residue_pos} and resn {mutated_residue}"
                        sel = cmd.select("mutation_site", resi_query)
                        chain = cmd.get_model("mutation_site").atom[0].chain if sel else "-"
                        residue = residue_pos
                        expected_key = f"{chain},{residue}"

                    # Collect all keys, even if WT and Variant have no contacts
                    all_keys = set(wt_contacts.keys()).union(var_contacts.keys())
                    if not all_keys:
                        all_keys = {expected_key}

                    for key in sorted(all_keys):
                        chain, residue = key.split(',')
                        wt_value = wt_contacts.get(key, 0)
                        var_value = var_contacts.get(key, 0)
                        contact_diff = var_value - wt_value

                        data.append({
                            'PDB File': mutation_name,
                            'WT Chain': chain if wt_value else '-',
                            'WT Residue': residue if wt_value else '-',
                            'WT Contacts': wt_value,
                            'Variant Chain': chain if var_value else '-',
                            'Variant Residue': residue if var_value else '-',
                            'Variant Contacts': var_value,
                            'Contact Difference': contact_diff
                        })

        # Create DataFrame
        df = pd.DataFrame(data)

        # Export to CSV
        csv_path = os.path.join(output_dir, 'contact_comparison_table.csv')
        df.to_csv(csv_path, index=False)

        # Export to HTML with borders and styling
        html_path = os.path.join(output_dir, 'contact_comparison_table.html')
        html_table = df.to_html(index=False, border=1)
        styled_html = f"""
        <html>
        <head>
        <style>
            table {{
                border-collapse: collapse;
                width: 100%;
            }}
            th, td {{
                border: 1px solid black;
                padding: 4px;
                text-align: left;
            }}
            th {{
                background-color: #f2f2f2;
            }}
        </style>
        </head>
        <body>
        {html_table}
        </body>
        </html>
        """
        with open(html_path, 'w') as f:
            f.write(styled_html)

        print(f"\n✅ CSV saved to: {csv_path}")
        print(f"✅ HTML saved to: {html_path}")


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

def main():
    parser = argparse.ArgumentParser(description="Analyze residue contacts.")
    parser.add_argument("--pdb_file", required=True, help="Path to the input PDB file")
    parser.add_argument("--output_dir", required=True, help="Directory to store results")
    args = parser.parse_args()

    print("🕒 Start time:", datetime.datetime.now())
    print("📄 PDB file:", args.pdb_file)
    print("📁 Output dir:", args.output_dir)

    analyzer = ContactAnalysis()
    analyzer.res_contacts(args.output_dir, args.pdb_file)
    analyzer.macromol_contacts(args.output_dir, args.pdb_file)

    print("✅ Finished at:", datetime.datetime.now())

if __name__ == "__main__":
    main()

