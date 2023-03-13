# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG"
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

# Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
# Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import os, sys, subprocess
mdtask = '{0}\{1}'.format(str(os.getcwd()), 'MDTASK')
sys.path.append('{0}\{1}'.format(str(os.getcwd()), 'MDTASK'))
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR, XTC
from MDTASK import calc_network

def traj_min_rin(traj_folders):
    for folder in os.listdir(traj_folders):
        top_file = ''
        traj_file = ''
        os.chdir('{0}\{1}'.format(str(traj_folders), folder))
        for file in os.listdir():
            if file.endswith('.xtc') and 'small' not in str(file):
                traj_file = file
            elif file.endswith('.tpr'):
                top_file = file
        u = mda.Universe(top_file, traj_file)
        ag = u.select_atoms("name CA or name CB")
        #ag.write(str(folder) + '_traj_small.xtc', frames='all')
        ag.write("small.pdb")

def rin(traj_folders):
    os.chdir(mdtask)
    for folder in os.listdir(traj_folders) :
        top_file = ''
        traj_file = ''
        for file in os.listdir('{0}\{1}'.format(str(traj_folders), str(folder))):
            if file.endswith('small.xtc'):
                traj_file = '{0}\{1}\{2}'.format(str(traj_folders), str(folder), str(file))
            elif file.endswith('small.pdb'):
                top_file = '{0}\{1}\{2}'.format(str(traj_folders), str(folder), str(file))
        os.chdir('{0}/{1}'.format(str(traj_folders), folder))
        os.system('python ' + mdtask + '/calc_network.py --topology ' + str(top_file) +
                  ' --threshold 7.0 --step 100 --generate-plots  --calc-L --discard-graphs --lazy-load ' +
                  str(traj_file))


rin(r'C:\Users\ingwe\MD_results')
#traj_min_rin(r'C:\Users\ingwe\MD_results')