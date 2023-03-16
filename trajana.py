# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG"
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

# Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
# Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.tests.datafiles import PDB, GRO, XTC, TPR
from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np, os

class TrajAna:

    def align_traj(self, u):
        """For visualization purposes, the trajectory needs to be aligned to itself in order to provide an accurate
        interpretation of the results. The script below aligns the c-alphas atoms within the protein to the reference
        frame (assigned as the 1st frame (0) of the trajectory and removes the translational and rotational data."""
        average = align.AverageStructure(u, u, select='protein and name CA', ref_frame=0).run(verbose=True, step=10)
        ref = average.results.universe
        aligner = align.AlignTraj(u, ref, select='protein and name CA',in_memory=False).run(verbose=True, step =10)

    def rmsd_calc(self, traj_file, top_file):
        """A trajectory universe is created from the topology file and trajectory files generated from the
        production run of the simulations. The script accepts the repaired trajectory (no atom jumping and
        all atoms centered) and then calculates the RMSD of the protein backbone of the entire system and returns it as
         a dataframe for downstream processing."""
        u = mda.Universe(top_file, traj_file)
        R = rms.RMSD(u, u, select='backbone', ref_frame=0).run(verbose=True)
        dataframe = pd.DataFrame(R.results.rmsd, columns=['Frame','Time (ns)', str(traj_file)])
        dataframe['Time (ns) '] = dataframe['Time (ns)']/1000
        return (dataframe)

    def rmsd_comp(self, traj_folder, output_dir):
        """Each of the topology and repaired trajectory files are passed to the script from their respective subfolders
        and the rmsd_calc function is called to calculate the systems RMSD. Each of the system RMSDs are then appended
        to a dataframe for comparison and to be stored as a .csv file for the user to analyze accordingly"""
        df = pd.DataFrame()
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}\{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            rmsd_df = self.rmsd_calc(traj_file, top_file)
            rmsd_df.to_csv(str(folder) + '_rmsd.csv')
            if df.empty:
                df = pd.concat([df, rmsd_df.iloc[:,2]], axis=1)
            else:
                df = pd.concat([df, rmsd_df.iloc[:,2]], axis=1)
        os.chdir(output_dir)
        df.to_csv('variant_rmsd.csv')

    def rmsd_plotting(self, output_dir):
        """Here the .csv output from the rmsd_comp function is retrieved and converted to a figure using matplotlib,
        which displays each of the systems RMSD for ease of comparison between WT and variants"""
        os.chdir(output_dir)
        for file in os.listdir(output_dir):
            if file.endswith('_rmsd.csv'):
                df = pd.read_csv(file)
                df['Unnamed: 0'] = df['Unnamed: 0'] / 100
                for (columnName, columnData) in df.iteritems():
                    if str(columnName) != 'Time (ns)':
                        df.rename(columns={str(columnName): str(columnName).split('.')[0]}, inplace=True)
                df.rename(columns={'Unnamed: 0': 'Time (ns)'}, inplace=True)
                df.plot(x='Time (ns)', linewidth=0.75)
                plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
                plt.ylabel('RMSD ($\AA$)')
                plt.savefig('{0}\{1}{2}'.format(str(output_dir), str(file).split('.')[0], '.tiff'), bbox_inches='tight', dpi=600)


    def rmsf_calc(self, traj_file, top_file):
        """A trajectory universe is created from the topology file and trajectory files generated from the
        production run of the simulations. The script accepts the repaired trajectory (no atom jumping and
        all atoms centered) and then calculates the RMSF of each chain of the protein backbone of the entire system and
        returns it as a dataframe for downstream processing."""
        structure_df = pd.DataFrame()
        file_name = str(top_file).split('.')[0]
        u = mda.Universe(top_file, traj_file)
        self.align_traj(u)
        for file in os.listdir():
            if file.startswith('rmsfit_' ):
                u2 = mda.Universe(top_file, file)
                for seg in u2.segments:
                    if 'PRO' in str(seg):
                        chain_CA = u2.select_atoms('segid ' + str(str(seg).split(' ')[1][:-1]) + ' and name CA')
                        R2 = rms.RMSF(chain_CA).run(verbose=True)
                        chain_df = pd.DataFrame(R2.results.rmsf, columns = [file_name + ' ' + str(str(seg).split(' ')[1][:-1])])
                        structure_df = pd.concat([structure_df, chain_df], axis=1)
        return structure_df, len(structure_df.columns)

    def rmsf_comp(self, traj_folder, output_dir):
        """Each of the topology and repaired trajectory files are passed to the script from their respective subfolders
        and the rmsd_calc function is called to calculate the systems RMSF. Each of the systems chains RMSFs are then
        appended to a dataframe for comparison and to be stored as a .csv file for the user to analyze accordingly"""
        df = pd.DataFrame()
        df_len = ''
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}\{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            rmsf_df = self.rmsf_calc(traj_file, top_file)
            rmsf_df[0].to_csv(str(top_file).split('.')[0] + '_rmsf.csv')
            df_len = rmsf_df[1]
        counter = 1
        while counter <= df_len:
            df1 = pd.DataFrame()
            for folder in os.listdir(traj_folder):
                os.chdir('{0}\{1}'.format(str(traj_folder), folder))
                for file in os.listdir():
                    if file.endswith('_rmsf.csv'):
                        df = pd.read_csv(file)
                        df2 = df.iloc[:,counter]
                        df1 = pd.concat([df1,df2], axis=1)
            os.chdir(output_dir)
            df1.to_csv(str(counter) + '_chain_rmsf.csv')
            counter += 1

    def rmsf_plotting(self,output_dir):
        """Here the .csv output from the rmsf_comp function is retrieved and converted to a figure using matplotlib,
        which displays each of the chains RMSF within the systems for ease of comparison between WT and variants"""
        os.chdir(output_dir)
        for file in os.listdir(output_dir):
            if file.endswith('chain_rmsf.csv'):
                df = pd.read_csv(file)
                df.iloc[:,1:len(df)].plot(linewidth=0.75)
                plt.xlabel('Residue Number')
                plt.ylabel('RMSF ($\AA$)')
                plt.legend(df.iloc[:,1:len(df)], bbox_to_anchor=(1.0, 1.0), loc='upper left')
                plt.savefig('{0}\{1}{2}'.format(str(output_dir), str(file).split('.')[0], '.tiff'),
                            bbox_inches='tight', dpi=600)

    def radgyr(self, atomgroup, masses, total_mass=None):
        """The radius of gyration is crucial to determining the compactsness of the systems around each of the axis.
        Here, the defined function calculates the rgyr over each of the axes as well as a sum for all the arrays and
        returns it for downstream processing."""
        # coordinates change for each frame
        coordinates = atomgroup.positions
        center_of_mass = atomgroup.center_of_mass()
        # get squared distance from center
        ri_sq = (coordinates - center_of_mass) ** 2
        # sum the unweighted positions
        sq = np.sum(ri_sq, axis=1)
        sq_x = np.sum(ri_sq[:, [1, 2]], axis=1)  # sum over y and z
        sq_y = np.sum(ri_sq[:, [0, 2]], axis=1)  # sum over x and z
        sq_z = np.sum(ri_sq[:, [0, 1]], axis=1)  # sum over x and y
        # make into array
        sq_rs = np.array([sq, sq_x, sq_y, sq_z])
        # weight positions
        rog_sq = np.sum(masses * sq_rs, axis=1) / total_mass
        # square root and return
        return np.sqrt(rog_sq)

    def rgyr_calc(self,traj_file, top_file):
        """Each of the topology and trajectory files from each of the systems are retrieved and undergoes rgyr
        calculations which are then returned in the form of a dataframe for further comaprison."""
        u = mda.Universe(top_file, traj_file)
        protein = u.select_atoms('protein')
        rog = AnalysisFromFunction(self.radgyr, u.trajectory, protein, protein.masses,
                                   total_mass=np.sum(protein.masses)).run(verbose=True)
        rgyr_df = pd.DataFrame(rog.results['timeseries'], columns=[str(top_file).split(',')[0], 'x-axis', 'y-axis', 'z-axis'])
        rgyr_df.to_csv(str(top_file).split('.')[0] + '_rgyr.csv')
        return rgyr_df

    def rgyr_comp(self, traj_folder, output_dir):
        """EAch of the dataframes are iterated over with and the rgyr for all of the axis (stored as the 1st column in
        the previous dataframe is extracted and appended to a new dataframe which is then stored in the specified output
        directory as a .csv file"""
        df = pd.DataFrame()
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}\{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            print(top_file, traj_file)
            rgyr_df = self.rgyr_calc(traj_file, top_file)
            df = pd.concat([df, rgyr_df.iloc[:,0]], axis=1)
        os.chdir(output_dir)
        df.to_csv('multi_variant_rgyr.csv')

    def rgyr_plotting(self, output_dir):
        """The .csv file for the rgyr is read from within the output directory and each dataset present within the file
        is plotted on the same set of axes for ease of comparison. The resultant figure is then resaved within the
        output directory for further analysis."""
        os.chdir(output_dir)
        for file in os.listdir(output_dir):
            if file.endswith('_rgyr.csv'):
                df = pd.read_csv(file)
                df.iloc[:,0] = df.iloc[:,0]/100
                df.rename(columns= {'Unnamed: 0':'Time (ns)'}, inplace=True)
                df.plot(x='Time (ns)', linewidth=0.5)
                plt.xlabel('Time (ns)')
                plt.ylabel('Rgyr ($\AA$)')
                plt.legend(df.iloc[:, 1:len(df)], bbox_to_anchor=(1.0, 1.0), loc='upper left')
                plt.savefig('{0}\{1}{2}'.format(str(output_dir), str(file).split('.')[0], '_rgyr.tiff'),
                            bbox_inches='tight', dpi=600)


def main(systems, output_dir):
    p = TrajAna()
    #p.rmsd_comp(systems, output_dir)
    #p.rmsd_plotting(output_dir)
    #p.rmsf_comp(systems, output_dir)
    #p.rmsf_plotting(output_dir)
    #p.rgyr_comp(systems, output_dir)
    p.rgyr_plotting(output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--systems", help="Path to the trajectory files have been simulated. ")
    parser.add_argument("--output_dir", help="Path to the directory that the variant systems will be stored in")
    args = parser.parse_args()
    systems = str(args.systems)
    output_dir = str(args.output_dir)
    main(systems, output_dir)