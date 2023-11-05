import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.tests.datafiles import PDB, GRO, XTC, TPR
from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction, analysis_class
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np, os
import seaborn as sns
import datetime
import multiprocessing
from MDAnalysis.analysis import pca, align
import datetime
import plotly.express as px

class TrajStat:
    def rmsd_calc(self, top_file, traj_file):
        u = mda.Universe(top_file, traj_file)
        R = rms.RMSD(u, u, select='protein', ref_frame=0).run(verbose=True)
        dataframe = pd.DataFrame(R.results.rmsd, columns=['Frame','Time (ns)', str(os.path.basename(top_file).split('.')[0])])
        dataframe['Time (ns)'] = dataframe['Time (ns)']/1000
        dataframe.to_csv(top_file.split('.')[0] + '_rmsd.csv')

    def rmsf_calc(self, top_file, traj_file, start_fr):
        os.chdir(os.path.dirname(top_file))
        u = mda.Universe(top_file, traj_file)
        # Select the protein atoms
        u = mda.Universe(top_file, traj_file)
        # Calculate the RMSF
        for seg in u.segments:
            if 'PRO' in str(seg):
                protein = u.select_atoms('segid ' + str(str(seg).split(' ')[1][:-1]) + ' and name CA')
                rmsf = rms.RMSF(protein).run(verbose=True, start=int(start_fr))
                rmsf_df = pd.DataFrame(rmsf.results.rmsf)
                rmsf_df.to_csv('{3}/{1}_{2}'.format(top_file.split('.')[0], str(str(seg).split(' ')[1][:-1]), 'rmsf.csv', os.path.dirname(top_file)))

    def rmsd_plot(self, systems, output_dir):
        multi_df = pd.DataFrame()
        for system in os.listdir(systems):
            os.chdir(system)
            for file in os.listdir():
                if file.endswith('rmsd.csv'):
                    df = pd.read_csv(file)
                    if multi_df.empty ==True:
                        multi_df = pd.concat([multi_df, df.iloc[:,2:]], axis=1)
                    else:
                        multi_df = pd.merge(multi_df, df.iloc[:,2:], on='Time (ns)', how='inner')
            os.chdir(systems)
        multi_df.plot(x='Time (ns)', linewidth=0.75)
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.ylabel('RMSD ($\AA$)')
        plt.savefig('{0}/{1}{2}'.format(str(output_dir), 'multi_rmsd', '.tiff'),
        bbox_inches='tight', dpi=900)

    def rmsd_na_calc(self, top_file, traj_file):
        u = mda.Universe(top_file, traj_file)
        R = rms.RMSD(u, u, select='nucleic', ref_frame=0).run(verbose=True)
        dataframe = pd.DataFrame(R.results.rmsd, columns=['Frame','Time (ns)', str(os.path.basename(top_file).split('.')[0])])
        dataframe['Time (ns)'] = dataframe['Time (ns)']/1000
        dataframe.to_csv(top_file.split('.')[0] + '_rmsd_na.csv')
        
    def rmsd_na_plot(self, systems, output_dir):
        multi_df = pd.DataFrame()
        for system in os.listdir(systems):
            os.chdir(system)
            for file in os.listdir():
                if file.endswith('rmsd_na.csv'):
                    df = pd.read_csv(file)
                    if multi_df.empty ==True:
                        multi_df = pd.concat([multi_df, df.iloc[:,2:]], axis=1)
                    else:
                        multi_df = pd.merge(multi_df, df.iloc[:,2:], on='Time (ns)', how='inner')
            os.chdir(systems)
        multi_df.plot(x='Time (ns)', linewidth=0.75)
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.ylabel('RMSD ($\AA$)')
        plt.savefig('{0}/{1}{2}'.format(str(output_dir), 'multi_rmsd_na', '.tiff'),
        bbox_inches='tight', dpi=900)

    def rmsf_plot(self, systems, output_dir):
        prot_chains = []
        chain_file_pair = {}
        for system in os.listdir(systems):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(systems), system))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            u = mda.Universe(top_file, traj_file)
            for seg in u.segments:
                if 'PRO' in str(seg):
                    chain = str(seg).split(' ')[1].split('_')[-1][:-1]
                    if chain not in prot_chains:
                        prot_chains.append(chain)
        for rmsf_chain in prot_chains:
            for system in os.listdir(systems):
                os.chdir('{0}/{1}'.format(str(systems), system))
                for rmsf_file in os.listdir():
                    if str(rmsf_chain) in str(rmsf_file):
                        chain_file_pair.setdefault(rmsf_chain, []).append('{0}/{1}'.format(str(os.getcwd()), rmsf_file))
        for key in chain_file_pair:
            df = pd.DataFrame()
            for rmsf_data in chain_file_pair[key]:
                rmsf_df = pd.read_csv(rmsf_data)
                rmsf_df.rename(columns={'Unnamed: 0': 'Residue Number'}, inplace=True)
                rmsf_df.rename(columns={'0': str(os.path.dirname(rmsf_data).split('/')[-1])}, inplace=True)
                if df.empty == True:
                    df = rmsf_df
                else:
                    df = pd.concat([df, rmsf_df.iloc[:,1]], axis=1)
            df.plot(x='Residue Number', linewidth=0.75)
            plt.title('Multi-system RMSF of ' + str(key) + ' chain')
            plt.ylabel('RMSF ($\AA$)')
            plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            plt.savefig('{0}/{1}_{2}'.format(str(output_dir),  str(key),  '_rmsf.tiff'), bbox_inches='tight', dpi=900)


    def rgyr(self, atomgroup, masses, total_mass=None):
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

    def rgyr_calc(self, top_file, traj_file, start_fr):
        u = mda.Universe(top_file, traj_file)
        protein = u.select_atoms('protein')
        rog = AnalysisFromFunction(self.rgyr, u.trajectory, protein, protein.masses,
                                   total_mass=np.sum(protein.masses))
        rog.run(verbose=True, start=int(start_fr))
        rgyr_df = pd.DataFrame(rog.results['timeseries'],
                               columns=[str(os.path.basename(top_file).split('.')[0]), 'x-axis', 'y-axis', 'z-axis'])
        rgyr_df.to_csv(str(top_file).split('.')[0] + '_rgyr.csv')

    def rgyr_plot(self, systems, output_dir, start_fr):
        multi_df = pd.DataFrame()
        for system in os.listdir(systems):
            os.chdir('{0}/{1}'.format(str(systems), system))
            for file in os.listdir():
                if file.endswith('_rgyr.csv'):
                    df = pd.read_csv(file)
                    df.iloc[:, 0] = (df.iloc[:, 0] + int(start_fr)) / 10
                    df.rename(columns={'Unnamed: 0': 'Time (ns)'}, inplace=True)
                    if multi_df.empty == True:
                        multi_df = pd.concat([multi_df, df.iloc[:,:2]], axis=1)
                    else:
                        multi_df = pd.merge(multi_df, df.iloc[:,:2], on='Time (ns)', how='inner')
                    os.chdir(systems)
            multi_df.plot(x='Time (ns)', linewidth=0.75)
            plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            plt.ylabel('Rgyr ($\AA$)')
            plt.title('Multi-system Radius of gyration')
            plt.savefig('{0}/{1}{2}'.format(str(output_dir), 'multi_rgyr', '.tiff'),
                            bbox_inches='tight', dpi=900)
            

    def pca_calc(self, top_file, traj_file, start_fr):
        cumulated_PCA = {}
        print(datetime.datetime.now())
        u = mda.Universe(top_file, traj_file)
        backbone = u.select_atoms('backbone')
        pc = pca.PCA(u, select='backbone', align=True, mean=None, n_components=None, verbose=True)
        pc.run(verbose=True, start=int(start_fr))
        transformed = pc.transform(backbone, n_components=3)
        cov_direction = sum(pc.results.variance)
        for i in range(3):
            print(f"Cumulated variance: {pc.results.variance[i]:.3f}")
            cumulated_PCA[pc.results.variance[i]] = str(top_file.split('.')[0])
            print(cumulated_PCA)
        file = open(str(top_file.split('.')[0]) + '.txt', 'w')
        file.write(str(cumulated_PCA))
        file.write(str(cov_direction))
        file.close()
        df = pd.DataFrame(transformed, columns=['PC{}'.format(i + 1) for i in range(3)])
        
        df['Time (ps)'] = df.index * u.trajectory.dt
        df.to_csv(str(top_file).split('.')[0] + '_PCA_multi_MDA.csv')
        return df

    def pca_comp(self, traj_folder, output_dir, start_fr):
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            pca_df = self.pca_calc(top_file, traj_file, start_fr)
            pca_df.to_csv(str(top_file).split('.')[0] + '_PCA_MDA.csv')

    def pca_plotting(self, traj_folder, output_dir, start_fr):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            broad_df = pd.DataFrame()
            for file in os.listdir():
                if file.endswith('_PCA_MDA.csv'):
                    df = pd.read_csv(file)
                    df = df.drop(columns=['PC3'])
                    df.rename(columns={'Time (ps)': 'Time (ns)'}, inplace=True)
                    df['Time (ns)'] = df['Time (ns)'] / 1000
                    df = df.iloc[int(start_fr):, 1:]
                    ax1.scatter(df['PC1'], df['PC2'], s=20, marker="o", label=str(folder))
                    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
                    plt.title('Multi-system PCA')
                    plt.xlabel('PC1 ($\AA$)')
                    plt.ylabel('PC2 ($\AA$)')
        plt.savefig('{0}/{1}{2}'.format(str(output_dir), 'multi_PCA2', '.tiff'),
                    bbox_inches='tight', dpi=900)

    def pca_time_cluster_plotting(self, traj_folder, output_dir, start_fr):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            broad_df = pd.DataFrame()
            for file in os.listdir():
                if file.endswith('_PCA_MDA.csv'):
                    df = pd.read_csv(file)
                    df = df.drop(columns=['PC3'])
                    df.rename(columns={'Time (ps)': 'Time (ns)'}, inplace=True)
                    df['Time (ns)'] = df['Time (ns)'] / 1000
                    df = df.iloc[int(start_fr):, 1:]
                    fig = px.scatter(df, x='PC1', y='PC2', color=df["Time (ns)"])
                    fig.write_image('{0}/{1}{2}'.format(str(output_dir), str(folder), '_PCA_timeclustering.png'),
                                    width=1 * 900,
                                    height=1 * 900, scale=1)
    def contacts_within_cutoff(self, u, group_a, group_b, radius=3.5):
        """The MDAnalysis package allows a user to iterate through every frame of the trajetcory and calculate the
        number of contacts between two selections of atoms within a specified distance. The contacts are then stored
        as a list of lists and returned for further processing."""
        timeseries = []
        for ts in u.trajectory:
            dist = contacts.distance_array(group_a.positions, group_b.positions)
            n_contacts = contacts.contact_matrix(dist, radius).sum()
            timeseries.append([ts.frame, n_contacts])
        return np.array(timeseries)

    def saltbridges(self, top_file, traj_file):
        """Saltbridges or ionic bonds occur between acidic and basic residues within the protein structure. Here the
        residues and the respective atoms involved are selected and passed to the contacts within cutoff function,
        and returned as a dataframe which is then stored within the original system folder as a .csv file
        for figure generation."""
        u = mda.Universe(top_file, traj_file)
        sel_basic = "(resname ARG LYS) and (name NH* NZ)"
        sel_acidic = "(resname ASP GLU) and (name OE* OD*)"
        basic = u.select_atoms(sel_basic)
        acidic = u.select_atoms(sel_acidic)
        ca = self.contacts_within_cutoff(u, basic, acidic, radius=4.5)
        ca_df = pd.DataFrame(ca, columns=['Time (ns)', str(top_file).split('.')[0]])
        ca_df['Time (ns)'] = ca_df['Time (ns)']/10
        ca_df.to_csv(str(top_file.split('.')[0] + '_salt_bridges.csv'))
        return ca_df

    def saltbridges_comp(self, traj_folder, output_dir):
        """This fucntion iterates through each of the systems that were simulated and calls the saltbridges
        function to determine the respective contacts. Each of the salt bridges for each of the system under
        investigation have their .csv retreived and appended to a new dataframe for ease of comparison between the
        systems, which is then stored in the output specified folder."""
        df = pd.DataFrame()
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            sb_df = self.saltbridges(top_file, traj_file)
            if df.empty:
                df = pd.concat([df, sb_df], axis=1)
            else:
                df = pd.merge(df, sb_df, on = 'Time (ns)', how='inner')
        os.chdir(output_dir)
        df.to_csv('multi_variant_saltbridges.csv')

    def multi_sb_plotting(self, outputs_dir, start_fr):
        """This function plots the multi-variant.csv file generated by the saltbridges_comp function from within the
        previously defined output directory and creates a figure with each dataset present within the file. This
        function was developed prior to the single plot generation, however too many datasets resulted in degraded
        visualization."""
        os.chdir(outputs_dir)
        for file in os.listdir(outputs_dir):
            if file.endswith('multi_variant_saltbridges.csv'):
                sb_df = pd.read_csv(file)
                sb_df['Time (ns)'] = sb_df['Time (ns)']
                sb_df.iloc[int(start_fr):, 1:len(sb_df)].plot(x='Time (ns)', linewidth=0.5)
                plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
                plt.title('Multi-system Protein-Protein Ionic Interactions')
                plt.ylabel(r"$N_{II}$")
                plt.savefig('{0}/{1}{2}'.format(str(outputs_dir), str(file).split('.')[0], '.tiff'),
                            bbox_inches='tight', dpi=600)


    def hbond_calc(self, top_file, traj_file, start_fr):
        """Currently the set of ions and solvent ions have been predefined for this function and will need to be
        expanded upon for further diversity in molecular analyses. Here the script iterates through each of the
        segments within a universe created using the trajectory and topology files. The script then looks for segments
        that are not metalic ions, solvent ions or other protein segments and calculates the h-bonds present between the
        entire protein atom selection as well as all the atoms present in the filtered segments. This results in a
        function which is currently abe to identify the ligand and Nucleic Acid present within the system and execute the
        respective anallysis and return the stored results as a dataframe for each of the systems."""
        solvent_mol = ['POT', 'TIP3', "CLA"]
        ions = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
         'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'ZN2','GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y',
         'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA',
         'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE',
         'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP',
         'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG',
         'CN', 'NH', 'FL', 'MC', 'LV', 'TS', 'OG']

        mols = []
        u = mda.Universe(traj_file, top_file)
        hbonds_df = pd.DataFrame()
        df = pd.DataFrame()
        for seg in u.segments:
            if 'PRO' not in str(seg) and str(seg)[1:-1].split(' ')[1].split('_')[2] not in solvent_mol and \
                    str(seg)[1:-1].split(' ')[1].split('_')[2] not in ions and 'DNA' not in str(seg) and 'RNA' not in str(seg):
                hbonds = HydrogenBondAnalysis(universe=u, d_a_cutoff=4.5,  update_selections=False)
                hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
                hbonds.acceptors_sel = hbonds.guess_acceptors("segid " + str(seg)[1:-1].split(' ')[1])
                hbonds.run(verbose=True, start=int(start_fr))
                df1 = pd.DataFrame(hbonds.times, columns=['Time (ns)'])
                df2 = pd.DataFrame(hbonds.count_by_time(), columns=[str(top_file).split('.')[0]])
                df = pd.concat([df1, df2], axis=1)
                df.to_csv(str(seg)[1:-1].split(' ')[1].split('_')[2] + '_hbonds.csv')
                hbonds_df = pd.concat([hbonds_df, df], axis=1)

    def nucleic_prot_hbonds(self, top_file, traj_file, start_fr):
        u = mda.Universe(traj_file, top_file)
        hbonds_df = pd.DataFrame()
        df = pd.DataFrame()
        hbonds = HydrogenBondAnalysis(universe=u, d_a_cutoff=3.5, update_selections=True)
        hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
        hbonds.acceptors_sel = hbonds.guess_acceptors("nucleic")
        hbonds.run(verbose=True, start=int(start_fr))
        df1 = pd.DataFrame(hbonds.times, columns=['Time (ns)'])
        df2 = pd.DataFrame(hbonds.count_by_time(), columns= [str(top_file).split('.')[0] +  str(' Nucleic Acid ')])
        nucleic_df = pd.concat([df1, df2], axis=1)
        nucleic_df.to_csv('nucleic_hbonds.csv')

    def hbond_comp(self, traj_folder, output_dir, start_fr):
        """This fucntion iterates through each of the systems that were simulated and calls the hbonds
        function to determine the respective contacts. Each of the hbonds for each of the systems under
        investigation have their hbonds.csv retrieved and each of the columns at the respective positions
        (each column correlates to the same molecule investigated) and appended to a new dataframe for ease of
        comparison between the systems, which is then stored in the output specified folder."""
        df = pd.DataFrame()
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            # sys_hbonds = self.hbond_calc(traj_file, top_file, start_fr)
            sys_nucleic_hbonds = self.nucleic_prot_hbonds(traj_file, top_file, start_fr)


    def multi_hbonds_plots(self, outputs_dir, start_df):
        """This function plots the 'multi-variant_hbonds'.csv file generated by the hbonds_comp function from within the
        previously defined output directory and creates a figure with each dataset present within the file. This
        function was developed prior to the single plot generation, however too many datasets resulted in degraded
        visualization."""
        solvent_mol = ['POT', 'TIP3', "CLA"]
        ions = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
         'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'ZN2','GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y',
         'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA',
         'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE',
         'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP',
         'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG',
         'CN', 'NH', 'FL', 'MC', 'LV', 'TS', 'OG']
        sys_hetatms = []
        hetatm_file_pair = {}
        for system in os.listdir(systems):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(systems), system))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            u = mda.Universe(top_file, traj_file)
            for seg in u.segments:
                if 'PRO' not in str(seg) and str(seg)[1:-1].split(' ')[1].split('_')[2] not in solvent_mol and str(seg)[1:-1].split(' ')[1].split('_')[2] not in ions and 'DNA' not in str(seg) and 'RNA' not in str(seg):
                    hetatm = str(seg).split(' ')[1].split('_')[-1][:-1]
                    if hetatm not in sys_hetatms:
                        sys_hetatms.append(hetatm)
        for hetatm in sys_hetatms:
            for system in os.listdir(systems):
                os.chdir('{0}/{1}'.format(str(systems), system))
                for sys_file in os.listdir():
                    if str(hetatm) in str(sys_file) and 'hbonds.csv' in str(sys_file):
                         hetatm_file_pair.setdefault(hetatm, []).append('{0}/{1}'.format(str(os.getcwd()), sys_file))
        for key in hetatm_file_pair:
            df = pd.DataFrame()
            for hetatm_data in hetatm_file_pair[key]:
                hetatm_df = pd.read_csv(hetatm_data)
                hetatm_df['Time (ns)'] =  hetatm_df['Time (ns)']/1000
                if df.empty == True:
                    df = hetatm_df.iloc[:,1:]
                else:
                    df = pd.concat([df, hetatm_df.iloc[:, -1]], axis=1)
            df.plot(x='Time (ns)', linewidth=0.5)
            plt.title('Multi-system Protein-'  + str(key) + ' H-bonds')
            plt.ylabel(r"$N_{HB}$")
            plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            plt.savefig('{0}/{1}'.format(str(output_dir), 'Multi-system Protein-'  + str(key) + '_hbonds.tiff'), bbox_inches='tight', dpi=900)

    def multi_nucleic_hbonds_plots(self, outputs_dir, start_df):
        """This function plots the 'multi-variant_hbonds'.csv file generated by the hbonds_comp function from within the
        previously defined output directory and creates a figure with each dataset present within the file. This
        function was developed prior to the single plot generation, however too many datasets resulted in degraded
        visualization."""
        solvent_mol = ['POT', 'TIP3', "CLA"]
        ions = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
         'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'ZN2','GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y',
         'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA',
         'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE',
         'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP',
         'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG',
         'CN', 'NH', 'FL', 'MC', 'LV', 'TS', 'OG']
        sys_hetatms = []
        hetatm_file_pair = {}
        df = pd.DataFrame()
        for system in os.listdir(systems):
            os.chdir('{0}/{1}'.format(str(systems), system))
            for file in os.listdir():
                if str('nucleic') in str(file) and 'hbonds.csv' in str(file):
                    nucleic_hbonds_df = pd.read_csv(file)
                    nucleic_hbonds_df['Time (ns)'] = nucleic_hbonds_df['Time (ns)']/1000
                    if df.empty == True:
                        df = nucleic_hbonds_df.iloc[:,1:]
                    else:
                        df = pd.concat([df, nucleic_hbonds_df.iloc[:,2:]], axis=1)
            os.chdir(systems)
        df.plot(x='Time (ns)', linewidth=0.5)
        plt.title('Multi-system Protein-Nucelic Acid H-bonds')
        plt.ylabel(r"$N_{HB}$")
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.savefig('{0}/{1}'.format(str(output_dir), 'Multi-system Protein-Nucleic_Acid_hbonds.tiff'), bbox_inches='tight', dpi=900)

    def hetatm_calc(self, top_file, traj_file, start_fr):
        print(datetime.datetime.now())
        mols = []
        solvent_mol = ['POT', 'TIP3', 'CLA']
        hbonds_df = pd.DataFrame()
        u = mda.Universe(top_file, traj_file)
        for seg in u.segments:
            if 'PRO' not in str(seg) and str(seg)[1:-1].split(' ')[1].split('_')[2] not in solvent_mol and 'DNA' not in str(seg) and 'RNA' not in str(seg):
                mols.append(str(seg).split(' ')[1][:-1])
        for mol in mols:
            for mol2 in mols:
                if mol != mol2:
                    hbonds = HydrogenBondAnalysis(universe=u, d_a_cutoff=4.5, update_selections=True)
                    hbonds.hydrogens_sel = hbonds.guess_hydrogens('segid ' + str(mol))
                    hbonds.acceptors_sel = hbonds.guess_acceptors('segid ' + str(mol2))
                    hbonds.run(verbose=True, start=int(start_fr))
                    df1 = pd.DataFrame(hbonds.times, columns=['Time (ns)'])
                    df2 = pd.DataFrame(hbonds.count_by_time(), columns=[str(top_file).split('.')[0]])
                    df = pd.concat([df1, df2], axis=1)
                    df.to_csv(str(mol) + str(mol2) + '_hetatm_hbonds.csv')
                    hbonds_df = pd.concat([hbonds_df, df], axis=1)
            mols.remove(mol)

    def hetatm_conts(self, traj_folder, start_fr):
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            self.hetatm_calc(top_file, traj_file, start_fr)

    def hetatm_plot(self, systems, outputs_dir, start_fr):
        os.chdir(systems)
        system_hetatm_conts = {}
        for folder in os.listdir(systems):
            os.chdir(folder)
            for file in os.listdir():
                if file.endswith('hetatm_contacts.csv'):
                    system_hetatm_conts.setdefault(folder, []).append('{0}/{1}'.format(str(os.getcwd()), file))
            os.chdir(systems)
        for key in system_hetatm_conts:
            contact_sys = len(system_hetatm_conts[key])
            for i in range(contact_sys):
                df = pd.DataFrame()
                for key in system_hetatm_conts:
                    title_hetatm = str(os.path.basename(system_hetatm_conts[key][i])).split('_')[1:3]
                    hetatm_df = pd.read_csv(system_hetatm_conts[key][i])
                    if df.empty == True:
                        df = hetatm_df.iloc[8000:,1:]
                    else:
                        df = pd.concat([df, hetatm_df.iloc[8000:, -1]], axis=1)
                df.plot(x='Time (ns)', linewidth=0.5)
                plt.title('Multi-system ' + str(title_hetatm[0]) + '-'+ str(title_hetatm[1]) + ' interactions')
                plt.ylabel(r"$N_{II}$")
                plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
                plt.savefig('{0}/{1}{2}'.format(outputs_dir, 'Multi-system ' + str(title_hetatm[0])+str(title_hetatm[1]), ' interactions.tiff'), bbox_inches='tight', dpi=900)

    def nucleic_saltbridges_calc(self, top_file, traj_file):
        """Saltbridges or ionic bonds occur between acidic and basic residues within the protein structure. Here the
        residues and the respective atoms involved are selected and passed to the contacts within cutoff function,
        and returned as a dataframe which is then stored within the original system folder as a .csv file
        for figure generation."""
        u = mda.Universe(top_file, traj_file)
        sel_basic = "(resname ARG LYS) and (name NH* NZ)"
        sel_acidic = "nucleic and (name *P)"
        group1 = u.select_atoms(sel_basic)
        group2 = u.select_atoms(sel_acidic)
        ca = self.contacts_within_cutoff(u, group1, group2, radius=4.5)
        ca_df = pd.DataFrame(ca, columns=['Time (ns)', str(top_file).split('.')[0]])
        ca_df['Time (ns)'] = ca_df['Time (ns)'] /10
        ca_df.to_csv(str(top_file.split('.')[0] + '_Nucleic_Acid_saltbridges.csv'))

    def nucleic_ionic_conts(self, traj_folder):
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}/{1}'.format(str(traj_folder), folder))
            for file in os.listdir():
                if file.endswith('.xtc') and 'rmsfit_' not in str(file):
                    traj_file = file
                elif file.endswith('.tpr'):
                    top_file = file
            self.nucleic_saltbridges_calc(top_file, traj_file)

    def ionic_nucleic_plot(self, systems, outputs_dir, start_fr):
        os.chdir(systems)
        multi_df = pd.DataFrame()
        for folder in os.listdir(systems):
            os.chdir(folder)
            for file in os.listdir():
                if file.endswith('Nucleic_Acid_saltbridges.csv'):
                    heading = str(file).split('.')[0].split('_')
                    title = heading[0] + ' ' + heading[1] + '-' + heading[2]
                    ionic_nucleic_df = pd.read_csv(file)
                    ionic_nucleic_df = ionic_nucleic_df.iloc[int(start_fr):, 1:]
                    if multi_df.empty == True:
                        multi_df = ionic_nucleic_df
                    else:
                        multi_df = pd.concat([multi_df, ionic_nucleic_df.iloc[:,1]], axis=1)
            os.chdir(systems)
        multi_df.plot(x='Time (ns)', linewidth=0.5)
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.ylabel(r"$N_{II}$")
        plt.title('Multi-system Protein-Nucleic Ionic Interactions')
        plt.savefig('{0}/{1}{2}'.format(outputs_dir, 'Multi-system Protein-Nucleic Ionic Interactions', '.tiff'), bbox_inches='tight', dpi=900)

    def traj_mean(self, systems):
        os.chdir(systems)
        for system in os.listdir():
            print(system)
            os.chdir(system)
            stats_file = open(f'{system}_mean_std_stats.txt', 'w')
            for file in os.listdir():
                if file.endswith('.csv'):
                    if 'pca' not in str(file).lower():
                        df = pd.read_csv(file)
                        df = df.iloc[:,1:]
                        columns = list(df)
                        for i in columns:
                            if i != "Time (ns)":
                                if i !="Frame":
                                    stats_file.write(file)
                                    stats_file.write('\n')
                                    stats_file.write(f"Mean {i}: {df[i].mean()/10}")
                                    stats_file.write('\n')
                                    stats_file.write(f"Standard deviation: {df[i].std()/10}")
                                    stats_file.write('\n')
                                    stats_file.write('\n')
            stats_file.close()
            os.chdir(systems)

def main(systems, output_dir, start_fr):
    p = Contacts()
    p.saltbridges_comp(systems, output_dir)
    p.multi_sb_plotting(output_dir, start_fr)
    p.hbond_comp(systems, output_dir, start_fr)
    p.multi_hbonds_plots(output_dir, start_fr)
    p.multi_nucleic_hbonds_plots(output_dir, start_fr)
    p.hetatm_conts(systems, start_fr)
    p.hetatm_plot(systems, output_dir, start_fr)
    p.nucleic_ionic_conts(systems)
    p.ionic_nucleic_plot(systems, output_dir, start_fr)

if __name__ =='__main__':
    p = TrajStat()
    parser = argparse.ArgumentParser()
    parser.add_argument("--systems", help="Path to the trajectory files have been simulated. ")
    parser.add_argument("--output_dir", help="Path to the directory that the variant systems will be stored in")
    args = parser.parse_args()
    systems = str(args.systems)
    output_dir = str(args.output_dir)
    pool = multiprocessing.Pool()
    pool2 = multiprocessing.Pool()
    pool4 = multiprocessing.Pool()
    os.chdir(systems)
    traj_files = []
    top_files = []
    for system in os.listdir(systems):
        os.chdir(system)
        for file in os.listdir():
            if file.endswith('.xtc') and 'rmsfit' not in str(file):
                traj_files.append('{0}/{1}'.format(os.getcwd(), file))
            if file.endswith('.tpr'):
                top_files.append('{0}/{1}'.format(os.getcwd(), file))
        os.chdir(systems)
    traj_files.sort()
    top_files.sort()
    for i,j in zip(top_files, traj_files):
        pool2.apply_async(p.rmsd_calc, args=(i,j,))
    pool2.close()
    pool2.join()
    p.rmsd_plot(systems, output_dir)
    for i,j in zip(top_files, traj_files):
        pool4.apply_async(p.rmsd_na_calc, args=(i,j,))
    pool4.close()
    pool4.join()
    p.rmsd_na_plot(systems, output_dir)
    start_fr = input("Enter starting frame of trajectory equilibration:")
    for i,j in zip(top_files, traj_files):
        pool.apply_async(p.rmsf_calc, args=(i,j,int(start_fr),))
    pool.close()
    pool.join()
    p.rmsf_plot(systems, output_dir)
    pool3 = multiprocessing.Pool()
    for i,j in zip(top_files, traj_files):
            pool3.apply_async(p.rgyr_calc, args=(i,j,int(start_fr),))
    pool3.close()
    pool3.join()
    p.rgyr_plot(systems, output_dir, start_fr)
    p.pca_comp(systems, systems, start_fr)
    p.pca_plotting(systems, output_dir, start_fr)
    p.pca_time_cluster_plotting(systems, output_dir, start_fr)
    p.saltbridges_comp(systems, output_dir)
    p.multi_sb_plotting(output_dir, start_fr)
    p.hbond_comp(systems, output_dir, start_fr)
    p.multi_hbonds_plots(output_dir, start_fr)
    p.multi_nucleic_hbonds_plots(output_dir, start_fr)
    p.hetatm_conts(systems, start_fr)
    p.hetatm_plot(systems, output_dir, start_fr)
    p.nucleic_ionic_conts(systems)
    p.ionic_nucleic_plot(systems, output_dir, start_fr)
    main(systems, output_dir, start_fr)
    p.traj_mean(systems)

