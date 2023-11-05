import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.tests.datafiles import PDB, GRO, XTC, TPR
from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction, analysis_class
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
        R = rms.RMSD(u, u, select='backbone', ref_frame=0).run(verbose=True)
        dataframe = pd.DataFrame(R.results.rmsd, columns=['Frame','Time (ns)', str(os.path.basename(top_file).split('.')[0])])
        dataframe['Time (ns)'] = dataframe['Time (ns)']/1000
        dataframe.to_csv(top_file.split('.')[0] + '_rmsd_2.csv')

    def rmsf_calc(self, top_file, traj_file, start_fr):
        os.chdir(os.path.dirname(top_file))
        u = mda.Universe(top_file, traj_file)
        average = align.AverageStructure(u, u, select='protein', ref_frame=0).run(verbose=True, start=int(start_fr))
        ref = average.results.universe
        aligner = align.AlignTraj(u, ref, select='protein', in_memory=False).run(verbose=True, start=int(start_fr))
        # Select the protein atoms
        u = mda.Universe(top_file, '{0}_{1}'.format('rmsfit', os.path.basename(traj_file)))
        # Calculate the RMSF
        for seg in u.segments:
            if 'PRO' in str(seg):
                protein = u.select_atoms('segid ' + str(str(seg).split(' ')[1][:-1]) + ' and name CA')
                rmsf = rms.RMSF(protein).run(verbose=True)
                rmsf_df = pd.DataFrame(rmsf.results.rmsf)
                rmsf_df.to_csv('{3}\{1}_{2}'.format(top_file.split('.')[0], str(str(seg).split(' ')[1][:-1]), 'rmsf.csv', os.path.dirname(top_file)))

    def rmsd_plot(self, systems, output_dir):
        multi_df = pd.DataFrame()
        for system in os.listdir(systems):
            os.chdir(system)
            for file in os.listdir():
                if file.endswith('rmsd_2.csv'):
                    df = pd.read_csv(file)
                    if multi_df.empty ==True:
                        multi_df = pd.concat([multi_df, df.iloc[:,2:]], axis=1)
                    else:
                        multi_df = pd.merge(multi_df, df.iloc[:,2:], on='Time (ns)', how='inner')
            os.chdir(systems)
        multi_df.plot(x='Time (ns)', linewidth=0.75)
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.ylabel('RMSD ($\AA$)')
        plt.savefig('{0}\{1}{2}'.format(str(output_dir), 'multi_rmsd_2', '.tiff'),
        bbox_inches='tight', dpi=900)


    def rmsf_plot(self, systems, output_dir):
        prot_chains = []
        chain_file_pair = {}
        for system in os.listdir(systems):
            top_file = ''
            traj_file = ''
            os.chdir('{0}\{1}'.format(str(systems), system))
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
                os.chdir('{0}\{1}'.format(str(systems), system))
                for rmsf_file in os.listdir():
                    if str(rmsf_chain) in str(rmsf_file):
                        chain_file_pair.setdefault(rmsf_chain, []).append('{0}\{1}'.format(str(os.getcwd()), rmsf_file))
        print(chain_file_pair)
        for key in chain_file_pair:
            df = pd.DataFrame()
            for rmsf_data in chain_file_pair[key]:
                rmsf_df = pd.read_csv(rmsf_data)
                rmsf_df.rename(columns={'Unnamed: 0': 'Residue Number'}, inplace=True)
                rmsf_df.rename(columns={'0': str(os.path.dirname(rmsf_data).split('\\')[-1])}, inplace=True)
                if df.empty == True:
                    df = rmsf_df
                else:
                    df = pd.concat([df, rmsf_df.iloc[:,1]], axis=1)
            df.plot(x='Residue Number', linewidth=0.75)
            plt.title('Multi-system RMSF of ' + str(key) + ' chain')
            plt.ylabel('RMSD ($\AA$)')
            plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            plt.savefig('{0}\{1}_{2}'.format(str(output_dir),  str(key),  '_rmsf.tiff'), bbox_inches='tight', dpi=900)


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
                                   total_mass=np.sum(protein.masses)).run(verbose=True, start=int(start_fr))
        rgyr_df = pd.DataFrame(rog.results['timeseries'],
                               columns=[str(os.path.basename(top_file).split('.')[0]), 'x-axis', 'y-axis', 'z-axis'])
        rgyr_df.to_csv(str(top_file).split('.')[0] + '_rgyr.csv')

    def rgyr_plot(self, systems, output_dir, start_fr):
        multi_df = pd.DataFrame()
        for system in os.listdir(systems):
            os.chdir(system)
            for file in os.listdir():
                if file.endswith('rgyr.csv'):
                    df = pd.read_csv(file)
                    df.iloc[:, 0] = (df.iloc[:, 0] + int(start_fr)) / 100
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
            plt.savefig('{0}\{1}{2}'.format(str(output_dir), 'multi_rgyr', '.tiff'),
                            bbox_inches='tight', dpi=900)

    def pca_calc(self, top_file, traj_file, start_fr):
        cumulated_PCA = {}
        print(datetime.datetime.now())
        u = mda.Universe(top_file, traj_file)
        backbone = u.select_atoms('backbone and not name O')
        pc = pca.PCA(u, select='backbone and not name O', align=True, mean=None, n_components=None, verbose=True)
        pc.run(verbose=True, start=int(start_fr), step=100)
        transformed = pc.transform(backbone)
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
            os.chdir('{0}\{1}'.format(str(traj_folder), folder))
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
            os.chdir('{0}\{1}'.format(str(traj_folder), folder))
            broad_df = pd.DataFrame()
            for file in os.listdir():
                if file.endswith('_PCA_MDA.csv'):
                    df = pd.read_csv(file)
                    df = df.drop(columns=['PC3'])
                    df.rename(columns={'Time (ps)': 'Time (ns)'}, inplace=True)
                    df['Time (ns)'] = df['Time (ns)'] / 1000
                    df = df.iloc[int(start_fr):, 1:]
                    ax1.scatter(df['PC1'], df['PC2'], s=1, marker="o", label=str(folder))
                    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
                    plt.title('Multi-system PCA')
                    plt.xlabel('PC1 ($\AA$)')
                    plt.ylabel('PC2 ($\AA$)')
        plt.savefig('{0}\{1}{2}'.format(str(output_dir), 'multi_PCA', '.tiff'),
                    bbox_inches='tight', dpi=900)

    def pca_time_cluster_plotting(self, traj_folder, output_dir, start_fr):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for folder in os.listdir(traj_folder):
            top_file = ''
            traj_file = ''
            os.chdir('{0}\{1}'.format(str(traj_folder), folder))
            broad_df = pd.DataFrame()
            for file in os.listdir():
                if file.endswith('_PCA_MDA.csv'):
                    df = pd.read_csv(file)
                    df = df.drop(columns=['PC3'])
                    df.rename(columns={'Time (ps)': 'Time (ns)'}, inplace=True)
                    df['Time (ns)'] = df['Time (ns)'] / 1000
                    df = df.iloc[int(start_fr):, 1:]
                    fig = px.scatter(df, x='PC1', y='PC2', color=df["Time (ns)"])
                    fig.write_image('{0}\{1}{2}'.format(str(output_dir), str(folder), '_PCA_timeclustering.png'),
                                    width=1 * 900,
                                    height=1 * 900, scale=1)


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
    os.chdir(systems)
    traj_files = []
    top_files = []
    for system in os.listdir(systems):
        os.chdir(system)
        for file in os.listdir():
            if file.endswith('.xtc') and 'rmsfit' not in str(file):
                traj_files.append('{0}\{1}'.format(os.getcwd(), file))
            if file.endswith('.tpr'):
                top_files.append('{0}\{1}'.format(os.getcwd(), file))
        os.chdir(systems)
    traj_files.sort()
    top_files.sort()
    for i,j in zip(top_files, traj_files):
        pool2.apply_async(p.rmsd_calc, args=(i,j,))
    pool2.close()
    pool2.join()
    p.rmsd_plot(systems, output_dir)
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




