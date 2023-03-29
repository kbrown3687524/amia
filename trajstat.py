import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.tests.datafiles import PDB, GRO, XTC, TPR
from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction, analysis_class
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np, os
import multiprocessing
import seaborn as sns
import datetime
import multiprocessing
from MDAnalysis.analysis import pca
from sklearn.decomposition import PCA

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
        for system in os.listdir(systems):
            os.chdir(system)
            for file in os.listdir():
                if file.endswith('rmsf.csv'):
                    df = pd.read_csv(file)
                    df.rename(columns={'Unnamed: 0': 'Residue Number'}, inplace=True)
                    df.plot(x='Residue Number', linewidth=0.75)
                    plt.legend([str(file).split('_rmsf.csv')[0]], bbox_to_anchor=(1.0, 1.0), loc='upper left')
                    plt.title(str(system))
                    plt.ylabel('RMSD ($\AA$)')
                    plt.savefig('{0}\{3}_{1}{2}'.format(str(output_dir), str(file).split('.')[0], '.tiff', str(system)),
                                bbox_inches='tight', dpi=900)
            os.chdir(systems)

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
            plt.savefig('{0}\{1}{2}'.format(str(output_dir), 'multi_rgyr_2', '.tiff'),
                            bbox_inches='tight', dpi=900)

    def pc_func(self, top_file, traj_file, output_dir):
        print(datetime.datetime.now())
        traj_file = '{2}\{0}_{1}'.format('rmsfit', os.path.basename(traj_file), os.path.dirname(traj_file))
        print(traj_file)
        u = mda.Universe(top_file, traj_file)
        sel = u.select_atoms('protein')
        pos = np.array([sel.positions for ts in u.trajectory])
        print(datetime.datetime.now())
        pca = PCA(n_components=2)
        pca.fit(pos.reshape(pos.shape[0], -1))
        trans_data = pca.transform(pos.reshape(pos.shape[0], -1))
        # Plot results
        plt.scatter(trans_data[:, 0], trans_data[:, 1])
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title(os.path.basename(top_file).split('.')[0])
        plt.savefig('{0}\{1}{2}'.format(str(output_dir), os.path.basename(top_file).split('.')[0], '_PCA2.tiff'), dpi=900)
        plt.clf()
        print(datetime.datetime.now())

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
    rmsd_plot(systems, output_dir)
    start_fr = input("Enter starting frame of trajectory equilibration:")
    for i,j in zip(top_files, traj_files):
        pool.apply_async(p.rmsf_calc, args=(i,j,int(start_fr),))
    pool.close()
    pool.join()
    rmsf_plot(systems, output_dir)
    pool3 = multiprocessing.Pool()
    for i,j in zip(top_files, traj_files):
            pool3.apply_async(p.rgyr_calc, args=(i,j,int(start_fr),))
    pool3.close()
    pool3.join()
    rgyr_plot(systems, output_dir, start_fr)
    for i,j in zip(top_files, traj_files):
        p.pc_func(i,j,output_dir)



