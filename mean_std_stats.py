import os
import pandas as pd
import argparse

class TrajOutput:
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

if __name__ == "__main__":
    p = TrajOutput()
    parser = argparse.ArgumentParser()
    parser.add_argument("--systems", help="Path to the trajectory files have been simulated. ")
    args = parser.parse_args()
    systems = str(args.systems)
    p.traj_mean(systems)
