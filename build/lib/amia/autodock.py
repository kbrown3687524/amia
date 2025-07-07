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

import os
import subprocess
from multiprocessing import Pool, cpu_count
import re
import plotly.graph_objs as go
import plotly.offline as pyo
import argparse
import datetime


class DockingPipeline:
    def convert_wt_pdb_to_pdbqt(self, input_pdb_path, output_dir):
        """
        Convert WT PDB to PDBQT and store in the specified variant_dir (output_dir).
        """
        if not os.path.isdir(output_dir):
            print(f"[ERROR] Output directory does not exist: {output_dir}")
            return None

        base_name = os.path.basename(input_pdb_path).replace('_auto.pdb', '')  # Extract 'WT' from 'WT_auto.pdb'
        output_pdbqt_path = os.path.join(output_dir, f"{base_name}_target.pdbqt")

        command = [
            "obabel", input_pdb_path,
            "-O", output_pdbqt_path,
            "-xr", "-p", "7.4", "--partialcharge", "eem"
        ]

        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"[INFO] Converted {os.path.basename(input_pdb_path)} -> {os.path.basename(output_pdbqt_path)} in {output_dir}")
            return output_pdbqt_path
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed to convert {input_pdb_path}:\n{e.stderr.decode()}")
            return None
    
    def _convert_wrapper(self, pdb_file_path):
        """
        Helper wrapper for multiprocessing. Converts a single *_auto.pdb file to *_target.pdbqt
        using the same directory as output.
        """
        output_dir = os.path.dirname(pdb_file_path)
        self.convert_wt_pdb_to_pdbqt(pdb_file_path, output_dir)

    def batch_convert_directory(self, directory_path, num_processes=None):
        if not os.path.isdir(directory_path):
            print(f"[ERROR] Provided path is not a directory: {directory_path}")
            return

        pdb_files = [
            os.path.join(directory_path, f)
            for f in os.listdir(directory_path)
            if f.endswith("_auto.pdb")
        ]

        if not pdb_files:
            print("[INFO] No *_auto.pdb files found.")
            return

        num_processes = num_processes or min(cpu_count(), len(pdb_files))
        print(f"[INFO] Converting {len(pdb_files)} files using {num_processes} processes...")

        with Pool(processes=num_processes) as pool:
            pool.map(self._convert_wrapper, pdb_files)

    def generate_conf_files(self, target_dir, center, size=(30.0, 30.0, 30.0)):
        """
        Generate *_conf.txt files for each *_target.pdbqt file in the directory.
        Includes both wild-type and variant receptors.
        """
        if not os.path.isdir(target_dir):
            print(f"[ERROR] Provided path is not a directory: {target_dir}")
            return

        target_files = sorted([
            f for f in os.listdir(target_dir)
            if f.endswith("_target.pdbqt")
        ])

        if not target_files:
            print("[INFO] No *_target.pdbqt files found for configuration.")
            return

        center_x, center_y, center_z = center
        size_x, size_y, size_z = size

        for target_file in target_files:
            variant_base = target_file.replace("_target.pdbqt", "")
            conf_filename = f"{variant_base}_conf.txt"
            conf_path = os.path.join(target_dir, conf_filename)

            with open(conf_path, 'w') as f:
                f.write(f"center_x = {center_x:.3f}\n")
                f.write(f"center_y = {center_y:.3f}\n")
                f.write(f"center_z = {center_z:.3f}\n")
                f.write("\n")  # spacer between center and size
                f.write(f"size_x = {size_x:.3f}\n")
                f.write(f"size_y = {size_y:.3f}\n")
                f.write(f"size_z = {size_z:.3f}\n")

            print(f"[INFO] Created config: {conf_filename}")

    def smiles_to_3d_structure(self, smiles, output_path, output_format="pdbqt"):
        """
        Convert a SMILES string into a 3D structure file using Open Babel.

        Args:
            smiles (str): The SMILES string.
            output_path (str): Full path to the output file (should match format).
            output_format (str): Output format (e.g., "pdbqt", "pdb", "mol2").
        
        Returns:
            bool: True if successful, False otherwise.
        """
        if not smiles.strip():
            print("[ERROR] Empty SMILES string provided.")
            return False

        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        cmd = [
            "obabel",
            f"-:{smiles}",
            "--gen3d",
            f"-o{output_format}",
            f"-O{output_path}"
        ]

        try:
            print(f"[INFO] Generating 3D structure at: {output_path}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"[SUCCESS] 3D structure written to {output_path}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Open Babel conversion failed.")
            print(e.stdout)
            print(e.stderr)
            return False

    def run_vina_docking(self, target_dir, ligand_path, vina_path="vina"):
        receptor_files = [
            f for f in os.listdir(target_dir)
            if f.endswith("_target.pdbqt")
        ]

        if not receptor_files:
            print("[INFO] No *_target.pdbqt files found.")
            return

        for receptor_file in receptor_files:
            variant_base = receptor_file.replace("_target.pdbqt", "")
            receptor_path = os.path.join(target_dir, receptor_file)
            config_path = os.path.join(target_dir, f"{variant_base}_conf.txt")
            output_path = os.path.join(target_dir, f"{variant_base}_docked.pdbqt")
            log_path = os.path.join(target_dir, f"{variant_base}_log.txt")

            if not os.path.exists(config_path):
                print(f"[WARNING] Missing config: {config_path}")
                continue

            cmd = [
                vina_path,
                "--receptor", receptor_path,
                "--ligand", ligand_path,
                "--config", config_path,
                "--out", output_path,
                "--log", log_path
            ]

            try:
                print(f"[INFO] Docking to {receptor_file}...")
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"[SUCCESS] Docked: {variant_base}_docked.pdbqt (log: {variant_base}_log.txt)")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] Docking failed for {variant_base}")
                print(e.stdout)
                print(e.stderr)
    
    def save_visualizations_with_pymol(self, target_dir):
        """
        For each receptor-ligand pair, load them in PyMOL, save .pse session and
        export one combined .pdb file with receptor and ligand together.
        """
        from pymol import cmd

        receptor_files = [
            f for f in os.listdir(target_dir)
            if f.endswith("_target.pdbqt")
        ]

        if not receptor_files:
            print("[INFO] No *_target.pdbqt files found.")
            return

        for receptor_file in receptor_files:
            variant_base = receptor_file.replace("_target.pdbqt", "")
            receptor_pdb = os.path.join(target_dir, receptor_file.replace("_target.pdbqt", "_auto.pdb"))
            ligand_pdbqt = os.path.join(target_dir, f"{variant_base}_docked.pdbqt")
            session_file = os.path.join(target_dir, f"{variant_base}.pse")
            combined_pdb = os.path.join(target_dir, f"{variant_base}_complex.pdb")

            if not (os.path.exists(receptor_pdb) and os.path.exists(ligand_pdbqt)):
                print(f"[WARNING] Missing receptor or ligand for {variant_base}")
                continue

            try:
                cmd.reinitialize()
                cmd.load(receptor_pdb, "receptor")
                cmd.load(ligand_pdbqt, "ligand")

                # Save PyMOL session
                cmd.save(session_file)
                print(f"[INFO] Saved PyMOL session: {session_file}")

                # Save combined receptor + ligand PDB
                cmd.create("complex", "receptor or ligand")  # Merge selections
                cmd.save(combined_pdb, "complex")
                print(f"[INFO] Saved combined PDB: {combined_pdb}")

            except Exception as e:
                print(f"[ERROR] PyMOL visualization failed for {variant_base}: {e}")


    def parse_vina_logs(self, target_dir, html_output="vina_results.html", plot_output="vina_affinity_plot.html"):
        """
        Parse *_log.txt files, tabulate all modes into HTML, and plot only mode 1 affinities.

        Returns:
            results_dict: {variant_base: [ {mode, affinity, rmsd_l, rmsd_u}, ... ] }
            Also writes:
            - HTML file with table of full results
            - HTML file with bar plot of mode 1 affinities
        """
        logs = [f for f in os.listdir(target_dir) if f.endswith("_log.txt")]
        results = {}

        pattern = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$")

        for log_file in logs:
            variant_base = log_file.replace("_log.txt", "")
            log_path = os.path.join(target_dir, log_file)

            modes = []

            with open(log_path, "r") as f:
                lines = f.readlines()

            capture = False
            for line in lines:
                line = line.strip()
                if line.startswith("mode |"):
                    capture = True
                    continue
                if capture:
                    if line == "" or line.startswith("Writing output"):
                        break
                    match = pattern.match(line)
                    if match:
                        modes.append({
                            "mode": int(match.group(1)),
                            "affinity": float(match.group(2)),
                            "rmsd_l": float(match.group(3)),
                            "rmsd_u": float(match.group(4))
                        })

            results[variant_base] = modes

        # --- HTML TABLE for all modes ---
        html = "<html><head><title>Vina Docking Results</title></head><body>"
        html += "<h2>Docking Results Summary</h2>"

        for variant, modes in results.items():
            html += f"<h3>{variant}</h3>"
            html += "<table border='1' cellpadding='5' cellspacing='0'>"
            html += "<tr><th>Mode</th><th>Affinity (kcal/mol)</th><th>RMSD l.b.</th><th>RMSD u.b.</th></tr>"
            for mode in modes:
                html += (
                    f"<tr><td>{mode['mode']}</td>"
                    f"<td>{mode['affinity']:.2f}</td>"
                    f"<td>{mode['rmsd_l']:.2f}</td>"
                    f"<td>{mode['rmsd_u']:.2f}</td></tr>"
                )
            html += "</table><br>"

        html += "</body></html>"

        html_path = os.path.join(target_dir, html_output)
        with open(html_path, "w") as f:
            f.write(html)
        print(f"[INFO] Full results HTML saved to: {html_path}")

        # --- Plot only Mode 1 ---
        variants = []
        affinities = []

        for variant, modes in results.items():
            mode1 = next((m for m in modes if m["mode"] == 1), None)
            if mode1:
                variants.append(variant)
                affinities.append(mode1["affinity"])

        trace = go.Bar(
            x=variants,
            y=affinities,
            marker=dict(color='steelblue'),
            name='Mode 1 Affinity'
        )

        layout = go.Layout(
            title='Mode 1 Affinity Scores',
            xaxis=dict(title='Variant'),
            yaxis=dict(title='Affinity (kcal/mol)', autorange='reversed'),
        )

        fig = go.Figure(data=[trace], layout=layout)
        plot_path = os.path.join(target_dir, plot_output)
        pyo.plot(fig, filename=plot_path, auto_open=False)
        print(f"[INFO] Plot saved to: {plot_path}")

        return results



# === Main Execution ===

def main():
    parser = argparse.ArgumentParser(description="Prepare files for AutoDock Vina docking.")
    parser.add_argument("--pdb_file", required=True, help="Path to the input PDB file (receptor)")
    parser.add_argument("--output_dir", required=True, help="Directory to store converted files and results")
    parser.add_argument("--smiles", required=False, default="CC1=C(C(=O)N2CCCCC2=N1)CCN3CCC(CC3)C4=NOC5=C4C=CC(=C5)F",
                        help="SMILES string of ligand (optional, default = Dolutegravir)")
    parser.add_argument("--center", nargs=3, type=float, required=True,
                        help="Center coordinates for docking box (x y z)")

    args = parser.parse_args()
    pdb_file = args.pdb_file
    output_dir = args.output_dir
    smiles = args.smiles
    center_coords = tuple(args.center)

    print("📄 Input PDB:", pdb_file)
    print("📁 Output Dir:", output_dir)

    pipeline = DockingPipeline()

    # Step 1: Convert PDB to PDBQT
    print("🔧 Converting receptor to PDBQT...")
    receptor_pdbqt = pipeline.convert_wt_pdb_to_pdbqt(pdb_file, output_dir)
    if not receptor_pdbqt:
        print("❌ Conversion failed.")
        return

    #Step 2: Batch convert directory if needed
    pipeline.batch_convert_directory(output_dir, num_processes=12)

    # Step 3: Generate config file
    print("🛠️ Generating docking config file...")
    pipeline.generate_conf_files(output_dir, center=center_coords)

    # Step 4: Generate ligand from SMILES
    ligand_path = os.path.join(output_dir, "ligand_from_smiles.pdbqt")
    print("🧪 Generating ligand structure from SMILES...")
    ligand_success = pipeline.smiles_to_3d_structure(smiles, ligand_path)
    if not ligand_success:
        print("❌ Ligand generation failed.")
        return

    # Step 5: Run docking
    print("🚀 Running docking...")
    pipeline.run_vina_docking(output_dir, ligand_path)

    # Step 6: Save PyMOL visualizations
    print("🧬 Generating visualizations with PyMOL...")
    pipeline.save_visualizations_with_pymol(output_dir)

    # Step 7: Parse logs and generate reports
    print("📊 Parsing docking logs and generating summary...")
    pipeline.parse_vina_logs(output_dir)

    print("✅ Finished at:", datetime.datetime.now())

if __name__ == "__main__":
    main()
