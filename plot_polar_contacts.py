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

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import argparse
import os

# Read data from CSV
def load_data(csv_path):
    df = pd.read_csv(csv_path)
    df = df.fillna("-")
    # Ensure Contact Difference is numeric
    df["Contact Difference"] = pd.to_numeric(df["Contact Difference"], errors="coerce").fillna(0).astype(int)
    return df

# Bar Plot of Contact Differences
def plot_contact_difference_bar(df, save=False, output_dir="."):
    summary = df.groupby("PDB File")["Contact Difference"].sum().reset_index()
    colors = summary["Contact Difference"].map(lambda x: 'blue' if x > 0 else 'red' if x < 0 else 'gray')

    plt.figure(figsize=(10, 6))
    plt.bar(summary["PDB File"], summary["Contact Difference"], color=colors)
    plt.axhline(0, color='black', linewidth=0.8)
    plt.title("Contact Differences by Variant (WT vs Mutant)")
    plt.xlabel("PDB Variant")
    plt.ylabel("Net Contact Difference")
    plt.xticks(rotation=45)
    plt.tight_layout()

    if save:
        path = os.path.join(output_dir, "contact_difference_barplot.png")
        plt.savefig(path)
        print(f"Saved bar plot to {path}")
        plt.close()
    else:
        plt.show()


# Heatmap of Contacts
def plot_contact_heatmap(df, save=False, output_dir="."):
    heatmap_data = df.groupby("PDB File")[["WT Contacts", "Variant Contacts", "Contact Difference"]].sum()
    plt.figure(figsize=(8, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="vlag", center=0, cbar_kws={'label': 'Contact Value'})
    plt.title("Heatmap of Contacts and Differences")
    plt.ylabel("PDB Variant")
    plt.tight_layout()

    if save:
        path = os.path.join(output_dir, "contact_heatmap.png")
        plt.savefig(path)
        print(f"Saved heatmap to {path}")
        plt.close()
    else:
        plt.show()


# Example Usage
def main():
    parser = argparse.ArgumentParser(description="Visualize WT and Variant contact differences from a CSV file.")
    parser.add_argument("-f", "--file", required=True, help="Path to the input CSV file")
    parser.add_argument("-s", "--save", action="store_true", help="Save plots instead of displaying")
    parser.add_argument("-o", "--output-dir", default=".", help="Directory to save plots (default: current directory)")

    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"Error: File not found - {args.file}")
        exit(1)

    if args.save and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    df = load_data(args.file)

    plot_contact_difference_bar(df, save=args.save, output_dir=args.output_dir)
    plot_contact_heatmap(df, save=args.save, output_dir=args.output_dir)

if __name__ == "__main__":
    main()
