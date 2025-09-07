# Automated Mutation Introduction and Analysis (AMIA) Workflow

The AMIA bioinformatics pipeline is an automated computational workflow designed for the effective prioritisation of potential drug resistance mutations by analysing their impact on protein folding and interactions, which is crucial for treatment success. To address this need, AMIA integrates a variety of structural analysis tools into a simplified and fully automated workflow, thereby optimising computational resources through automated data transformations to enhance scalability and reproducibility. This open-source pipeline automates key steps such as mutation introduction into protein structures, calculation of polar interaction changes, docking of ligand to WT and variant structures and analysis of protein fold energy using pre-established software tools. Furthermore, AMIA includes automated molecular dynamics analysis, which reduces the need for constant user input and output management often required by standalone tools. By facilitating the visualisation of mutation effects on protein structure and dynamic states, AMIA aids in prioritising variants for experimental validation and contributes to the development of improved treatment regimens against drug-resistant mutations.  

Detailed documentation can be found at: [https://kbrown3687524.github.io/amia/](https://kbrown3687524.github.io/amia/)

---

## Table of Contents
- [Installation](#installation)  
- [Usage](#usage)  
  - [Configuration File](#configuration-file)  
  - [Test Case](#test-case)  
- [Trajectory Analyses](#trajectory-analyses)  
- [Queries](#queries)  
- [Authors](#authors)  

---

## Installation

AMIA is installed via the command line using **conda** for environment management.  

### 1. Clone Repository
```bash
git clone https://github.com/kbrown3687524/amia
cd amia
```

### 2. Create Environment
```bash
conda env create -f amia_environment.yml
conda activate amia
```

### 3. Install Package
```bash
pip install .
```

### 4. Install MAESTRO
MAESTRO (v1.2.35) is a required standalone tool. Download and extract it into the main AMIA directory:

```
AMIA/
 ├─ MAESTRO_linux_x64/
 │   └─ maestro
```

---

## Usage

AMIA is executed using a **YAML configuration file** instead of passing long command-line arguments.  

### Configuration File

A typical `config.yaml` looks like this:  

```yaml
pdb_file: "/path/to/structure.pdb"
output_dir: "/path/to/output/"
mutations: "/path/to/mutations.csv"
mode: "multiple"

run_passer: true

passer_dir: "/path/to/output/"
passer_txt: "/path/to/output/passer_all_results.txt"
passer_html: "/path/to/output/passer_summary.html"

run_docking: true

smiles: "SMILES_STRING"
compound_name: "LigandName"
center: [X, Y, Z]
```

Run the workflow:

```bash
amia --config config.yaml
```

---

### Test Case

A ready-to-use test case is provided. Create a file called `config.yaml` with the following content:

```yaml
pdb_file: "/home/user/amia/test/HIV-1C_ZA/HIV_IN_1C_ZA_5U1C_model.pdb"
output_dir: "/home/user/amia/test/variant_outputs/"
mutations: "/home/user/amia/test/HIV-1C_ZA/mutations.csv"
mode: "multiple"

run_passer: true

passer_dir: "/home/user/amia/test/variant_outputs2"
passer_txt: "/home/user/amia/test/variant_outputs/passer_all_results.txt"
passer_html: "/home/user/amia/test/variant_outputs/passer_summary.html"

run_docking: true

smiles: "CC1=NN=C(O1)C(=O)NC(C)(C)C2=NC(=C(C(=O)N2C)O)C(=O)NCC3=CC=C(C=C3)F"
compound_name: "Aspirin"
center: [116.516, 139.229, 142.900]
```

Run with:

```bash
amia --config config.yaml
```

Outputs will be written to the directory specified under `output_dir`.

---

## Trajectory Analyses

After Phase 1 is complete, Molecular Dynamics (MD) simulations of WT and variant systems can be performed externally. Store repaired trajectories in this structure:

```
Trajectories/
 ├─ System1/
 │   ├─ System1.xtc (repaired)
 │   └─ System1.tpr
 ├─ System2/
     ├─ System2.xtc (repaired)
     └─ System2.tpr
```

Run trajectory statistics:

```bash
python3 trajstat.py --systems path/to/Trajectories                     --output_dir path/to/output_directory
```

Outputs include:
- RMSD, RMSF, Radius of Gyration  
- H-bond and salt bridge changes  
- PCA plots  

---

## Queries
For questions or support, contact:  
- **Keaghan Brown** — 3687524@myuwc.ac.za  
- **Dr. Ruben Cloete** — ruben@sanbi.ac.za  

---

## Authors
Keaghan Brown & Ruben Cloete  
