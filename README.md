# Automated Mutation Introduction and Analysis (AMIA) Workflow

The AMIA bioinformatics pipeline is an automated computational workflow designed for the effective prioritisation of potential drug resistance mutations by analysing their impact on protein folding and interactions, which is crucial for treatment success. To address this need, AMIA integrates a variety of structural analysis tools into a simplified and fully automated workflow, thereby optimising computational resources through automated data transformations to enhance scalability and reproducibility. This open-source pipeline automates key steps such as mutation introduction into protein structures, calculation of polar interaction changes, docking of ligands to WT and variant structures, and analysis of protein fold energy using pre-established software tools. Furthermore, AMIA includes automated molecular dynamics analysis, which reduces the need for constant user input and output management often required by standalone tools. By facilitating the visualisation of mutation effects on protein structure and dynamic states, AMIA aids in prioritising variants for experimental validation and contributes to the development of improved treatment regimens against drug-resistant mutations.  

Detailed documentation (under development and updates) can be found at: [https://kbrown3687524.github.io/amia/](https://kbrown3687524.github.io/amia/)

---

## Table of Contents
- [Installation](#installation)  
- [Pipeline Execution](#pipeline-execution)  
  - [Configuration File](#configuration-file)  
  - [Test Case](#test-case)  
- [Checkpointing and Resuming](#checkpointing-and-resuming)  
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

### 4. Install MAESTRO (Optional for Docking/Analysis)
MAESTRO (v1.2.35) is a required standalone tool for some analyses. Download and extract it into the main AMIA directory:

```
AMIA/
 ├─ MAESTRO_linux_x64/
 │   └─ maestro
```

---

## Pipeline Execution

AMIA is executed using the **`run_pipeline.py`** script, which reads a YAML configuration file and manages all pipeline steps automatically.

```bash
python run_pipeline.py --config config.yaml
```

### Options

| Option        | Description |
|---------------|-------------|
| `--config`    | Path to the YAML configuration file (required) |
| `--force`     | Ignore checkpoints and rerun all steps (optional) |

---

## Configuration File (`config.yaml`)

Below is a description of each parameter used in the configuration file:

- **`pdb_file`**: Path to the input protein structure in **PDB format**.  
- **`output_dir`**: Directory where all workflow results will be stored.  
- **`mutations`**: CSV file listing mutations to introduce.  
- **`mode`**: `"single"` = introduce each mutation individually, `"multiple"` = introduce all mutations together.  

Optional steps:

- **`run_maestroana`**: Boolean to run Maestro analysis.  
- **`run_passer`**: Boolean to run PASSER analysis for protein fold stability.  
- **`passer_dir`**: Output directory for PASSER results.  
- **`passer_txt`**: Tabulated PASSER results file path.  
- **`passer_html`**: PASSER summary HTML report.  
- **`run_docking`**: Boolean to run ligand docking.  
- **`smiles`**: SMILES string of the ligand.  
- **`compound_name`**: Descriptive ligand name.  
- **`center`**: `[X, Y, Z]` docking grid center coordinates.

#### Example Config

```yaml
pdb_file: "/home/user/amia/test/HIV-1C_ZA/HIV_IN_1C_ZA_5U1C_model.pdb"
output_dir: "/home/user/variant_outputs/"
mutations: "/home/user/amia/test/HIV-1C_ZA/mutations.csv"
mode: "multiple"

run_maestroana: false
run_passer: true
passer_dir: "/home/user/variant_outputs"
passer_txt: "/home/user/variant_outputs/passer_all_results.txt"
passer_html: "/home/user/variant_outputs/passer_summary.html"

run_docking: true
smiles: "CC1=NN=C(O1)C(=O)NC(C)(C)C2=NC(=C(C(=O)N2C)O)C(=O)NCC3=CC=C(C=C3)F"
compound_name: "Aspirin"
center: [116.516, 139.229, 142.900]
```

---

## Test Case

Use the example config above to test the pipeline. Run with:

```bash
python run_pipeline.py --config config.yaml
```

Outputs will be stored in `output_dir`, and checkpoints will allow the workflow to resume from the last completed step if interrupted.

---

## Checkpointing and Resuming

- Checkpoints are stored in `.checkpoints/last_completed.txt` inside the `output_dir`.  
- The pipeline will resume automatically from the next unfinished step.  
- To ignore checkpoints and rerun all steps, use:

```bash
python run_pipeline.py --config config.yaml --force
```

---

## Trajectory Analyses

After Phase 1 (mutation introduction, contacts, docking, PASSER) is complete, Molecular Dynamics (MD) simulations of WT and variant systems can be performed externally. Store repaired trajectories in this structure:

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
python3 trajstat.py --systems path/to/Trajectories --output_dir path/to/output_directory
```

Outputs include:

- RMSD, RMSF
- Radius of Gyration  
- H-bond and salt bridge changes  
- PCA plots
- Future updates may include SASA & MM-GBSA/MM-PBSA

---

## Queries

For questions or support, contact:  

- **Keaghan Brown** — 3687524@myuwc.ac.za  
- **Dr. Ruben Cloete** — ruben@sanbi.ac.za  

---

## Authors

Keaghan Brown & Ruben Cloete

