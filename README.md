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
- [Troubleshooting](#troubleshooting)
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
conda env create -f amia_environment_portable.yml
conda activate amia
```

The original `amia_environment.yml` is a Linux-specific export with Linux ABI
pins and an absolute environment prefix. Use the portable file on Windows,
macOS, and new Linux installations.

### 3. Install Package
```bash
python -m pip install .
```

### 4. Install optional native tools

Open Babel and AutoDock Vina are installed by the portable Conda environment.
FoldX and MAESTRO are separate native applications and are only required when
their corresponding pipeline steps are enabled. Install a build matching your
operating system.

```
Linux:
AMIA/
 ├─ MAESTRO_linux_x64/
 │   └─ maestro

OR

Windows:
AMIA/
 ├─ MAESTRO_win_x64/
 │   └─ maestro
```

The bundled `MAESTRO_linux_x64` executable is Linux-only. Windows and macOS
users must install a native MAESTRO build and set `AMIA_MAESTRO` to its full
path. Set `AMIA_FOLDX` to the native FoldX executable as well.

In Windows PowerShell:

```powershell
$env:AMIA_MAESTRO = "C:\Users\<user>\amia\MAESTRO_win_x64\maestro.exe"
```

In Windows Command Prompt:

```cmd
set "AMIA_MAESTRO=C:\Users\<user>\amia\MAESTRO_win_x64\maestro.exe"
```

In macOS/Linux shells:

```bash
export AMIA_MAESTRO="$HOME/Tools/maestro/maestro"
```

These variables apply only to the current terminal session. To confirm a
configured executable on Windows, run `Test-Path $env:AMIA_MAESTRO` in
PowerShell or `if exist "%AMIA_MAESTRO%" echo OK` in Command Prompt.

---

## Pipeline Execution

AMIA is executed using the installed **`amia`** command, which reads a YAML
configuration file and manages all pipeline steps automatically.

```bash
amia --config config.yaml
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
- **`run_trajstat`**: Boolean controlling whether trajectory analysis is run as
  an optional checkpointed AMIA stage.
- **`trajstat_systems`**: Directory containing system subdirectories with
  topology and trajectory files. Defaults to `output_dir` when omitted.
- **`trajstat_start_fr`**: Starting trajectory frame for TrajStat analyses.
  Defaults to `0` and avoids interactive input when managed by AMIA.

Relative paths in the configuration file are resolved relative to the directory
containing that configuration file. This makes the same configuration portable
between Windows, macOS, and Linux. Avoid machine-specific paths such as
`/home/user/...` or `C:\Users\...` when sharing a configuration.

#### Example Config

```yaml
pdb_file: "test/HIV-1C_ZA/HIV_IN_1C_ZA_5U1C_model.pdb"
output_dir: "variant_outputs"
mutations: "test/HIV-1C_ZA/mutations.csv"
mode: "multiple"

run_maestroana: false
run_passer: true
passer_dir: "variant_outputs"
passer_txt: "passer_all_results.txt"
passer_html: "passer_summary.html"

run_docking: true
smiles: "CC1=NN=C(O1)C(=O)NC(C)(C)C2=NC(=C(C(=O)N2C)O)C(=O)NCC3=CC=C(C=C3)F"
compound_name: "Aspirin"
center: [116.516, 139.229, 142.900]
```

---

## Test Case

Use the example config above to test the pipeline. Run with:

```bash
amia --config config.yaml
```

Outputs will be stored in `output_dir`, and checkpoints will allow the workflow to resume from the last completed step if interrupted.

---

## Checkpointing and Resuming

- Checkpoints are stored in `.checkpoints/last_completed.txt` inside the `output_dir`.  
- The pipeline will resume automatically from the next unfinished step.  
- To ignore checkpoints and rerun all steps, use:

```bash
amia --config config.yaml --force
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

## Troubleshooting

### `FileNotFoundError` contains `C:\home\...`

The configuration still contains a Linux path such as `/home/...`. Replace it
with a relative path, or use a valid Windows path. The repository test
configuration is portable:

```cmd
amia --config test\HIV-1C_ZA\config.yaml
```

### `KeyError: 0` in `mutintro.py`

This indicates that an older installed AMIA package is being used. Reinstall
the current checkout after pulling updates:

```cmd
python -m pip install --upgrade --force-reinstall .
```

### MAESTRO cannot be started

Check that `AMIA_MAESTRO` points to the executable file, not only its folder,
and use the syntax for the shell you are running. PowerShell uses `$env:NAME`,
Command Prompt uses `set "NAME=value"`, and Git Bash uses `export NAME=value`.
The Linux MAESTRO executable cannot run on Windows.

---

## Queries

For questions or support, contact:  

- **Keaghan Brown** — 3687524@myuwc.ac.za  
- **Dr. Ruben Cloete** — ruben@sanbi.ac.za  

---

## Authors

Keaghan Brown & Ruben Cloete

