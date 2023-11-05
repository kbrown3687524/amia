# Automated Mutation Introduction and Analysis (AMIA) Workflow
This workflow is designed to automatically introduce mutation datasets into a specidied protein model and analyze the resultant effects by comparing the WT and Variant structures as well as to automate some the analyses involved in Molecular Dynamics Simulation Trajectory Statistics. This workflow is designed to be installed and executed on a Linux Device. This project was developed in fulfillment of the MSc Project:  Automated computational workflow to prioritize potential resistance variants in HIV Integrase Subtype C and CRF02_AG.



## Table of Contents

 * [Installation](#installation)
 * [Usage](#usage)
 * [Queries](#queries)
## Installation
This project was developed in a virtual environment on the command line and managed using conda for the respective packages and is recommended for installation purposes. 

### AMIA Download
The current release an be downloaded as follows:
```
git clone https://github.com/kbrown3687524/amia
```

### Conda venv Setup
Open a terminal with mamba and create a new env:
```
conda env create -f amia_environment.yml
conda activate amia_main
```
FoldX(4.0) is a standalone software tool that is required for this package to run successfully. It should be downloaded and extracted within the main AMIA direcotry to ensure successful integration with the workflow.
```
|-- AMIA Folder:
  |  
  |-- foldx:
      |-- foldx_4
      |-- yasaraPlugin.zip
      |-- rotabase.txt
```
Once the dependencies have been installed, a successful installation can be tested by using the dataset provided within the /test folder. The test can be initialized by providing chmod u+x access to the AMIA.py script and running the following:

```
./AMiA.py --pdb_file ~/test/HIV-1C_ZA/RLT_Model_Repair.pdb --mutations ~/test/HIV-1C_ZA/mutations.csv --output_dir ~/test/variant_outputs
```

## Usage

### Phase 1: Mutation Introduction
After successful installation, the scripts should now be available to execute individually or imported into other projects. The normal method to execute this workflow involves calling the AMIA script supplied with its respective arguments. The files provided to the workflow are in designated formats - see test files for further clarification on formats: 

* Structure File: Protein Data Bank (.pdb)
* Mutations File: Comma Separated Variable (.csv)

```
./AMIA.py --mode *single or multiple* --pdb_file *path/to/pdb_structure* --mutations *path/to/mutations_file* --output_dir *path/to/output_directory*
```
The ```--mode```  specifies whether the mutations from each subset present withint the mutation file should be introduced individually or together into the supplied protein structure. 
The ```--pdb_file``` specifies the path to the Protein File that the mutations will be introduced to. 
The ```--mutations``` specifies the path to the Mutations File that the mutations will be introduced to. 
The ```--output_dir``` specifies the directory that the respective output files will be stored in. 

This will automatically introduce the mutations from each subset into the protein file and store the respective outputs in the defined directory. From there the changes in contacts before and after mutation introductuin as well as the stability of the systems will be tabulated and stored for the user to visualize (.html).

### Phase 2: Trajectory Analyses
Once all the respective output files have been generated from the first phase of the workflow, the WT and Variant systems should then undergo Molecular Dynamis Simulations, afterwhich the repaired trajecotry and topology files shoudl be stored in a new directory under their respective subfolders. See below for simulation storage:

```
|-- Trajectories Folder:
  |
  |-- System 1:
  |   |-- System1.xtc (repaired)
  |   |-- System1.tpr
  |
  |-- System 2:
      |-- System2.xtc (repaired)
      |-- System2.tpr

``` 
The trajectory analyses and bond type changes are then calculated automatically using the trajectory file (.xtc) and the associated topology files (.tpr):

```
python3 trajana.py --systems *Main Trajectories Folder* --output_dir  *path/to/output_directory*
python3 bond_types.py --systems *Main Trajectories Folder* --output_dir  *path/to/output_directory*
```
The trajana.py script calculates and plots the Root Mean Square Deviation (RMSD) for the entire protein, Root Mean Square Fluctuation (RMSF) per protein chain and Radius of Gyration (rgyr) for the entire protein while the bond_types.py calculates Hydrogen Bond changes between DNA and Ligands within the system, the changes in ionic bonds (saltbridges) within the protein structure and the changes in hydrophic residue interactions. These outputs are then saved to the specified output directory.
## Queries
For any related queries please contact Mr. Keaghan Brown on 3687524@myuwc.ac.za
## Authors

- [@kbrown3687524](https://www.github.com/kbrown3687524)

