# Automated OneOpes Simulation Pipeline

A comprehensive pipeline for automating the setup and execution of OneOpes simulations. This pipeline handles everything from initial structure preparation to PLUMED files generation.

## Overview

This pipeline automates the following processes:
1. Structure file conversion and preparation
2. Gaussian input preparation and processing
3. Playmolecule output processing
4. Topology generation and system setup
5. Water model integration
6. Molecular docking for initial binding poses
7. System preparation for molecular dynamics
8. PLUMED configuration and analysis

## Prerequisites

### Conda Environment
A conda environment file (`environment.yml`) is provided to replicate the exact environment used for this pipeline. Create the environment using:

```bash
conda env create -f environment.yml
```

This will create a new environment called `OneOpes` with all required Python packages and dependencies.

### Additional Required Software
The following software must be installed separately and added to your system PATH:
- GROMACS with MPI support (gmx_mpi)
- PLUMED (version 2.9 or higher)
- Gaussian16
- Any additional external dependencies required by your specific system

### Required Software (included in conda environment)
- Python 3.7+
- PyMOL
- BioPython
- AutoDock Vina
- MGLTools
- AmberTools
- RDKit
- MDAnalysis
- Numpy/Scipy/Pandas

## Directory Structure

```
.
├── gaussian_results/                 # Gaussian results example
├── playmolecule_results/             # Playmolecule results example
├── structure_files/
│   ├── host/
│   └── guests/
├── scripts/
│   └── automated_protocol/           # All scripts must be run from this directory
│       ├── 1_prepare_gaussian.py
│       ├── 2_process_gaussian_folders.py
│       ├── 3_process_playmolecule.py
│       └── ...
├── system_parameters/
│   ├── charges/
│   ├── topologies/
│   ├── mdp_files/
│   └── water_force_fields/
└── system_preparation/
```

## Important Note on Script Usage

All scripts must be executed from within the `scripts/automated_protocol/` directory to ensure correct relative paths:

```bash
cd path/to/project/scripts/automated_protocol/
```

## Pipeline Steps

### 1. Structure Preparation (`1_prepare_gaussian.py`)
```bash
python 1_prepare_gaussian.py --calculation_type [vacuum|dielectric]
```
**Input:**
- Structure files (.cif, .pdb, .mol) in `structure_files/host/` and `structure_files/guests/`

**Output:**
- Gaussian input files (.com) in separate directories for each molecule
- SLURM submission scripts (gaussian_1.sh, gaussian_2.sh)
- Intermediate files (.xyz)
- Directory structure: `gaussian/{molecule_name}_gaussian/`

### 2. Process Gaussian Results (`2_process_gaussian_folders.py`)
```bash
python 2_process_gaussian_folders.py
```
**Input:**
- Gaussian output files (.log, .esp) in `gaussian_results/{molecule_name}_gaussian/`

**Output:**
- RESP charges (ANTECHAMBER_RESP2.OUT)
- Prep files (.prep)
- Force field modification files (frcmod)
- Directory: `system_parameters/charges/{output_folder_name}/`

### 3. Process PlayMolecule Results (`3_process_playmolecule.py`)
```bash
python 3_process_playmolecule.py
```
**Input:**
- PlayMolecule output folders in `playmolecule_results/`
- Each folder should contain parameters/GAFF2/ directory with force field files

**Output:**
- GROMACS topology files (.top)
- PDB files
- Directories: 
  - `system_parameters/topologies/{output_folder_name}/`
  - `system_parameters/pdb_files/`

### 4. Topology File Processing
#### 4.1 Host Topology Preprocessing (`4_1_preprocess_host_top_file.py`)
```bash
python 4_1_preprocess_host_top_file.py
```
**Input:**
- Host topology file (.top) from `system_parameters/topologies/{folder}/`

**Output:**
- Modified host topology file with:
  - Updated atom names (prefixed with 'y')
  - Modified residue names
  - Adjusted atom numbering

#### 4.2 Merge Topology Files (`4_merge_top_files.py`)
```bash
python 4_merge_top_files.py
```
**Input:**
- Processed host topology (.top)
- Guest topology files (.top)
- PDB files in `system_parameters/pdb_files/`

**Output:**
- Combined topology file (topol.top)
- Individual .itp files
- Directory: `system_preparation/{output_folder_name}/{host}_{guest}/`

### 5. Water Model Integration (`5_get_water_model.py`)
```bash
python 5_get_water_model.py
```
**Input:**
- Water model files from `system_parameters/water_force_fields/`
- Topology files in system preparation directories

**Output:**
- Modified topology file (topol_water.top) with integrated water model parameters
- Copied water model .itp file in each system directory
- Directory structure preserved with water model files in each relevant subfolder

### 6. Molecular Docking (`6_docking.py`)
```bash
python 6_docking.py
```
**Input:**
- Host PDB file
- Guest PDB file
- System preparation folder

**Output:**
- Docked guest-host complex
- PDBQT files for receptor and ligand
- Docking configuration file
- Docked poses in PDB format

**Features:**
- Automated search box calculation
- Structure preparation using MGLTools
- AutoDock Vina docking
- Conversion of docked poses to PDB format

### 7. System Preparation (`7_auto_prepare.py`)
```bash
python 7_auto_prepare.py [--water_points {3,4,5}]
```

**Input:**
- Host structure file (.pdb/.gro)
- Guest structure file (.pdb/.gro)
- Topology file (.top)
- Water model points specification

**Output:**
- Solvated system (solv.gro)
- Ion-added system (solv_ions.gro)
- Energy minimized structure (em.gro)
- NVT equilibrated structure (nvt.gro)
- NPT equilibrated structure (npt.gro)
- Index files (index.ndx)
- Position restraint files (posre.itp, posre_ligand.itp)
- MDP files copied from system_parameters/mdp_files/

**Features:**
- Runs from scripts/automated_protocol/ directory
- Automatically copies MDP files from system_parameters
- Handles multiple systems in parallel
- Interactive system selection

### 8. PLUMED Configuration and Analysis (`8_auto_plumed.py`)
```bash
python 8_auto_plumed.py system.pdb npt.gro [--double_funnel] [--other_side]
```
**Input:**
- Template PDB file (system.pdb)
- Equilibrated structure (npt.gro)
- Options for funnel configuration

**Output:**
- Final PLUMED input file (plumed.dat) containing:
  - Atom definitions and alignments
  - Collective variables and descriptors
  - Funnel restraints and wall definitions
  - OPES parameters
  - Coordination numbers and analysis metrics
- PyMOL visualization session containing:
  - Molecule representation
  - Center of mass indication
  - Principal axes visualization
  - Funnel and cylinder representations
  - Virtual atoms for analysis

## Usage Example

```bash
# Change to the scripts/automated_protocol directory
cd path/to/project/scripts/automated_protocol/

# 1. Prepare Gaussian calculations
python 1_prepare_gaussian.py --calculation_type vacuum

# 2. Process Gaussian results
python 2_process_gaussian_folders.py

# 3. Process PlayMolecule results
python 3_process_playmolecule.py

# 4. Process topology files
python 4_1_preprocess_host_top_file.py
python 4_merge_top_files.py

# 5. Add water model
python 5_get_water_model.py

# 6. Perform molecular docking
python 6_docking.py

# 7. Prepare system
python 7_auto_prepare.py --water_points 3

# 8. Configure PLUMED and generate analysis
python 8_auto_plumed.py system.pdb npt.gro --double_funnel
```

## Notes

- All scripts must be run from the `scripts/automated_protocol/` directory
- The pipeline assumes specific directory structure and file naming conventions
- Each script includes interactive prompts for necessary user input
- The pipeline is designed for host-guest systems but can be adapted for other cases
- All relative paths in scripts use `../../` to reference directories outside of scripts/automated_protocol/
- Ensure all external dependencies (GROMACS, PLUMED, etc.) are properly installed and accessible in your PATH before running the pipeline
