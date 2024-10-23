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
8. Initial PLUMED setup
9. Final PLUMED configuration and analysis

## Prerequisites

### Required Software
- GROMACS (with MPI support)
- Gaussian16
- PLUMED
- AmberTools
- Python 3.7+
- PyMOL
- BioPython
- AutoDock Vina
- MGLTools

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

Most scripts must be executed from within the `scripts/automated_protocol/` directory to ensure correct relative paths. The exception is `7_auto_prepare.py`, which must be run inside each individual system folder within `system_preparation/`. Before running scripts, navigate to:

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
- Topology file (topol_0.top) in system preparation directories

**Output:**
- Modified topology file (topol_water.top) with integrated water model parameters

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
# Must be run inside each system folder in system_preparation/
cd path/to/project/system_preparation/host_guest_folder/
python ../../scripts/automated_protocol/7_auto_prepare.py protein_file ligand_file topol_file [--water_points {3,4,5}]
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

### 8. Initial PLUMED Setup (`8_plumed_generator_start.py`)
```bash
python 8_plumed_generator_start.py
```
**Input:**
- NPT equilibrated structure (npt.gro)

**Output:**
- Initial PLUMED info file (plumed_info.txt) containing:
  - Group definitions
  - Selected atoms for collective variables
  - Water oxygen atoms information
  - Ligand selected atoms

### 9. PLUMED Configuration and Analysis (`9_auto_plumed.py`)
```bash
python 9_auto_plumed.py system.pdb npt.gro [--double_funnel] [--other_side]
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

**Analysis Features:**
- Principal component analysis of molecular structure
- Center of mass calculations
- Radius and height optimizations for funnel shape
- Coordination number analysis
- Angle and position restraints
- Virtual atom placement for path definition

## Usage Example

```bash
# First, change to the scripts/automated_protocol directory
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

# 7. Prepare system (run inside each system folder)
cd ../../system_preparation/host_guest_folder/
python ../../scripts/automated_protocol/7_auto_prepare.py protein.pdb ligand.pdb topol.top --water_points 3

# Return to scripts directory for final steps
cd ../../scripts/automated_protocol/

# 8. Initial PLUMED setup
python 8_plumed_generator_start.py

# 9. Configure PLUMED and generate analysis
python 9_auto_plumed.py system.pdb npt.gro --double_funnel
```

## Notes

- All scripts must be run from the `scripts/automated_protocol/` directory
- The pipeline assumes specific directory structure and file naming conventions
- Each script includes interactive prompts for necessary user input
- The pipeline is designed for host-guest systems but can be adapted for other cases
- All relative paths in scripts use `../../` to reference directories outside of scripts/automated_protocol/
