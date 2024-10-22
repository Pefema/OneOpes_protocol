# Automated OneOpes Simulation Pipeline

A comprehensive pipeline for automating the setup and execution of OneOpes simulations. This pipeline handles everything from initial structure preparation to PLUMED files generation.

## Overview

This pipeline automates the following processes:
1. Structure file conversion and preparation
2. Gaussian input preparation and processing
3. Playmolecule output processing
4. Topology generation and system setup
5. System preparation for molecular dynamics
6. PLUMED files generation and multireplica generation

## Prerequisites

### Required Software
- GROMACS (with MPI support)
- Gaussian16
- PLUMED
- AmberTools
- Python 3.7+
- PyMOL
- BioPython

## Directory Structure

```
.
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
│   └── water_force_fields/
└── system_preparation/
```

## Important Note on Script Usage

All scripts must be executed from within the `scripts/automated_protocol/` directory. This ensures that relative paths to other directories are correctly maintained. Before running any script, make sure to:

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

### 6. PLUMED Configuration
#### 6.1 Initial PLUMED Setup (`6_plumed_generator_start.py`)
```bash
python 6_plumed_generator_start.py
```
**Input:**
- NPT equilibrated structure (npt.gro)

**Output:**
- Initial PLUMED info file (plumed_info.txt) containing:
  - Group definitions
  - Selected atoms for collective variables
  - Water oxygen atoms information

#### 6.2 Template Structure Renumbering (`6_2_renumber.py`)
```bash
python 6_2_renumber.py
```
**Input:**
- Template PDB files (guest.pdb, host_template.pdb)

**Output:**
- Renumbered template PDB file (template_updated.pdb) with consistent atom numbering

### 7. System Preparation (`7_auto_prepare.py`)
```bash
python 7_auto_prepare.py protein_file ligand_file topol_file [--water_points {3,4,5}]
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

### 8. PLUMED Configuration Finalization (`8_auto_plumed.py`)
```bash
python 8_auto_plumed.py system.pdb npt.gro [--double_funnel] [--other_side]
```
**Input:**
- Template PDB file (system.pdb)
- Equilibrated structure (npt.gro)
- Options for funnel configuration

**Output:**
- Final PLUMED input file (plumed.dat) containing:
  - Atom definitions
  - Collective variables
  - Funnel restraints
  - OPES parameters
- PyMOL session file for visualization

## Usage Example

Make sure you're in the correct scripts directory before running any commands:

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

# 6. Set up PLUMED
python 6_plumed_generator_start.py
python 6_2_renumber.py

# 7. Prepare system
python 7_auto_prepare.py protein.pdb ligand.pdb topol.top --water_points 3

# 8. Finalize PLUMED configuration
python 8_auto_plumed.py system.pdb npt.gro --double_funnel
```

## Notes

- All scripts must be run from the `scripts/automated_protocol/` directory
- The pipeline assumes specific directory structure and file naming conventions
- Each script includes interactive prompts for necessary user input
- The pipeline is designed for host-guest systems but can be adapted for other cases
- All relative paths in scripts use `../../` to reference directories outside of scripts/automated_protocol/
