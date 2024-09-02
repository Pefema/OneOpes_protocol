## Scripts Usage (WIP):

### 1. prepare_gaussian.py
Prepares SDF, CIF, MOL, or PDB files for use with Gaussian.
1. Place all files you want to prepare for Gaussian (in SDF, CIF, or PDB format) in a folder. If multiple formats are provided for the same file name, only the first found (in the order: SDF, CIF, PDB) will be processed. Add prepare_gaussian.py to this folder.
2. Modify the functions `write_first_sbatch_file` and `write_second_sbatch_file` to generate the appropriate SLURM files for your desired cluster.
3. Run the script with `python3 prepare_gaussian.py`. It will create a subfolder for each processed file.
4. Run Gaussian on your desired cluster using the generated input files.

Requirements: PyMol cmd must be installed (https://anaconda.org/conda-forge/pymol-open-source)

### 2. process_gaussian_folders.py
Processes Gaussian results to generate RESP2.OUT files.
1. Prepare the results in subfolders in the same format as the output from prepare_gaussian.py.
2. Run `process_gaussian_folders.py` in the folder containing the Gaussian results.

Requirements: AmberTools (antechamber) must be installed (https://ambermd.org/GetAmber.php#ambertools)

### 3. process_playmolecule.py
Processes PlayMolecule results to obtain topology files.
1. Place the output folder from PlayMolecule in a designated folder.
2. Run `process_playmolecule.py`.

Requirements: parmed and tleap (AmberTools) must be installed

### 4. charges_to_top.py
Incorporates charges from Gaussian RESP2.OUT files into the .top files from PlayMolecule.
1. Place the results folder from process_gaussian_folders.py and the folder with results from process_playmolecule.py in the same directory.
2. Run the script: `python3 charges_to_top.py gaussian_folder playmolecule_folder`.
3. The resulting topologies will be created in a new folder called "topologies".

### 5. preprocess_host_top_file.py
Modifies atom type names and molecule type names from PlayMolecule to ensure compatibility.
1. Run it with: `python3 preprocess_host_top_file.py <host.top>`

This script:
- Adds a prefix to all atom names to avoid compatibility errors with guest files.
- Modifies the molecule name from MOL to the file name (e.g., CB8, TEMOA, etc.).

### 6. merge_top_files.py
Merges host and guest topology files.
1. Run it with: `python3 merge_top_files.py <host_file.top> <guest_file.top>`

This script:
- Creates .itp files for both the host and the guest.
- Combines the [ atomtypes ] and other sections of the topology in a topol.top file.
- Generates three files: guest_file.itp, host_file.itp, and topol.top.

