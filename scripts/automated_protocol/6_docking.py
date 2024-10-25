import os
import sys
import subprocess
import numpy as np
import MDAnalysis as mda
import glob
import argparse

def prepare_receptor(pdb_file):
    """Prepare receptor using MGLTools prepare_receptor4 script with minimal processing."""
    # Save original file with _not_docked suffix
    not_docked_file = pdb_file.replace('.pdb', '_not_docked.pdb')
    if not os.path.exists(not_docked_file):
        subprocess.run(['cp', pdb_file, not_docked_file], check=True)

    # Get paths
    abs_pdb_file = os.path.abspath(pdb_file)
    output_pdbqt = abs_pdb_file.replace('.pdb', '.pdbqt')
    
    # Get directory and filename
    working_dir = os.path.dirname(abs_pdb_file)
    pdb_filename = os.path.basename(abs_pdb_file)
    output_filename = os.path.basename(output_pdbqt)
    
    # Store current directory
    original_dir = os.getcwd()
    
    try:
        # Change to working directory
        os.chdir(working_dir)
        
        # Run prepare_receptor4 with local filenames
        cmd = f"prepare_receptor4.py -r {pdb_filename} -o {output_filename}"
        print(f"Running command in {working_dir}: {cmd}")
        
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
            
        # Verify the output file exists
        if not os.path.exists(output_filename):
            raise Exception(f"PDBQT file was not created: {output_filename}")
            
        return output_pdbqt
        
    finally:
        # Always return to original directory
        os.chdir(original_dir)

def prepare_ligand(pdb_file):
    """Prepare ligand using MGLTools prepare_ligand4 script while strictly preserving structure."""
    # Save original file with _not_docked suffix
    not_docked_file = pdb_file.replace('.pdb', '_not_docked.pdb')
    if not os.path.exists(not_docked_file):
        subprocess.run(['cp', pdb_file, not_docked_file], check=True)

    # Get paths
    abs_pdb_file = os.path.abspath(pdb_file)
    output_pdbqt = abs_pdb_file.replace('.pdb', '.pdbqt')
    
    # Get directory and filename
    working_dir = os.path.dirname(abs_pdb_file)
    pdb_filename = os.path.basename(abs_pdb_file)
    output_filename = os.path.basename(output_pdbqt)
    
    # Store current directory
    original_dir = os.getcwd()
    
    try:
        # Change to working directory
        os.chdir(working_dir)
        
        # Run prepare_ligand4 with local filenames
        cmd = f"prepare_ligand4.py -l {pdb_filename} -o {output_filename} -U -C -p -k -Z"
        print(f"Running command in {working_dir}: {cmd}")
        
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
            
        # Verify the output file exists
        if not os.path.exists(output_filename):
            raise Exception(f"PDBQT file was not created: {output_filename}")
            
        return output_pdbqt
        
    finally:
        # Always return to original directory
        os.chdir(original_dir)

def read_pdbqt_coordinates(pdbqt_file):
    """Read coordinates from PDBQT file to calculate box dimensions."""
    coords = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)

def calculate_box(pdbqt_file, margin=3.0):
    """Calculate search box dimensions based on molecule size."""
    coords = read_pdbqt_coordinates(pdbqt_file)
    
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    
    center = (min_coords + max_coords) / 2
    size = (max_coords - min_coords) + margin * 2
    
    return center, size

def create_config_file(center, size, receptor, ligand, output_dir, config_file="config.txt"):
    """Create Vina configuration file with parameters optimized for rigid docking."""
    config_path = os.path.join(output_dir, config_file)
    
    # Use basenames for receptor and ligand paths in config file
    receptor_name = os.path.basename(receptor)
    ligand_name = os.path.basename(ligand)
    
    with open(config_path, 'w') as f:
        f.write(f"receptor = {receptor_name}\n")
        f.write(f"ligand = {ligand_name}\n")

        f.write(f"\ncenter_x = {center[0]:.3f}")
        f.write(f"\ncenter_y = {center[1]:.3f}")
        f.write(f"\ncenter_z = {center[2]:.3f}")

        f.write(f"\nsize_x = {size[0]:.3f}")
        f.write(f"\nsize_y = {size[1]:.3f}")
        f.write(f"\nsize_z = {size[2]:.3f}")

        f.write("\nexhaustiveness = 32")
        f.write("\nnum_modes = 1")
        f.write("\nenergy_range = 1")
        f.write("\ncpu = 4")

    return config_path

def run_vina(config_file, output_dir, output_pdbqt="docked.pdbqt"):
    """Run AutoDock Vina using command line."""
    output_path = os.path.join(output_dir, output_pdbqt)
    working_dir = os.path.dirname(config_file)
    config_name = os.path.basename(config_file)
    output_name = os.path.basename(output_path)
    
    original_dir = os.getcwd()
    
    try:
        os.chdir(working_dir)
        cmd = f"vina --config {config_name} --out {output_name}"
        print(f"Running Vina in {working_dir}: {cmd}")
        
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
            
        if not os.path.exists(output_name):
            raise Exception(f"Docked PDBQT file was not created: {output_name}")
            
        return output_path
        
    finally:
        os.chdir(original_dir)

def convert_pdbqt_to_pdb_mda(pdbqt_file, output_pdb):
    """Convert PDBQT to PDB using MDAnalysis."""
    try:
        u = mda.Universe(pdbqt_file)
        with mda.Writer(output_pdb) as writer:
            writer.write(u)
    except Exception as e:
        print(f"Error during MDAnalysis conversion: {str(e)}")
        raise

def select_input_folder():
    """Select input folder from system_preparation directory."""
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    folders = [d for d in os.listdir(prep_dir)
               if os.path.isdir(os.path.join(prep_dir, d))]
    
    if not folders:
        print("Error: No folders found in system_preparation directory.")
        return None
    
    print("\nAvailable input folders:")
    for i, folder in enumerate(folders, 1):
        print(f"{i}. {folder}")
    
    while True:
        try:
            choice = int(input("Enter the number of the folder you want to process: "))
            if 1 <= choice <= len(folders):
                return folders[choice - 1]
            print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def select_host_pdb(subfolder_path):
    """Select host PDB file from the first subfolder."""
    pdb_files = [f for f in os.listdir(subfolder_path) if f.endswith('.pdb')]
    
    if not pdb_files:
        print(f"Error: No PDB files found in {subfolder_path}")
        return None
    
    print("\nAvailable PDB files:")
    for i, file in enumerate(pdb_files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Select the host PDB file to use for all docking runs: "))
            if 1 <= choice <= len(pdb_files):
                return pdb_files[choice - 1]
            print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def find_guest_pdb(subfolder_path, host_basename):
    """Find the guest PDB file based on the folder name pattern (HOST_GX)."""
    # Extract the guest number from the folder name
    folder_name = os.path.basename(subfolder_path)
    
    try:
        # Split the folder name by '_G' and take the last part as the guest number
        guest_number = folder_name.split('_G')[-1]
        guest_filename = f"G{guest_number}.pdb"
        
        if os.path.exists(os.path.join(subfolder_path, guest_filename)):
            return guest_filename
        else:
            print(f"Warning: Expected guest file {guest_filename} not found in {subfolder_path}")
            # List all PDB files in directory for debugging
            pdb_files = [f for f in os.listdir(subfolder_path) if f.endswith('.pdb')]
            print(f"Available PDB files in directory: {pdb_files}")
            return None
            
    except IndexError:
        print(f"Warning: Folder name {folder_name} doesn't match expected pattern HOST_GX")
        return None

def dock_molecules_in_folder(input_folder):
    """Process all subfolders using the selected host against each guest."""
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    base_dir = os.path.join(prep_dir, input_folder)

    # Get all subfolders
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not subfolders:
        print("Error: No subfolders found")
        return

    # Select host PDB from the first subfolder
    first_subfolder = os.path.join(base_dir, subfolders[0])
    host_basename = select_host_pdb(first_subfolder)
    if not host_basename:
        print("Error: No host PDB file selected")
        return

    host_file = os.path.join(first_subfolder, host_basename)
    print(f"\nSelected host file: {host_basename}")

    # Prepare the receptor once since it will be used for all docking runs
    try:
        print(f"\nPreparing receptor: {host_file}")
        receptor_pdbqt = prepare_receptor(host_file)
        print(f"Receptor preparation successful: {receptor_pdbqt}")
    except Exception as e:
        print(f"Error preparing receptor: {str(e)}")
        return

    # Process each subfolder
    for subfolder in subfolders:
        subfolder_path = os.path.join(base_dir, subfolder)
        print(f"\nProcessing subfolder: {subfolder}")

        # Find guest PDB file
        guest_basename = find_guest_pdb(subfolder_path, host_basename)
        if not guest_basename:
            print(f"Skipping {subfolder}: No guest PDB file found")
            continue

        guest_file = os.path.join(subfolder_path, guest_basename)
        print(f"Found guest file: {guest_basename}")

        try:
            # Copy host files to current subfolder if not already there
            if not os.path.exists(os.path.join(subfolder_path, host_basename)):
                subprocess.run(['cp', host_file, subfolder_path], check=True)
            if not os.path.exists(os.path.join(subfolder_path, os.path.basename(receptor_pdbqt))):
                subprocess.run(['cp', receptor_pdbqt, subfolder_path], check=True)

            print(f"Preparing ligand: {guest_file}")
            ligand_pdbqt = prepare_ligand(guest_file)

            print("Calculating search box...")
            center, size = calculate_box(os.path.join(subfolder_path, os.path.basename(receptor_pdbqt)))

            print("Creating configuration file...")
            config_path = create_config_file(center, size, 
                                          os.path.join(subfolder_path, os.path.basename(receptor_pdbqt)), 
                                          os.path.join(subfolder_path, os.path.basename(ligand_pdbqt)), 
                                          subfolder_path)

            print("Running AutoDock Vina...")
            docked_pdbqt = run_vina(config_path, subfolder_path)

            # Save docked complex with original filename
            output_pdb = os.path.join(subfolder_path, guest_basename)

            print("Converting docked complex to PDB...")
            convert_pdbqt_to_pdb_mda(docked_pdbqt, output_pdb)

            print(f"Docking completed successfully for {subfolder}")
            print(f"Complex saved as: {guest_basename}")

        except Exception as e:
            print(f"Error processing {subfolder}: {str(e)}")
            continue

def main():
    parser = argparse.ArgumentParser(description="Perform automated docking for all molecule pairs in the system.")
    args = parser.parse_args()

    print("\n=== Starting Automated Docking Process ===\n")
    
    input_folder = select_input_folder()
    if input_folder:
        dock_molecules_in_folder(input_folder)
    else:
        print("No valid input folder selected. Exiting.")
        return

    print("\nDocking process completed!")

if __name__ == "__main__":
    main()
