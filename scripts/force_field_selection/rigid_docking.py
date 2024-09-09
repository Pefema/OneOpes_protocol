import os
import glob
from rdkit import Chem
import numpy as np
from pymol import cmd

def convert_cif_to_pdb(cif_file):
    pdb_file = os.path.splitext(cif_file)[0] + '.pdb'
    try:
        cmd.reinitialize()
        cmd.load(cif_file, "molecule")
        cmd.save(pdb_file, "molecule", format="pdb")
        print(f"Converted {cif_file} to {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Error converting {cif_file} to PDB: {e}")
        return None

def load_and_rename_molecule(file_path, new_name):
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return None
    
    mol = Chem.MolFromPDBFile(file_path, removeHs=False)
    
    if mol is None:
        print(f"Error: Failed to load molecule from file: {file_path}")
        return None
    
    # Rename residues
    for atom in mol.GetAtoms():
        atom.GetPDBResidueInfo().SetResidueName(new_name[:3])
    
    print(f"Successfully loaded and renamed molecule from: {file_path}")
    return mol

def rigid_docking(ligand, host):
    """Perform a simple rigid docking by translating the ligand to the host's center."""
    host_coords = host.GetConformer().GetPositions()
    host_center = np.mean(host_coords, axis=0)
    
    ligand_coords = ligand.GetConformer().GetPositions()
    ligand_center = np.mean(ligand_coords, axis=0)
    
    translation = host_center - ligand_center
    
    new_ligand = Chem.Mol(ligand)
    conf = new_ligand.GetConformer()
    for i in range(conf.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_pos = pos.x + translation[0], pos.y + translation[1], pos.z + translation[2]
        conf.SetAtomPosition(i, new_pos)
    
    return new_ligand

def combine_molecules(mol1, mol2):
    """Combine two RDKit molecules into one."""
    combined = Chem.CombineMols(mol1, mol2)
    return combined

def save_molecule(mol, output_file):
    """Save an RDKit molecule to a PDB file."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    writer = Chem.PDBWriter(output_file)
    writer.write(mol)
    writer.close()
    print(f"Saved to {output_file}")

def process_host_guest_pair(host_file, guest_file, output_folder):
    host_name = os.path.splitext(os.path.basename(host_file))[0]
    guest_name = os.path.splitext(os.path.basename(guest_file))[0]
    
    output_dir = os.path.join(output_folder, f"{host_name}_{guest_name}")
    
    host_output = os.path.join(output_dir, f"{host_name}.pdb")
    guest_output = os.path.join(output_dir, f"{guest_name}.pdb")
    combined_output = os.path.join(output_dir, f"{host_name}_{guest_name}.pdb")
    
    # Load molecules and rename residues
    host = load_and_rename_molecule(host_file, host_name[:3])
    guest = load_and_rename_molecule(guest_file, guest_name[:3])
    
    if host is None or guest is None:
        print(f"Docking cannot proceed for {host_name} and {guest_name} due to errors in loading molecules.")
        return
    
    # Perform rigid docking
    print(f"Performing rigid docking of {guest_name} to {host_name}...")
    docked_guest = rigid_docking(guest, host)
    
    # Save the combined docked pose
    combined = combine_molecules(docked_guest, host)
    save_molecule(combined, combined_output)
    
    # Save the docked guest pose
    save_molecule(docked_guest, guest_output)
    
    # Save the host (unchanged)
    save_molecule(host, host_output)
    
    print(f"Rigid docking and file saving completed for {host_name} and {guest_name}.")

def get_output_folder(base_dir):
    system_prep_dir = os.path.join(base_dir, 'system_preparation')
    
    # Get existing folders
    existing_folders = [f for f in os.listdir(system_prep_dir) if os.path.isdir(os.path.join(system_prep_dir, f))]
    
    while True:
        print("\nAvailable folders in system_preparation:")
        for i, folder in enumerate(existing_folders, 1):
            print(f"{i}. {folder}")
        print(f"{len(existing_folders) + 1}. Create a new folder")
        
        choice = input("\nEnter the number of your choice or the name of a new folder: ")
        
        if choice.isdigit():
            choice = int(choice)
            if 1 <= choice <= len(existing_folders):
                return os.path.join(system_prep_dir, existing_folders[choice - 1])
            elif choice == len(existing_folders) + 1:
                new_folder = input("Enter the name for the new folder: ")
                new_path = os.path.join(system_prep_dir, new_folder)
                os.makedirs(new_path, exist_ok=True)
                return new_path
        else:
            new_path = os.path.join(system_prep_dir, choice)
            os.makedirs(new_path, exist_ok=True)
            return new_path
        
        print("Invalid choice. Please try again.")

def get_molecule_files(directory):
    pdb_files = glob.glob(os.path.join(directory, '*.pdb'))
    cif_files = glob.glob(os.path.join(directory, '*.cif'))
    
    all_files = {}
    
    # First, add all PDB files
    for pdb_file in pdb_files:
        base_name = os.path.splitext(os.path.basename(pdb_file))[0]
        all_files[base_name] = pdb_file
    
    # Then, convert CIF files to PDB if a PDB doesn't exist for that molecule
    for cif_file in cif_files:
        base_name = os.path.splitext(os.path.basename(cif_file))[0]
        if base_name not in all_files:
            pdb_file = convert_cif_to_pdb(cif_file)
            if pdb_file:
                all_files[base_name] = pdb_file
    
    return list(all_files.values())

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    host_dir = os.path.join(base_dir, 'structure_files', 'host')
    guest_dir = os.path.join(base_dir, 'structure_files', 'guests')
    
    # Get output folder from user
    output_base = get_output_folder(base_dir)
    print(f"Results will be saved in: {output_base}")
    
    # Get all host and guest files (converting CIF to PDB if necessary)
    host_files = get_molecule_files(host_dir)
    guest_files = get_molecule_files(guest_dir)
    
    # Process all combinations
    for host_file in host_files:
        for guest_file in guest_files:
            process_host_guest_pair(host_file, guest_file, output_base)
    
    print("All docking processes completed.")

if __name__ == "__main__":
    main()
