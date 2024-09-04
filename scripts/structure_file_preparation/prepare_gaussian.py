from rdkit import Chem
from rdkit.Chem import AllChem
import CifFile
import os
import glob
import shutil
import argparse
from pymol import cmd

def cif_to_rdkit_mol(cif_file_path):
    cif = CifFile.ReadCif(cif_file_path)
    block = cif.first_block()
    
    # Using fields based on the provided structure
    atom_symbol_key = '_chem_comp_atom.type_symbol'
    x_coord_key = '_chem_comp_atom.model_cartn_x'
    y_coord_key = '_chem_comp_atom.model_cartn_y'
    z_coord_key = '_chem_comp_atom.model_cartn_z'
    
    if atom_symbol_key not in block:
        raise KeyError(f"No suitable atom symbol field found in CIF file. Available fields: {list(block.keys())}")
    
    # Extract atomic information
    atoms = block[atom_symbol_key]
    x_coords = block[x_coord_key]
    y_coords = block[y_coord_key]
    z_coords = block[z_coord_key]
    
    # Create RDKit molecule
    mol = Chem.RWMol()
    atom_indices = {}
    
    for i, atom_symbol in enumerate(atoms):
        atom = Chem.Atom(atom_symbol)
        idx = mol.AddAtom(atom)
        atom_indices[i] = idx
    
    # Add formal charges if available in CIF file
    if '_chem_comp_atom.charge' in block:
        charges = block['_chem_comp_atom.charge']
        for i, charge in enumerate(charges):
            mol.GetAtomWithIdx(atom_indices[i]).SetFormalCharge(int(float(charge)))
    
    # Set coordinates
    conf = Chem.Conformer(len(atoms))
    for i, (x, y, z) in enumerate(zip(x_coords, y_coords, z_coords)):
        conf.SetAtomPosition(atom_indices[i], (float(x), float(y), float(z)))
    
    mol.AddConformer(conf)
    
    return mol

def save_mol_as_sdf(mol, output_sdf_path):
    w = Chem.SDWriter(output_sdf_path)
    w.write(mol)
    w.close()

def convert_all_cif_to_sdf(directory):
    # Iterate over all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith('.cif'):
            cif_file_path = os.path.join(directory, filename)
            output_sdf_path = os.path.splitext(cif_file_path)[0] + '.sdf'
            
            if os.path.exists(output_sdf_path):
                print(f"Skipping conversion for {cif_file_path} as {output_sdf_path} already exists.")
            else:
                try:
                    mol = cif_to_rdkit_mol(cif_file_path)
                    save_mol_as_sdf(mol, output_sdf_path)
                    print(f"Successfully converted {cif_file_path} to {output_sdf_path}")
                except Exception as e:
                    print(f"Failed to convert {cif_file_path}: {e}")

def extract_total_charge(sdf_file_path):
    total_charge = 0
    with open(sdf_file_path, 'r') as f:
        for line in f:
            if line.startswith("M  CHG"):
                line = line.replace("M  CHG", "")
                parts = line.strip().split()[2::2]
                charges = sum([int(i) for i in parts])
                total_charge += charges
    return total_charge

def convert_sdf_to_xyz(sdf_file_path):
    cmd.reinitialize()
    total_charge = extract_total_charge(sdf_file_path)
    base_name = os.path.splitext(sdf_file_path)[0]
    xyz_file_path = f"{base_name}.xyz"
    cmd.load(sdf_file_path, "molecule")
    cmd.save(xyz_file_path, "molecule", format="xyz")
    print(f"Conversion complete. XYZ file saved at {xyz_file_path}")
    return xyz_file_path, total_charge

def convert_xyz_to_gaussian(input_file, charge):
    output_file = input_file.replace('.xyz', '.com')
    with open(output_file, 'w') as f_out:
        f_out.write("%Chk={}.chk\n".format(input_file.replace('.xyz', '')))
        f_out.write("%nprocshared=8\n")
        f_out.write("%mem=1GB\n")
        f_out.write("# B3LYP/6-31G** opt=Modredundant\n")
        f_out.write("\n")
        f_out.write("{}\n".format(input_file.replace('.xyz', '')))
        f_out.write("\n")
        f_out.write("{}  1\n".format(charge))
        with open(input_file, 'r') as f_in:
            lines = f_in.readlines()[2:]
            f_out.writelines(lines)
        f_out.write("\n")
        f_out.write(" * * * * F\n")
    print(f"Conversion complete. Gaussian input file saved as {output_file}")

def write_second_com_file(input_file, calculation_type):
    output_file = input_file.replace('.xyz', '_2.com')
    with open(output_file, 'w') as f_out:
        f_out.write("%Chk={}.chk\n".format(input_file.replace('.xyz', '')))
        f_out.write("%nprocshared=8\n")
        f_out.write("%mem=1GB\n")
        
        if calculation_type == "dielectric":
            f_out.write("# HF/6-31G** SCRF=(Solvent=Water) POP=MK iop(6/50=1) geom=allcheck guess=read\n")
        else:  # Default to "vacuum"
            f_out.write("# HF/6-31G** POP=MK iop(6/50=1) geom=allcheck guess=read\n")
        
        f_out.write("\n")
        f_out.write("{}.esp\n".format(input_file.replace('.xyz', '')))
    print("Second com file successfully written.")

def write_first_sbatch_file(input_file):
    output_file = "gaussian_1.sh"
    with open(output_file, 'w') as f_out:
        f_out.write("#!/bin/env bash\n")
        f_out.write("#SBATCH -J Gaus{}\n".format(input_file.replace('.xyz', '')))
        f_out.write("#SBATCH -e gaus_%j.e\n")
        f_out.write("#SBATCH -o gaus_%j.o\n")
        f_out.write("#SBATCH --nodes=1\n")
        f_out.write('#SBATCH --constraint="V4|V5|V6|V9"\n')
        f_out.write("#SBATCH --ntasks=1\n")
        f_out.write("#SBATCH --cpus-per-task=8\n")
        f_out.write("#SBATCH -t 12:00:00\n")
        f_out.write("#SBATCH -p shared-cpu\n")
        f_out.write("#SBATCH --mem=10G\n")
        f_out.write("\n")
        f_out.write("ml gaussian\n")
        f_out.write("\n")
        f_out.write("g16 < {}.com > {}.log".format(input_file.replace('.xyz', ''), input_file.replace('.xyz', '')))
    print("First sbatch file successfully written.")

def write_second_sbatch_file(input_file):
    output_file = "gaussian_2.sh"
    with open(output_file, 'w') as f_out:
        f_out.write("#!/bin/env bash\n")
        f_out.write("#SBATCH -J Gaus_2_{}\n".format(input_file.replace('.xyz', '')))
        f_out.write("#SBATCH -e gaus_%j.e\n")
        f_out.write("#SBATCH -o gaus_%j.o\n")
        f_out.write("#SBATCH --nodes=1\n")
        f_out.write('#SBATCH --constraint="V4|V5|V6|V9"\n')
        f_out.write("#SBATCH --ntasks=1\n")
        f_out.write("#SBATCH --cpus-per-task=8\n")
        f_out.write("#SBATCH -t 12:00:00\n")
        f_out.write("#SBATCH -p shared-cpu\n")
        f_out.write("#SBATCH --mem=10G\n")
        f_out.write("\n")
        f_out.write("ml gaussian\n")
        f_out.write("\n")
        f_out.write("g16 < {}_2.com > {}_2.log".format(input_file.replace('.xyz', ''), input_file.replace('.xyz', '')))
    print("Second sbatch file successfully written.")

def process_single_sdf_file(sdf_file_path, calculation_type):
    xyz_file, total_charge = convert_sdf_to_xyz(sdf_file_path)
    convert_xyz_to_gaussian(xyz_file, total_charge)
    write_second_com_file(xyz_file, calculation_type)
    write_first_sbatch_file(xyz_file)
    write_second_sbatch_file(xyz_file)
    
    # Create a directory for the sdf file
    folder_name = os.path.splitext(os.path.basename(sdf_file_path))[0] + "_gaussian"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    # Move all generated files and the sdf file to the new directory
    for generated_file in [xyz_file, xyz_file.replace('.xyz', '.com'), xyz_file.replace('.xyz', '_2.com'), "gaussian_1.sh", "gaussian_2.sh"]:
        shutil.move(generated_file, os.path.join(folder_name, generated_file))
    shutil.move(sdf_file_path, os.path.join(folder_name, os.path.basename(sdf_file_path)))

def pdb_to_rdkit_mol(pdb_file_path):
    # Read PDB file and create RDKit molecule
    mol = Chem.MolFromPDBFile(pdb_file_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read PDB file: {pdb_file_path}")
    return mol

def convert_all_to_sdf(directory):
    # Iterate over all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith(('.cif', '.pdb')):
            input_file_path = os.path.join(directory, filename)
            output_sdf_path = os.path.splitext(input_file_path)[0] + '.sdf'
            
            if os.path.exists(output_sdf_path):
                print(f"Skipping conversion for {input_file_path} as {output_sdf_path} already exists.")
            else:
                try:
                    if filename.endswith('.cif'):
                        mol = cif_to_rdkit_mol(input_file_path)
                    elif filename.endswith('.pdb'):
                        mol = pdb_to_rdkit_mol(input_file_path)
                    
                    save_mol_as_sdf(mol, output_sdf_path)
                    print(f"Successfully converted {input_file_path} to {output_sdf_path}")
                except Exception as e:
                    print(f"Failed to convert {input_file_path}: {e}")

def copy_structure_files(source_dir, destination_dir):
    file_types = ['*.sdf', '*.cif', '*.pdb', '*.mol']
    for file_type in file_types:
        for file in glob.glob(os.path.join(source_dir, file_type)):
            shutil.copy(file, destination_dir)
            print(f"Copied {file} to {destination_dir}")

def remove_structure_files(destination_dir):
    file_types = ['*.sdf', '*.cif', '*.pdb', '*.mol']
    for file_type in file_types:
        for file in glob.glob(os.path.join(destination_dir, file_type)):
            os.remove(file)
            print(f"Removed {file}")

def main():
    parser = argparse.ArgumentParser(description="Process CIF, PDB, SDF, and MOL files to prepare Gaussian inputs.")
    parser.add_argument("--calculation_type", type=str, choices=["vacuum", "dielectric"], default="vacuum", help="Type of calculation: vacuum or dielectric (default: vacuum)")
    
    args = parser.parse_args()
    
    current_directory = os.getcwd()
    
    # Create the gaussian folder two levels up
    gaussian_folder = os.path.abspath(os.path.join(current_directory, '..', '..', 'gaussian'))
    os.makedirs(gaussian_folder, exist_ok=True)
    
    # Define paths for host and guest structure files
    host_dir = os.path.abspath(os.path.join(current_directory, '..', '..', 'structure_files', 'host'))
    guests_dir = os.path.abspath(os.path.join(current_directory, '..', '..', 'structure_files', 'guests'))
    
    # Copy structure files from host and guests directories to the gaussian folder
    copy_structure_files(host_dir, gaussian_folder)
    copy_structure_files(guests_dir, gaussian_folder)
    
    # Change working directory to the gaussian folder
    os.chdir(gaussian_folder)
    
    # First, convert CIF, PDB, and MOL files to SDF if needed
    convert_all_to_sdf(gaussian_folder)
    
    # Process all SDF files in the gaussian folder
    for sdf_file in glob.glob("*.sdf"):
        print(f"Processing {sdf_file}")
        process_single_sdf_file(sdf_file, args.calculation_type)

    remove_structure_files(gaussian_folder)
if __name__ == "__main__":
    main()