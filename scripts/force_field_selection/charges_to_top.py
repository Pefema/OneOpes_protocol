import os
import re
import shutil
import sys

# Function to extract q(init) values from the .OUT file
def extract_q_init_values(out_file):
    q_init_values = []
    q_init_section_started = False
    
    with open(out_file, 'r') as file:
        for line in file:
            if 'no.  At.no.' in line:
                q_init_section_started = True
                continue
            
            if q_init_section_started and re.match(r'\s*\d+\s+\d+\s+[+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?', line):
                parts = line.split()
                q_init_values.append(parts[2])  # Index 2 for q(init)
    
    return q_init_values

# Function to count the number of valid atom entries in the MOL.top file's [ atoms ] section
def count_atoms_in_top(mol_top_file):
    atom_count = 0
    in_atoms_section = False
    
    with open(mol_top_file, 'r') as file:
        for line in file:
            if '[ atoms ]' in line:
                in_atoms_section = True
                continue
            
            if in_atoms_section:
                if line.strip() == "" or line.startswith('['):
                    break
                
                # Count only lines that start with a number (indicating an atom entry)
                if re.match(r'^\s*\d+', line):
                    atom_count += 1
    
    return atom_count

# Function to replace charges in the MOL.top file while preserving the formatting
def replace_charges_in_top(mol_top_file, q_init_values, output_file):
    mol_top_modified = []
    in_atoms_section = False
    q_init_values_iter = iter(q_init_values)
    
    with open(mol_top_file, 'r') as file:
        for line in file:
            if '[ atoms ]' in line:
                in_atoms_section = True
                mol_top_modified.append(line)
                continue
            
            if in_atoms_section:
                if line.strip() == "" or line.startswith('['):
                    in_atoms_section = False
                
                if in_atoms_section:
                    # Process only lines that start with a number (indicating an atom entry)
                    if re.match(r'^\s*\d+', line):
                        parts = re.split(r'(\s+)', line)  # Split by whitespace, keeping separators
                        if len(parts) > 13:  # Ensure there are enough parts to modify
                            try:
                                parts[14] = next(q_init_values_iter)  # Replace charge value (9th field, accounting for separators)
                                mol_top_modified.append("".join(parts))
                            except StopIteration:
                                print(f"Error: The number of q(init) values is less than the number of atoms in {mol_top_file}.")
                                return
                        else:
                            mol_top_modified.append(line)
                    else:
                        mol_top_modified.append(line)
                else:
                    mol_top_modified.append(line)
            else:
                mol_top_modified.append(line)
    
    with open(output_file, 'w') as file:
        file.writelines(mol_top_modified)

# Main function to process all folders
def process_folders(gaussian_folder, playmolecule_folder):
    # Create the topologies folder if it doesn't exist
    topologies_folder = os.path.join(os.getcwd(), "topologies")
    os.makedirs(topologies_folder, exist_ok=True)
    
    # Iterate over all subfolders in the gaussian directory
    for folder_name in os.listdir(gaussian_folder):
        gaussian_subfolder = os.path.join(gaussian_folder, folder_name)
        if os.path.isdir(gaussian_subfolder):
            out_file_path = os.path.join(gaussian_subfolder, 'ANTECHAMBER_RESP2.OUT')
            if os.path.exists(out_file_path):
                # Extract the q(init) values from the .out file
                q_init_values = extract_q_init_values(out_file_path)
                
                # Find the corresponding subfolder in playmolecule_results
                base_folder_name = folder_name.replace('_gaussian', '')
                corresponding_pm_subfolder = os.path.join(playmolecule_folder, base_folder_name)
                top_file_path = os.path.join(corresponding_pm_subfolder, 'parameters', 'GAFF2', 'MOL.top')
                
                if os.path.exists(top_file_path):
                    # Check the number of atoms in the .top file
                    atom_count = count_atoms_in_top(top_file_path)
                    
                    if len(q_init_values) != atom_count:
                        print(f"Error: The number of q(init) values in {out_file_path} ({len(q_init_values)}) "
                              f"does not match the number of atoms in {top_file_path} ({atom_count}).")
                    else:
                        output_file = top_file_path.replace('MOL.top', 'MOL_modified.top')
                        replace_charges_in_top(top_file_path, q_init_values, output_file)
                        print(f"Modified topology file saved as {output_file}.")
                        
                        # Copy the modified file to the topologies folder with the subfolder name
                        destination_file = os.path.join(topologies_folder, f"{base_folder_name}.top")
                        shutil.copy(output_file, destination_file)
                        print(f"Copied {output_file} to {destination_file}.")
                else:
                    print(f"Warning: Corresponding MOL.top file not found for {out_file_path}")
            else:
                print(f"Warning: ANTECHAMBER_RESP2.OUT file not found in {gaussian_subfolder}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <gaussian folder> <playmolecule_results folder>")
    else:
        gaussian_folder = sys.argv[1]
        playmolecule_folder = sys.argv[2]
        process_folders(gaussian_folder, playmolecule_folder)
