import sys
import os
import shutil

def extract_sections(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        writing = False
        for line in infile:
            if '[ moleculetype ]' in line:
                writing = True
            
            if writing:
                outfile.write(line)
            
            if '[ dihedrals ]' in line:
                # Continue writing until we hit a blank line or a new section
                for next_line in infile:
                    if next_line.strip() == '' or next_line.strip().startswith('['):
                        break
                    outfile.write(next_line)
                break

def get_section_content(file_content, section_name):
    content = []
    in_section = False
    for line in file_content:
        if f'[ {section_name} ]' in line:
            in_section = True
            content.append(line)
        elif in_section:
            if line.strip().startswith('['):
                break
            content.append(line)
    return content

def merge_files(file1, file2, output_file):
    with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output_file, 'w') as out:
        content1 = f1.readlines()
        content2 = f2.readlines()

        # Write [ defaults ] section from file1
        defaults = get_section_content(content1, 'defaults')
        out.writelines(defaults)

        # Write [ atomtypes ] section
        atomtypes1 = get_section_content(content1, 'atomtypes')
        atomtypes2 = [line for line in get_section_content(content2, 'atomtypes') if not line.strip().startswith(';') and not line.startswith("[")]
        
        # Remove trailing newlines from atomtypes1 and leading newlines from atomtypes2
        while atomtypes1 and atomtypes1[-1].strip() == '':
            atomtypes1.pop()
        while atomtypes2 and atomtypes2[0].strip() == '':
            atomtypes2.pop(0)
        
        out.writelines(atomtypes1)
        out.writelines(atomtypes2)

        # Write include statements
        out.write(f'#include "{os.path.basename(file1).split(".")[0]}.itp"\n\n')
        out.write('; Include Position restraint file\n')
        out.write('#ifdef POSRES\n')
        out.write('#include "posre.itp"\n')
        out.write('#endif\n\n')
        out.write(f'#include "{os.path.basename(file2).split(".")[0]}.itp"\n\n')
        out.write('; Include Position restraint file\n')
        out.write('#ifdef POSRES\n')
        out.write('#include "posre_ligand.itp"\n')
        out.write('#endif\n\n')

        # Write [ system ] section from file1
        system = get_section_content(content1, 'system')
        out.writelines(system)

        # Write [ molecules ] section
        molecules2 = get_section_content(content2, 'molecules')
        molecules1 = [line for line in get_section_content(content1, 'molecules') if not line.strip().startswith(';') and not line.startswith("[")]
        out.writelines(molecules2)
        out.writelines(molecules1)

def select_folder(base_path):
    folders = [f for f in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, f))]
    
    if not folders:
        print("No folders found in the specified directory.")
        sys.exit(1)
    
    print("Available folders:")
    for i, folder in enumerate(folders, 1):
        print(f"{i}. {folder}")
    
    while True:
        try:
            choice = int(input("Enter the number of the folder you want to select: "))
            if 1 <= choice <= len(folders):
                return folders[choice - 1]
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def select_host_file(folder_path):
    files = [f for f in os.listdir(folder_path) if f.endswith('.top')]
    
    if not files:
        print("No .top files found in the selected folder.")
        sys.exit(1)
    
    print("Available .top files:")
    for i, file in enumerate(files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Enter the number of the host file: "))
            if 1 <= choice <= len(files):
                return files[choice - 1]
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def get_output_folder_name():
    while True:
        name = input("Enter a name for the output folder: ")
        output_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation', name))
        
        if os.path.exists(output_dir):
            print(f"Warning: The folder '{name}' already exists.")
            confirm = input("Do you want to overwrite its contents? (yes/no): ").lower()
            if confirm == 'yes':
                return name
            else:
                print("Please choose a different name.")
        else:
            return name

def main():
    # Set the base directory to ../../system_parameters/topologies/
    base_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_parameters', 'topologies'))
    
    # Select folder
    selected_folder = select_folder(base_dir)
    folder_path = os.path.join(base_dir, selected_folder)
    
    # Select host file
    host_file = select_host_file(folder_path)
    host_file_path = os.path.join(folder_path, host_file)
    
    # Get all other .top files in the folder
    guest_files = [f for f in os.listdir(folder_path) if f.endswith('.top') and f != host_file]
    
    # Get output folder name from user
    output_folder_name = get_output_folder_name()
    
    # Create output directory
    output_base_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation', output_folder_name))
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Set the pdb files directory
    pdb_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_parameters', 'pdb_files'))
    
    for guest_file in guest_files:
        guest_file_path = os.path.join(folder_path, guest_file)
        output_dir = os.path.join(output_base_dir, f"{host_file.split('.')[0]}_{guest_file.split('.')[0]}")
        os.makedirs(output_dir, exist_ok=True)
        
        # Create .itp files
        extract_sections(host_file_path, os.path.join(output_dir, f"{host_file.split('.')[0]}.itp"))
        extract_sections(guest_file_path, os.path.join(output_dir, f"{guest_file.split('.')[0]}.itp"))
        
        # Create combined topol.top
        merge_files(host_file_path, guest_file_path, os.path.join(output_dir, "topol.top"))
        
        print(f"Created {host_file.split('.')[0]}.itp, {guest_file.split('.')[0]}.itp, and topol.top in {output_dir}")
        
        # Copy corresponding PDB files
        host_pdb = os.path.join(pdb_dir, f"{host_file.split('.')[0]}.pdb")
        guest_pdb = os.path.join(pdb_dir, f"{guest_file.split('.')[0]}.pdb")
        
        if os.path.exists(host_pdb):
            shutil.copy(host_pdb, os.path.join(output_dir, f"{host_file.split('.')[0]}.pdb"))
            print(f"Copied {host_pdb} to {output_dir}")
        else:
            print(f"Warning: {host_pdb} not found")
        
        if os.path.exists(guest_pdb):
            shutil.copy(guest_pdb, os.path.join(output_dir, f"{guest_file.split('.')[0]}.pdb"))
            print(f"Copied {guest_pdb} to {output_dir}")
        else:
            print(f"Warning: {guest_pdb} not found")

    print(f"All files have been created in {output_base_dir}")

if __name__ == "__main__":
    main()
