import os
import sys
import shutil

def modify_topology_file(input_file, new_name):
    # Read the entire input file
    with open(input_file, 'r') as f:
        content = f.readlines()

    modified_content = []
    current_section = None

    # Process the file line by line
    for line in content:
        # Check if this line starts a new section
        if line.strip().startswith('['):
            current_section = line.strip()
            modified_content.append(line)
        elif current_section == '[ atomtypes ]':
            # Process the atomtypes section
            if line.strip() and not line.strip().startswith(';'):
                parts = line.split()
                if len(parts) > 1 and not parts[0].startswith(';'):
                    parts[0] = 'y' + parts[0]  # Add 'y' prefix to atom name
                    parts.pop(1)  # Remove at.num column
                    modified_line = '\t'.join(parts) + '\n'
                    modified_content.append(modified_line)
                else:
                    modified_content.append(line)  # Keep comment lines as is
            else:
                modified_content.append(line)
        elif current_section == '[ moleculetype ]':
            # Replace 'MOL' with the new name in the moleculetype section
            modified_line = line.replace('MOL', new_name)
            modified_content.append(modified_line)
        elif current_section == '[ atoms ]':
            # Process the atoms section
            if line.strip() and not line.strip().startswith(';'):
                parts = line.split()
                if len(parts) > 4 and not parts[0].startswith(';'):
                    parts[3] = new_name  # Change residue to new_name
                    parts[1] = 'y' + parts[1]  # Add 'y' prefix to atom type
                    modified_line = '\t'.join(parts) + '\n'
                    modified_content.append(modified_line)
                else:
                    modified_content.append(line)  # Keep comment lines as is
            else:
                modified_content.append(line)
        elif current_section == '[ molecules ]':
            # Replace 'MOL' with the new name in the molecules section
            if 'MOL' in line:
                modified_line = line.replace('MOL', new_name)
                modified_content.append(modified_line)
            else:
                modified_content.append(line)
        else:
            # For all other sections, keep the lines unchanged
            modified_content.append(line)

    # Create a temporary file with the modified content
    temp_file = input_file + '.temp'
    with open(temp_file, 'w') as f:
        f.writelines(modified_content)

    return temp_file

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

def select_file(folder_path):
    files = [f for f in os.listdir(folder_path) if f.endswith('.top')]
    
    if not files:
        print("No .top files found in the selected folder.")
        sys.exit(1)
    
    print("Available .top files:")
    for i, file in enumerate(files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Enter the number of the file you want to process: "))
            if 1 <= choice <= len(files):
                return files[choice - 1]
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def main():
    # Get the current working directory
    current_dir = os.getcwd()
    
    # Set the base directory to ../../system_parameters/topologies/
    base_dir = os.path.abspath(os.path.join(current_dir, '..', '..', 'system_parameters', 'topologies'))
    
    print(f"Working with topologies in: {base_dir}")
    
    # Select folder
    selected_folder = select_folder(base_dir)
    folder_path = os.path.join(base_dir, selected_folder)
    
    # Select file
    selected_file = select_file(folder_path)
    input_file = os.path.join(folder_path, selected_file)
    
    # Prompt user for the new molecule name
    new_name = input("Enter the new molecule name: ")
    
    # Call the function to modify the topology file
    temp_file = modify_topology_file(input_file, new_name)
    
    # Replace the original file with the modified file
    shutil.move(temp_file, input_file)
    
    print(f"Modified topology file has replaced the original file: {input_file}")
    print(f"The molecule name has been changed to: {new_name}")

# This condition is true if the script is run directly (not imported as a module)
if __name__ == "__main__":
    main()
