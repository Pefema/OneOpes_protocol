import sys
import os

def modify_topology_file(input_file, output_file, new_name):
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

    # Write the modified content to the output file
    with open(output_file, 'w') as f:
        f.writelines(modified_content)

def main():
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_filename>")
        sys.exit(1)

    # Get the input filename from command-line argument
    input_file = sys.argv[1]
    
    # Extract the part before the first underscore to use as the new name
    new_name = input_file.split('_')[0]
    
    # Generate the output filename
    output_file = f"{new_name}_modified.top"

    # Call the function to modify the topology file
    modify_topology_file(input_file, output_file, new_name)
    
    print(f"Modified topology file saved as {output_file}")

# This condition is true if the script is run directly (not imported as a module)
if __name__ == "__main__":
    main()
