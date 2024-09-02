import sys
import os

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

def main():
    if len(sys.argv) != 3:
        print("Usage: python merge_top_files.py <host_file.top> <guest_file.top>")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]
    
    # Create .itp files
    extract_sections(input_file1, f"{os.path.splitext(input_file1)[0]}.itp")
    extract_sections(input_file2, f"{os.path.splitext(input_file2)[0]}.itp")
    
    # Create combined topol.top
    merge_files(input_file1, input_file2, "topol.top")

    print(f"Created {os.path.splitext(input_file1)[0]}.itp, {os.path.splitext(input_file2)[0]}.itp, and topol.top")

if __name__ == "__main__":
    main()
