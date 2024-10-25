import os
import sys
import subprocess
import numpy as np
from Bio import PDB
import argparse
import pymol
from pymol import cmd, stored
from scipy.linalg import svd
from itertools import combinations
import re

def get_last_atom_number(filename):
    """Get the highest HETATM number from a PDB file."""
    last_number = 0
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                try:
                    number = int(line[6:11].strip())
                    last_number = max(last_number, number)
                except ValueError:
                    continue
    return last_number

def update_atom_numbers(input_filename, output_filename, increment):
    """
    Update atom numbers in PDB file by adding increment, excluding CONECT rows.
    Also sets occupancy and temperature factors to 1.00.
    
    Args:
        input_filename (str): Input PDB file path
        output_filename (str): Output PDB file path
        increment (int): Number to add to atom numbers
    """
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            # Skip CONECT lines completely
            if line.startswith('CONECT'):
                continue
            
            if line.startswith('HETATM'):
                # Parse the current atom number
                current_num = int(line[6:11].strip())
                # Calculate new number
                new_num = current_num + increment
                
                # Reconstruct the line with updated values:
                # First part: Keep original content up to column 54 (coordinates)
                # Then add 1.00 for both occupancy and B-factor
                new_line = (
                    f"{line[:6]}{new_num:5d}{line[11:54]}"  # Original content up to coordinates
                    f"  1.00  1.00"                         # New occupancy and B-factor
                    f"{line[66:]}"                         # Rest of the line (chain ID, etc.)
                )
                outfile.write(new_line)
            else:
                outfile.write(line)

def find_guest_pdb(subfolder_name):
    """Extract guest number and construct guest PDB filename."""
    try:
        # Split the folder name by '_G' and take the last part as the guest number
        guest_number = subfolder_name.split('_G')[-1]
        return f"G{guest_number}.pdb"
    except IndexError:
        print(f"Warning: Folder name {subfolder_name} doesn't match expected pattern HOST_GX")
        return None

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

def get_file_patterns(base_dir):
    """Get file patterns from user based on files in first subfolder."""
    # Get first subfolder
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not subfolders:
        print("Error: No subfolders found")
        return None, None
    
    first_subfolder = os.path.join(base_dir, subfolders[0])
    print(f"\nLooking at files in: {subfolders[0]}")
    
    # Get all PDB files
    pdb_files = [f for f in os.listdir(first_subfolder) if f.endswith('.pdb') 
                 and not f.endswith('_not_docked.pdb')]
    
    # Get all GRO files
    gro_files = [f for f in os.listdir(first_subfolder) if f.endswith('.gro')]
    
    if not pdb_files or not gro_files:
        print("Error: Required files not found in first subfolder")
        return None, None
    
    # Let user select PDB file
    print("\nAvailable PDB files:")
    for i, file in enumerate(pdb_files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Select the PDB file to use as HOST: "))
            if 1 <= choice <= len(pdb_files):
                pdb_host = pdb_files[choice - 1]
                break
            print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")
    
    # Let user select GRO file pattern
    print("\nAvailable GRO files:")
    for i, file in enumerate(gro_files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Select the GRO file to use to generate the plumed.dat file: "))
            if 1 <= choice <= len(gro_files):
                gro_pattern = gro_files[choice - 1]
                break
            print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")
    
    print(f"\nSelected files:")
    print(f"PDB host: {pdb_host}")
    print(f"GRO file: {gro_pattern}")
    
    return pdb_host, gro_pattern

def find_required_files(subfolder_path, pdb_host, gro_pattern):
    """Find the required files using the selected files."""
    pdb_path = os.path.join(subfolder_path, pdb_host)
    gro_path = os.path.join(subfolder_path, gro_pattern)
    
    if not os.path.exists(pdb_path):
        print(f"PDB file not found: {pdb_host}")
        return None, None
    
    if not os.path.exists(gro_path):
        print(f"GRO file not found: {gro_pattern}")
        return None, None
    
    return pdb_path, gro_path

def calculate_center_of_mass(coords, masses):
    total_mass = np.sum(masses)
    center_of_mass = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
    return center_of_mass

def calculate_principal_axes(coords, masses):
    com = calculate_center_of_mass(coords, masses)
    centered_coords = coords - com
    
    inertia_tensor = np.zeros((3, 3))
    for coord, mass in zip(centered_coords, masses):
        inertia_tensor += mass * np.outer(coord, coord)
    
    _, _, Vt = svd(inertia_tensor)
    return com, Vt.T

def calculate_rmax(max_h, r_min, r_molecule):
    closest_rmax = 0
    closest_x = 0
    closest_y = 0

    x_values = [i / 1000.0 for i in range(1, 100001)]

    for x in x_values:
        if x == 0:
            continue

        y = max_h + (r_min / x)

        x_rounded = round(x, 5)
        y_rounded = round(y, 5)

        if abs(x - x_rounded) < 1e-9 and abs(y - y_rounded) < 1e-9:
            rmax = x * y
            rmax_rounded = round(rmax, 5)

            if abs(rmax_rounded - r_molecule) < abs(closest_rmax - r_molecule):
                closest_rmax = rmax_rounded
                closest_x = x_rounded
                closest_y = y_rounded

    return closest_x, closest_y, r_min, closest_rmax, max_h

def calculate_radius(x,y,z):
    if abs(z) <= 0.5:
        return -x * (-y +abs(z))
    else:
        return 0.2

def draw_funnel_and_cylinder(com, principal_axis, max_radius, height, double_funnel=False, other_side=False, vec_step=0.05, angle_sample=36):
    nm_to_angstrom = 10
    max_radius_ang = max_radius * nm_to_angstrom
    height_ang = height * nm_to_angstrom
    
    if other_side:
        principal_axis = -principal_axis
    
    z_values = np.arange(-height, height + vec_step, vec_step) if double_funnel else np.arange(0, height + vec_step, vec_step)
    max_h = 0.5
    r_min = 0.2
    r_molecule = max_radius
    x_formula, y_formula, r_min, rmax, max_h = calculate_rmax(max_h, r_min, r_molecule)
    
    for z in z_values:
        radius = calculate_radius(x_formula, y_formula, z) * nm_to_angstrom
        for angle in np.linspace(0, 2*np.pi, angle_sample, endpoint=False):
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            
            point = np.dot(principal_axis, [x, y, z * nm_to_angstrom])
            pos = com + point
            
            obj_name = 'funnel_cylinder' if z >= 0 or not double_funnel else 'funnel_cylinder_2'
            cmd.pseudoatom(obj_name, pos=pos.tolist())

    cmd.color('orange', 'funnel_cylinder')
    cmd.show_as('nonbonded', 'funnel_cylinder')
    if double_funnel:
        cmd.color('red', 'funnel_cylinder_2')
        cmd.show_as('nonbonded', 'funnel_cylinder_2')
    return x_formula, y_formula, r_min, rmax, max_h

def create_virtual_atoms(com, principal_axis, double_funnel=False, other_side=False):
    nm_to_angstrom = 10
    if double_funnel:
        z_positions = np.arange(-1.25, 1.26, 0.25)
    else:
        z_positions = np.arange(-0.25, 1.76, 0.25)
    
    if other_side:
        principal_axis = -principal_axis
    virtual_atoms = []
    for i, z in enumerate(z_positions):
        point = np.dot(principal_axis, [0, 0, z * nm_to_angstrom])
        pos = com + point
        cmd.pseudoatom(f'virtual_atom_{i}', pos=pos.tolist())
        cmd.show('spheres', f'virtual_atom_{i}')
        cmd.set('sphere_scale', 0.25, f'virtual_atom_{i}')
        cmd.color('yellow', f'virtual_atom_{i}')
        virtual_atoms.append(pos)
    return virtual_atoms

def visualize_molecule(pdb_file, com, principal_axis, max_atom_distance, avg_radius, double_funnel, other_side):
    """Create and save PyMOL visualization without displaying GUI"""
    # Initialize PyMOL in quiet mode without GUI
    pymol.finish_launching(['pymol', '-cq'])

    cmd.load(pdb_file, "molecule")
    cmd.hide("everything", "molecule")
    cmd.show("cartoon", "molecule")
    cmd.color("cyan", "molecule")

    cmd.pseudoatom("com", pos=com.tolist())
    cmd.show("spheres", "com")
    cmd.set("sphere_scale", 0.3, "com")
    cmd.color("red", "com")

    z_axis_length = 12.5 if double_funnel else 17.5
    z_direction = -principal_axis[:, 2] if other_side else principal_axis[:, 2]
    if double_funnel:
        z_start = com - principal_axis[:, 2] * z_axis_length
    else:
        z_start = com - principal_axis[:, 2] * 2.5
    z_end = com + z_direction * z_axis_length
    cmd.pseudoatom("z_start", pos=z_start.tolist())
    cmd.pseudoatom("z_end", pos=z_end.tolist())
    cmd.distance("z_axis", "z_start", "z_end")
    cmd.hide("labels", "z_axis")
    cmd.color("green", "z_axis")

    max_radius_nm = avg_radius / 10
    total_height_nm = 1.5 if double_funnel else 1.8

    print(f"Maximum radius: {max_radius_nm} nm")
    print(f"Total height: {total_height_nm} nm")

    x_formula, y_formula, r_min, rmax, max_h = draw_funnel_and_cylinder(com, principal_axis, max_radius_nm, total_height_nm, double_funnel, other_side)
    virtual_atoms = create_virtual_atoms(com, principal_axis, double_funnel, other_side)

    cmd.label("com", '"Center of Mass"')
    cmd.zoom("all", 1.2)

    # Save session and clean up
    session_file = os.path.splitext(pdb_file)[0] + "_visualization.pse"
    cmd.save(session_file)
    print(f"PyMOL session saved as: {session_file}")
    
    # Clean up PyMOL
    cmd.delete('all')
    cmd.reinitialize()
    
    return virtual_atoms, x_formula, y_formula, r_min, rmax, max_h

def analyze_molecule(pdb_file, double_funnel, other_side):
    parser = PDB.PDBParser()
    structure = parser.get_structure("molecule", pdb_file)
    
    atoms = list(structure.get_atoms())
    coords = np.array([atom.coord for atom in atoms])
    masses = np.array([atom.mass for atom in atoms])
    
    com, principal_axes = calculate_principal_axes(coords, masses)
    
    centered_coords = coords - com
    xy_coords = np.dot(centered_coords, principal_axes[:, :2])
    xy_distances = np.sqrt(np.sum(xy_coords**2, axis=1))
    
    avg_radius = np.mean(xy_distances)
    std_dev = np.std(xy_distances)
    max_atom_distance = np.max(xy_distances)
    
    z_axis_length = 12.5 if double_funnel else 17.5
    z_direction = -principal_axes[:, 2] if other_side else principal_axes[:, 2]
    z_axis_start = com - z_direction * (2.5 if not double_funnel else z_axis_length)
    z_axis_end = com + z_direction * z_axis_length
    
    print(f"Center of Mass: {com}")
    print(f"Principal axis (Z-axis): {z_direction}")
    print(f"Average radius of the ring: {avg_radius:.2f} Angstroms")
    print(f"Standard deviation of radii: {std_dev:.2f} Angstroms")
    print(f"Maximum atom distance from COM in XY plane: {max_atom_distance:.2f} Angstroms")
    print(f"Z-axis start: {z_axis_start}")
    print(f"Z-axis end: {z_axis_end}")
    
    return visualize_molecule(pdb_file, com, principal_axes, max_atom_distance, avg_radius, double_funnel, other_side)

def parse_gro_file(file_path):
    atoms_1mol = []
    first_heavy_atom_1mol = None
    last_heavy_atom_1mol = None
    first_atom_2host = None
    last_atom_2host = None
    first_water_oxygen = None
    last_water_oxygen = None
    water_model = None
    
    molecule_1_pattern = re.compile(r'^1[A-Za-z]+')
    molecule_2_pattern = re.compile(r'^2[A-Za-z]+')
    
    with open(file_path, 'r') as f:
        next(f)  # Skip title
        next(f)  # Skip atom count
        prev_sol_number = None
        atoms_since_last_ow = 0
        for line in f:
            if line.strip() == '':
                break
            parts = line.split()
            if len(parts) >= 6:
                molecule = parts[0]
                atom_label = parts[1]
                atom_number = int(parts[2])
                
                if molecule_1_pattern.match(molecule) and not atom_label.startswith('H'):
                    atoms_1mol.append((atom_label, atom_number, np.array([float(x) for x in parts[3:6]])))
                    if first_heavy_atom_1mol is None:
                        first_heavy_atom_1mol = atom_number
                    last_heavy_atom_1mol = atom_number
                
                elif molecule_2_pattern.match(molecule):
                    if first_atom_2host is None:
                        first_atom_2host = atom_number
                    last_atom_2host = atom_number
                
                elif molecule.endswith('SOL'):
                    current_sol_number = int(molecule[:-3])
                    if atom_label == 'OW':
                        if first_water_oxygen is None:
                            first_water_oxygen = atom_number
                        last_water_oxygen = atom_number
                        
                        if prev_sol_number is not None and current_sol_number != prev_sol_number:
                            if water_model is None:
                                water_model = atoms_since_last_ow
                        
                        atoms_since_last_ow = 0
                        prev_sol_number = current_sol_number
                    else:
                        atoms_since_last_ow += 1
    
    return (atoms_1mol, first_heavy_atom_1mol, last_heavy_atom_1mol, 
            first_atom_2host, last_atom_2host, 
            first_water_oxygen, last_water_oxygen, water_model)

def distance(atom1, atom2):
    return np.linalg.norm(atom1[2] - atom2[2])

def find_furthest_pair(atoms):
    return max(combinations(atoms, 2), key=lambda pair: distance(*pair))

def find_perpendicular_atoms(atoms, furthest_pair):
    atom1, atom2 = furthest_pair
    vector = atom2[2] - atom1[2]
    midpoint = (atom1[2] + atom2[2]) / 2
    
    def find_furthest_perpendicular(point, exclude_atoms):
        return max((atom for atom in atoms if (atom[0], atom[1]) not in exclude_atoms),
                   key=lambda a: np.linalg.norm(np.cross(a[2] - point, vector)))

    exclude_set = {(atom[0], atom[1]) for atom in furthest_pair}
    perp_atom1 = find_furthest_perpendicular(midpoint, exclude_set)
    
    center = (atom1[2] + atom2[2] + perp_atom1[2]) / 3
    exclude_set.add((perp_atom1[0], perp_atom1[1]))
    perp_atom2 = min((atom for atom in atoms if (atom[0], atom[1]) not in exclude_set),
                     key=lambda a: np.linalg.norm(a[2] - center))
    
    return perp_atom1, perp_atom2

def write_plumed_input(pdb_file, gro_file, virtual_atoms, double_funnel, x_formula, y_formula, r_min, rmax, max_h):
    output_file = "plumed.dat"
    
    (atoms_1mol, first_heavy_atom_1mol, last_heavy_atom_1mol, 
     first_atom_2cb8, last_atom_2cb8, 
     first_water_oxygen, last_water_oxygen, water_model) = parse_gro_file(gro_file)
    
    furthest_pair = find_furthest_pair(atoms_1mol)
    perp_atoms = find_perpendicular_atoms(atoms_1mol, furthest_pair)
    
    selected_atoms = furthest_pair + perp_atoms
    first_atom_lig = selected_atoms[0][1]
    last_atom_lig = selected_atoms[-1][1]
    with open(output_file, 'w') as f:
        f.write("# --- (1) ATOMS DEFINITIONS and ALIGNMENT ---\n\n")
        f.write(f"HOST: GROUP ATOMS={first_atom_2cb8}-{last_atom_2cb8}      #host atoms\n")
        f.write(f"LIGC: GROUP ATOMS={first_heavy_atom_1mol}-{last_heavy_atom_1mol}  #heavy atoms in the ligand\n")
        for i, atom in enumerate(selected_atoms, start=1):
            f.write(f"l{i}: GROUP ATOMS={atom[1]}             #ligand selected atoms\n")
        f.write(f"WO: GROUP ATOMS={first_water_oxygen}-{last_water_oxygen}:{water_model+1}    #water oxygen atoms\n\n")

        f.write("WHOLEMOLECULES ENTITY0=HOST\n")
        f.write(f"FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL #coordinates alignment\n")
        f.write("lig: CENTER ATOMS=LIGC\n")
        
        # Write virtual atoms
        for i, atom in enumerate(virtual_atoms, start=1):
            f.write(f"v{i}: FIXEDATOM AT={abs(atom[0]/10):.4f},{abs(atom[1]/10):.4f},{atom[2]/10:.4f}\n")
        
        f.write("cyl: DISTANCE ATOMS=v1,lig COMPONENTS\n")
        f.write("radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO\n")
        
        f.write("\n# --- (2) DESCRIPTORS ---\n")
        for i, atom in enumerate(selected_atoms, start=1):
            f.write(f"L{i}: COORDINATION GROUPA=l{i} GROUPB=WO SWITCH={{RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8}} NLIST NL_CUTOFF=1.5 NL_STRIDE=20\n")
        
        for i in range(1, len(virtual_atoms) + 1):
            f.write(f"V{i}: COORDINATION GROUPA=v{i} GROUPB=WO SWITCH={{RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8}} NLIST NL_CUTOFF=1.5 NL_STRIDE=20\n")
        
        for i in range(1, 5):
            f.write(f"d{i}: MATHEVAL ARG=L{i} FUNC=(x/2.5)-1.0 PERIODIC=NO\n")
        
        for i in range(5, len(virtual_atoms) + 5):
            f.write(f"d{i}: MATHEVAL ARG=V{i-4} FUNC=(x/2.8)-1.0 PERIODIC=NO\n")

        f.write("\n# --- (3) FUNNEL AND WALLS DEFINITION ---\n")
        if double_funnel:
            f.write(f"funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+{x_formula}*(-{y_formula}+z))*step(-z+{max_h})+(r-{r_min})*step(z-{max_h}))*step(z-0.0)+((r-{x_formula}*({y_formula}+z))*step(z+{max_h})+(r-{r_min})*step(-z-{max_h}))*step(-z+0.0) PERIODIC=NO\n")
        else:
            f.write(f"funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+{x_formula}*(-{y_formula}+z))*step(-z+{max_h})+(r-{r_min})*step(z-{max_h}))*step(z-0.0) PERIODIC=NO\n")
        f.write("UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall  #funnel restraint\n")
        f.write("UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall\n")
        if double_funnel:
            f.write("LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall\n")
        f.write(f"ang: ANGLE ATOMS=v3,v5,{first_atom_lig},{last_atom_lig}   #angle of a ligand's axis with z\n")
        f.write("cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO\n")
        f.write("ene: ENERGY\n")
        
        f.write("\n# --- (4) OPES  ---\n")
        f.write("\nPRINT ARG=cyl.z,radius,cosang,")
        for i, atom in enumerate(selected_atoms, start=1):
            f.write(f"L{i},")
        for i, atom in enumerate(virtual_atoms, start=1):
            if i != len(virtual_atoms):
                f.write(f"V{i},")
            else:
                f.write(f"V{i}")
        f.write(" STRIDE=50 FILE=COLVAR FMT=%8.4f\n\n")
        f.write("ENDPLUMED")

    print(f"PLUMED input file saved as: {output_file}")

def process_subfolder(subfolder_path, pdb_host, gro_pattern, double_funnel, other_side):
    """Process a single subfolder."""
    print(f"\nProcessing subfolder: {os.path.basename(subfolder_path)}")
    
    # Find required files
    pdb_file, gro_file = find_required_files(subfolder_path, pdb_host, gro_pattern)
    
    if not pdb_file or not gro_file:
        print(f"Skipping {os.path.basename(subfolder_path)}: Required files not found")
        return False
    
    # Store current directory and change to working directory
    original_dir = os.getcwd()
    os.chdir(subfolder_path)
    
    try:
        # Process the files using original functionality
        virtual_atoms, x_formula, y_formula, r_min, rmax, max_h = analyze_molecule(
            pdb_file, double_funnel, other_side)
        write_plumed_input(pdb_file, gro_file, virtual_atoms, double_funnel,
                          x_formula, y_formula, r_min, rmax, max_h)
        
        # Additional renumbering functionality
        subfolder_name = os.path.basename(subfolder_path)
        guest_pdb = find_guest_pdb(subfolder_name)
        if guest_pdb and os.path.exists(os.path.join(subfolder_path, guest_pdb)):
            print(f"\nGenerating template PDB file...")
            last_atom_number = get_last_atom_number(guest_pdb)
            template_pdb = "host_template.pdb"
            update_atom_numbers(pdb_file, template_pdb, last_atom_number)
            print(f"Created template PDB file: {template_pdb}")
        else:
            print(f"\nSkipping renumbering: Guest PDB file not found")
        
        print(f"Processing completed successfully for {os.path.basename(subfolder_path)}")
        return True
        
    except Exception as e:
        print(f"Error processing {os.path.basename(subfolder_path)}: {str(e)}")
        print(f"Full error details:")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # Return to original directory
        os.chdir(original_dir)

def main():
    parser = argparse.ArgumentParser(description="Analyze molecules and create PLUMED input files for all systems")
    parser.add_argument("--double_funnel", action="store_true", 
                       help="Create a double funnel (default is single funnel)")
    parser.add_argument("--other_side", action="store_true", 
                       help="Flip the direction of the Z-axis and funnel")
    parser.add_argument("--skip_renumber", action="store_true",
                       help="Skip the generation of template PDB files")
    args = parser.parse_args()

    # Get the base directory path
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    
    # Select the input folder
    input_folder = select_input_folder()
    if not input_folder:
        print("No valid input folder selected. Exiting.")
        return
    
    base_dir = os.path.join(prep_dir, input_folder)
    print(f"\nProcessing all subfolders in: {base_dir}")
    
    # Get file patterns from user
    pdb_host, gro_pattern = get_file_patterns(base_dir)
    if not pdb_host or not gro_pattern:
        print("Failed to get file patterns. Exiting.")
        return
    
    # Process all subfolders
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not subfolders:
        print("Error: No subfolders found")
        return
    
    successful = 0
    failed = 0
    
    for subfolder in sorted(subfolders):
        subfolder_path = os.path.join(base_dir, subfolder)
        if process_subfolder(subfolder_path, pdb_host, gro_pattern, 
                           args.double_funnel, args.other_side):
            successful += 1
        else:
            failed += 1
    
    print("\nProcessing Summary:")
    print(f"Total subfolders processed: {len(subfolders)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")

if __name__ == "__main__":
    main()
