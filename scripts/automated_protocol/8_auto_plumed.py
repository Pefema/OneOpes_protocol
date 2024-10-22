import sys
import numpy as np
from Bio import PDB
import argparse
import os
import pymol
from pymol import cmd, stored
from scipy.linalg import svd
from itertools import combinations
import re

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
    pymol.finish_launching()

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

    session_file = os.path.splitext(pdb_file)[0] + "_visualization.pse"
    cmd.save(session_file)
    print(f"PyMOL session saved as: {session_file}")
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
        f.write(f"FIT_TO_TEMPLATE STRIDE=1 REFERENCE={os.path.basename(pdb_file)} TYPE=OPTIMAL #coordinates alignment\n")
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
def main():
    parser = argparse.ArgumentParser(description="Analyze a molecule from a PDB file and create a PLUMED input file")
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("gro_file", help="Path to the GRO file")
    parser.add_argument("--double_funnel", action="store_true", help="Create a double funnel (default is single funnel)")
    parser.add_argument("--other_side", action="store_true", help="Flip the direction of the Z-axis and funnel")
    args = parser.parse_args()

    virtual_atoms, x_formula, y_formula, r_min, rmax, max_h = analyze_molecule(args.pdb_file, args.double_funnel, args.other_side)
    write_plumed_input(args.pdb_file, args.gro_file, virtual_atoms, args.double_funnel, x_formula, y_formula, r_min, rmax, max_h)

if __name__ == "__main__":
    main()
