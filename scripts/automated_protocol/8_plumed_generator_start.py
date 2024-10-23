import numpy as np
from itertools import combinations

def parse_gro_file(file_path):
    atoms_1mol = []
    first_heavy_atom_1mol = None
    last_heavy_atom_1mol = None
    first_atom_2cb8 = None
    last_atom_2cb8 = None
    first_water_oxygen = None
    last_water_oxygen = None
    water_model = None
    
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
                
                if molecule == '1MOL' and not atom_label.startswith('H'):
                    atoms_1mol.append((atom_label, atom_number, np.array([float(x) for x in parts[3:6]])))
                    if first_heavy_atom_1mol is None:
                        first_heavy_atom_1mol = atom_number
                    last_heavy_atom_1mol = atom_number
                
                elif molecule == '2CB8':
                    if first_atom_2cb8 is None:
                        first_atom_2cb8 = atom_number
                    last_atom_2cb8 = atom_number
                
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
            first_atom_2cb8, last_atom_2cb8, 
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

def main(file_path):
    (atoms_1mol, first_heavy_atom_1mol, last_heavy_atom_1mol, 
     first_atom_2cb8, last_atom_2cb8, 
     first_water_oxygen, last_water_oxygen, water_model) = parse_gro_file(file_path)
    
    furthest_pair = find_furthest_pair(atoms_1mol)
    perp_atoms = find_perpendicular_atoms(atoms_1mol, furthest_pair)
    
    selected_atoms = furthest_pair + perp_atoms
    
    with open("plumed_info.txt", "w") as f:
        for i, atom in enumerate(selected_atoms, 1):
            f.write(f"l{i}: GROUP ATOMS={atom[1]}             #ligand selected atoms\n")
        f.write(f"\nLIGC: GROUP ATOMS={first_heavy_atom_1mol}-{last_heavy_atom_1mol}  #heavy atoms in the ligand\n")
        f.write(f"\nHOST: GROUP ATOMS={first_atom_2cb8}-{last_atom_2cb8}      #host atoms\n")
        f.write(f"\nWO: GROUP ATOMS={first_water_oxygen}-{last_water_oxygen}:{water_model+1}    #water oxygen atoms\n")

if __name__ == "__main__":
    file_path = "npt.gro"
    main(file_path)
