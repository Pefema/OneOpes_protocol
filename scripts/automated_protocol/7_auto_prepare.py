import os
import sys
import subprocess
import argparse

def run_command(command, error_message, inputs=None):
    print(f"\n{'#'*20} COMMAND START {'#'*20}")
    print(command if isinstance(command, str) else ' '.join(command))
    print(f"{'#'*20} COMMAND END {'#'*20}\n")
    
    try:
        if inputs:
            if isinstance(command, str):
                process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            else:
                process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate(input=inputs)
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, command, stdout, stderr)
        else:
            if isinstance(command, str):
                result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            else:
                result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = result.stdout, result.stderr
        
        print(stdout)
        if stderr:
            print(stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error: {error_message}")
        print(f"Command output: {e.output}")
        sys.exit(1)

def convert_pdb_to_gro(file_path):
    if file_path.endswith('.pdb'):
        gro_file = file_path.rsplit('.', 1)[0] + '.gro'
        run_command(f"gmx_mpi editconf -f {file_path} -o {gro_file}",
                    f"Failed to convert {file_path} to .gro format")
        return gro_file
    return file_path

def main():
    parser = argparse.ArgumentParser(description="Setup Gromacs simulation with optional water model selection.")
    parser.add_argument("protein_file", help="Protein .gro or .pdb file")
    parser.add_argument("ligand_file", help="Ligand .gro or .pdb file")
    parser.add_argument("topol_file", help="Topology .top file")
    parser.add_argument("--water_points", type=int, choices=[3, 4, 5], default=3, 
                        help="Number of points in water model (default: 3)")
    
    args = parser.parse_args()

    # Convert PDB files to GRO if necessary
    protein_file = convert_pdb_to_gro(args.protein_file)
    ligand_file = convert_pdb_to_gro(args.ligand_file)
    topol_file = args.topol_file
    water_points = args.water_points

    # Get the number of atoms in each file
    with open(ligand_file, 'r') as f:
        n_atoms_ligand = int(f.readlines()[1])
    with open(protein_file, 'r') as f:
        n_atoms_protein = int(f.readlines()[1])

    # Calculate the total number of atoms
    total_atoms = n_atoms_ligand + n_atoms_protein

    # Create complex_file.gro
    with open('complex_file.gro', 'w') as complex_file:
        with open(ligand_file, 'r') as f:
            complex_file.write(f.readline())  # Write header
        complex_file.write(f"{total_atoms}\n")  # Write total atoms
        
        # Write atoms from ligand file
        with open(ligand_file, 'r') as f:
            complex_file.writelines(f.readlines()[2:n_atoms_ligand+2])
        
        # Write atoms from protein file
        with open(protein_file, 'r') as f:
            complex_file.writelines(f.readlines()[2:n_atoms_protein+2])
        
        # Write box vectors from protein file
        with open(protein_file, 'r') as f:
            complex_file.write(f.readlines()[-1])

    # Copy topol.top
    os.system(f"cp {topol_file} topol.top")

    # Set OMP_NUM_THREADS
    os.environ['OMP_NUM_THREADS'] = '2'
    print("export OMP_NUM_THREADS=2")

    # Run gmx_mpi commands
    
    # Create a new box
    print("Creating a new box...")
    run_command("gmx_mpi editconf -f complex_file.gro -o newbox.gro -c -d 1.0 -bt cubic",
                "Failed to create a new box")

    # Solvate the system
    print("Solvating the system...")
    water_model = {3: "spc216.gro", 4: "tip4p.gro", 5: "tip5p.gro"}[water_points]
    run_command(f"gmx_mpi solvate -cp newbox.gro -cs {water_model} -o solv.gro -p topol.top",
                "Failed to solvate the system")

    # Prepare for adding ions
    print("Preparing for adding ions...")
    run_command("gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2",
                "Failed to prepare for adding ions")

    # Add ions to neutralize the system
    print("Adding ions to neutralize the system...")
    run_command("gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname Na+ -nname Cl- -neutral -conc 0.25",
                "Failed to add ions",
                inputs="5\n")

    # Correct ion names in solv_ions.gro
    print("Correcting ion names...")
    run_command("sed -i 's/ Na /Na+ /g' solv_ions.gro", "Failed to correct Na+ names")
    run_command("sed -i 's/ Cl /Cl- /g' solv_ions.gro", "Failed to correct Cl- names")

    # Energy minimization
    print("Running energy minimization...")
    run_command("gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2",
                "Failed to prepare for energy minimization")
    run_command("gmx_mpi mdrun -v -deffnm em",
                "Energy minimization failed")

    # Generate index files and position restraints
    print("Generating index files and position restraints...")
    run_command(['gmx_mpi', 'make_ndx', '-f', ligand_file, '-o', 'index_ligand.ndx'],
                "Failed to generate ligand index file",
                inputs="0 & ! a H*\nq\n")
    
    run_command(['gmx_mpi', 'genrestr', '-f', ligand_file, '-n', 'index_ligand.ndx', '-o', 'posre_ligand.itp', '-fc', '1000', '1000', '1000'],
                "Failed to generate ligand position restraints",
                inputs="3\n")
    
    run_command(['gmx_mpi', 'make_ndx', '-f', protein_file, '-o', 'index_protein.ndx'],
                "Failed to generate protein index file",
                inputs="0 & ! a H*\nq\n")
    
    run_command(['gmx_mpi', 'genrestr', '-f', protein_file, '-n', 'index_protein.ndx', '-o', 'posre.itp', '-fc', '1000', '1000', '1000'],
                "Failed to generate protein position restraints",
                inputs="3\n")
    
    run_command(['gmx_mpi', 'make_ndx', '-f', 'em.gro', '-o', 'index.ndx'],
                "Failed to generate combined index file",
                inputs="3|2\n4|5|7\nq\n")

    # NVT equilibration
    print("Running NVT equilibration...")
    run_command("gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2",
                "Failed to prepare for NVT equilibration")
    run_command("gmx_mpi mdrun -deffnm nvt",
                "NVT equilibration failed")

    # NPT equilibration
    print("Running NPT equilibration...")
    run_command("gmx_mpi grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2",
                "Failed to prepare for NPT equilibration")
    run_command("gmx_mpi mdrun -deffnm npt",
                "NPT equilibration failed")

    # Prepare final run
    print("Preparing for production run...")
    run_command("gmx_mpi grompp -f NVT.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o topol.tpr -maxwarn 2",
                "Failed to prepare for production run")

    print(f"Setup completed successfully with {water_points}-point water model!")

if __name__ == "__main__":
    main()
