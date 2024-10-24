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
                process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                        text=True, bufsize=1)
            else:
                process = subprocess.Popen(command, stdin=subprocess.PIPE, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                        text=True, bufsize=1)
            
            stdout, stderr = process.communicate(input=inputs)
            print(stdout)
            if stderr:
                print("STDERR output:")
                print(stderr)
                
            if process.returncode != 0:
                print("\nCommand failed with return code:", process.returncode)
                print("\nFull command output:")
                print("STDOUT:")
                print(stdout)
                print("\nSTDERR:")
                print(stderr)
                raise subprocess.CalledProcessError(process.returncode, command, stdout, stderr)
        else:
            if isinstance(command, str):
                process = subprocess.Popen(command, shell=True, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                        text=True, bufsize=1)
            else:
                process = subprocess.Popen(command, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                        text=True, bufsize=1)
            
            # Real-time output processing
            while True:
                stdout_line = process.stdout.readline()
                stderr_line = process.stderr.readline()
                
                if stdout_line:
                    print(stdout_line.rstrip())
                if stderr_line:
                    print("STDERR:", stderr_line.rstrip())
                
                if not stdout_line and not stderr_line:
                    break
            
            process.wait()
            if process.returncode != 0:
                print("\nCommand failed with return code:", process.returncode)
                # Get any remaining output
                remaining_stdout, remaining_stderr = process.communicate()
                if remaining_stdout:
                    print("\nRemaining STDOUT:")
                    print(remaining_stdout)
                if remaining_stderr:
                    print("\nRemaining STDERR:")
                    print(remaining_stderr)
                raise subprocess.CalledProcessError(process.returncode, command)
            
    except subprocess.CalledProcessError as e:
        print(f"\nError: {error_message}")
        if hasattr(e, 'output') and e.output:
            print("\nDetailed error output:")
            if hasattr(e, 'stdout') and e.stdout:
                print("\nSTDOUT:")
                print(e.stdout)
            if hasattr(e, 'stderr') and e.stderr:
                print("\nSTDERR:")
                print(e.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {str(e)}")
        sys.exit(1)

def convert_pdb_to_gro(file_path):
    if file_path.endswith('.pdb'):
        gro_file = file_path.rsplit('.', 1)[0] + '.gro'
        run_command(f"gmx_mpi editconf -f {file_path} -o {gro_file}",
                    f"Failed to convert {file_path} to .gro format")
        return gro_file
    return file_path

def select_input_folder():
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    folders = [d for d in os.listdir(prep_dir)
               if os.path.isdir(os.path.join(prep_dir, d))]
    
    if not folders:
        print("Error: No folders found in system_preparation directory.")
        sys.exit(1)
    
    print("\nAvailable input folders:")
    for i, folder in enumerate(folders, 1):
        print(f"{i}. {folder}")
    
    while True:
        try:
            choice = int(input("Enter the number of the folder you want to process: "))
            if 1 <= choice <= len(folders):
                return folders[choice - 1]
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def select_host_pdb(subfolder_path):
    pdb_files = [f for f in os.listdir(subfolder_path) if f.endswith('.pdb')]
    
    if not pdb_files:
        print(f"Error: No PDB files found in {subfolder_path}")
        sys.exit(1)
    
    print("\nAvailable PDB files for host:")
    for i, file in enumerate(pdb_files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Select the host PDB file to use for all preparations: "))
            if 1 <= choice <= len(pdb_files):
                return pdb_files[choice - 1]
            print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def select_topology_file(subfolder_path):
    top_files = [f for f in os.listdir(subfolder_path) if f.endswith('.top')]
    
    if not top_files:
        print(f"Error: No topology files found in {subfolder_path}")
        sys.exit(1)
    
    print("\nAvailable topology files:")
    for i, file in enumerate(top_files, 1):
        print(f"{i}. {file}")
    
    while True:
        try:
            choice = int(input("Select the topology file to use for all preparations: "))
            if 1 <= choice <= len(top_files):
                return top_files[choice - 1]
            print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def process_subfolders(base_folder, water_points):
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    base_dir = os.path.join(prep_dir, base_folder)
    
    # Get all subfolders
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not subfolders:
        print("Error: No subfolders found")
        return
    
    # Select host PDB and topology files from the first subfolder
    first_subfolder = os.path.join(base_dir, subfolders[0])
    host_basename = select_host_pdb(first_subfolder)
    original_topol = select_topology_file(first_subfolder)
    
    # Process each subfolder
    for subfolder in subfolders:
        subfolder_path = os.path.join(base_dir, subfolder)
        print(f"\nProcessing subfolder: {subfolder}")
        
        # Extract guest number from folder name (assuming format CB8_GX)
        try:
            guest_number = subfolder.split('_G')[-1]
            guest_file = f"G{guest_number}.pdb"
        except IndexError:
            print(f"Warning: Folder name {subfolder} doesn't match expected pattern CB8_GX")
            continue
        
        # Check if required files exist
        host_path = os.path.join(subfolder_path, host_basename)
        guest_path = os.path.join(subfolder_path, guest_file)
        original_topol_path = os.path.join(subfolder_path, original_topol)
        
        missing_files = []
        if not os.path.exists(host_path):
            missing_files.append(host_basename)
        if not os.path.exists(guest_path):
            missing_files.append(guest_file)
        if not os.path.exists(original_topol_path):
            missing_files.append(original_topol)
            
        if missing_files:
            print(f"Missing files in {subfolder}: {', '.join(missing_files)}. Skipping...")
            continue
        
        # Change to subfolder directory
        os.chdir(subfolder_path)
        
        # Copy the selected topology file to topol.top
        print(f"Copying {original_topol} to topol.top")
        run_command(f"cp {original_topol} topol.top", "Failed to copy topology file")
        
        print(f"Processing system with:")
        print(f"Host: {host_basename}")
        print(f"Guest: {guest_file}")
        print(f"Topology: topol.top (copied from {original_topol})")
        
        # Convert PDB files to GRO if needed
        host_gro = convert_pdb_to_gro(host_path)
        guest_gro = convert_pdb_to_gro(guest_path)
        
        # Get the number of atoms in each file
        with open(guest_gro, 'r') as f:
            n_atoms_ligand = int(f.readlines()[1])
        with open(host_gro, 'r') as f:
            n_atoms_protein = int(f.readlines()[1])
        
        # Create complex_file.gro
        with open('complex_file.gro', 'w') as complex_file:
            with open(guest_gro, 'r') as f:
                complex_file.write(f.readline())  # Write header
            complex_file.write(f"{n_atoms_ligand + n_atoms_protein}\n")  # Write total atoms
            
            # Write atoms from ligand file
            with open(guest_gro, 'r') as f:
                complex_file.writelines(f.readlines()[2:n_atoms_ligand+2])
            
            # Write atoms from protein file
            with open(host_gro, 'r') as f:
                complex_file.writelines(f.readlines()[2:n_atoms_protein+2])
            
            # Write box vectors from protein file
            with open(host_gro, 'r') as f:
                complex_file.write(f.readlines()[-1])
        
        # Run the preparation commands
        os.environ['OMP_NUM_THREADS'] = '2'
        print("export OMP_NUM_THREADS=2")

        print("Creating a new box...")
        run_command("gmx_mpi editconf -f complex_file.gro -o newbox.gro -c -d 1.0 -bt cubic",
                    "Failed to create a new box")

        print("Solvating the system...")
        water_model = {3: "spc216.gro", 4: "tip4p.gro", 5: "tip5p.gro"}[water_points]
        run_command(f"gmx_mpi solvate -cp newbox.gro -cs {water_model} -o solv.gro -p topol.top",
                    "Failed to solvate the system")

        print("Preparing for adding ions...")
        run_command("gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2",
                    "Failed to prepare for adding ions")

        print("Adding ions to neutralize the system...")
        run_command("gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname Na+ -nname Cl- -neutral -conc 0.25",
                    "Failed to add ions",
                    inputs="5\n")

        print("Correcting ion names...")
        run_command("sed -i 's/ Na /Na+ /g' solv_ions.gro", "Failed to correct Na+ names")
        run_command("sed -i 's/ Cl /Cl- /g' solv_ions.gro", "Failed to correct Cl- names")

        print("Running energy minimization...")
        run_command("gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2",
                    "Failed to prepare for energy minimization")
        run_command("gmx_mpi mdrun -v -deffnm em",
                    "Energy minimization failed")

        print("Generating index files and position restraints...")
        run_command(['gmx_mpi', 'make_ndx', '-f', guest_gro, '-o', 'index_ligand.ndx'],
                    "Failed to generate ligand index file",
                    inputs="0 & ! a H*\nq\n")
        
        run_command(['gmx_mpi', 'genrestr', '-f', guest_gro, '-n', 'index_ligand.ndx', '-o', 'posre_ligand.itp', '-fc', '1000', '1000', '1000'],
                    "Failed to generate ligand position restraints",
                    inputs="3\n")
        
        run_command(['gmx_mpi', 'make_ndx', '-f', host_gro, '-o', 'index_protein.ndx'],
                    "Failed to generate protein index file",
                    inputs="0 & ! a H*\nq\n")
        
        run_command(['gmx_mpi', 'genrestr', '-f', host_gro, '-n', 'index_protein.ndx', '-o', 'posre.itp', '-fc', '1000', '1000', '1000'],
                    "Failed to generate protein position restraints",
                    inputs="3\n")
        
        run_command(['gmx_mpi', 'make_ndx', '-f', 'em.gro', '-o', 'index.ndx'],
                    "Failed to generate combined index file",
                    inputs="3|2\n4|5|7\nq\n")

        print("Running NVT equilibration...")
        run_command("gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2",
                    "Failed to prepare for NVT equilibration")
        run_command("gmx_mpi mdrun -deffnm nvt",
                    "NVT equilibration failed")

        print("Running NPT equilibration...")
        run_command("gmx_mpi grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2",
                    "Failed to prepare for NPT equilibration")
        run_command("gmx_mpi mdrun -deffnm npt",
                    "NPT equilibration failed")

        print("Preparing for production run...")
        run_command("gmx_mpi grompp -f NVT.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o topol.tpr -maxwarn 2",
                    "Failed to prepare for production run")

        print(f"Setup completed successfully for {subfolder}!")

def main():
    parser = argparse.ArgumentParser(description="Setup Gromacs simulation with optional water model selection.")
    parser.add_argument("--water_points", type=int, choices=[3, 4, 5], default=3, 
                        help="Number of points in water model (default: 3)")
    
    args = parser.parse_args()
    
    # Select the input folder
    input_folder = select_input_folder()
    if input_folder:
        process_subfolders(input_folder, args.water_points)
    else:
        print("No valid input folder selected. Exiting.")
        return

if __name__ == "__main__":
    main()
