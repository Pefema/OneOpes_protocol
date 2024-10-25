import os
import sys
import subprocess
import argparse
import re
import shutil
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from Bio import PDB

@dataclass
class ReplicaSetup:
    """Container for replica setup data"""
    original_content: str
    sigma_values: Dict[str, float]
    reference_pdb: str
    base_path: str
    folder_name: str

def run_command(command: str, error_message: str, inputs: Optional[str] = None) -> None:
    """Execute a command and handle its output."""
    print(f"\n{'#'*20} COMMAND START {'#'*20}")
    print(command)
    print(f"{'#'*20} COMMAND END {'#'*20}\n")
    
    try:
        if inputs:
            process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, 
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                    text=True, bufsize=1)
            stdout, stderr = process.communicate(input=inputs)
            print(stdout)
            if stderr:
                print("STDERR output:")
                print(stderr)
                
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, command, stdout, stderr)
        else:
            process = subprocess.Popen(command, shell=True, 
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                    text=True, bufsize=1)
            
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
        raise RuntimeError(f"Command failed: {error_message}")
    except Exception as e:
        print(f"\nUnexpected error: {str(e)}")
        raise

def select_input_folder() -> str:
    """Select input folder from system_preparation directory."""
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    folders = [d for d in os.listdir(prep_dir)
               if os.path.isdir(os.path.join(prep_dir, d))]
    
    if not folders:
        raise RuntimeError("No folders found in system_preparation directory.")
    
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

def check_required_files(subfolder_path: str) -> List[str]:
    """Check if required files exist in the subfolder."""
    required_files = ['topol.tpr', 'plumed.dat']
    return [file for file in required_files 
            if not os.path.exists(os.path.join(subfolder_path, file))]

def run_production(subfolder_path: str, nsteps: int, ncores: int) -> bool:
    """Run the production simulation in the specified subfolder."""
    missing_files = check_required_files(subfolder_path)
    if missing_files:
        print(f"Error: Missing required files in {subfolder_path}:")
        for file in missing_files:
            print(f"  - {file}")
        return False
    
    original_dir = os.getcwd()
    os.chdir(subfolder_path)
    
    try:
        os.environ['OMP_NUM_THREADS'] = '2'
        command = f"mpirun -n {ncores} gmx_mpi mdrun -deffnm topol -plumed plumed.dat -pin on -pinoffset 0 --nsteps {nsteps}"
        run_command(command, "Production run failed")
        return True
        
    except Exception as e:
        print(f"Error running production in {os.path.basename(subfolder_path)}: {str(e)}")
        return False
        
    finally:
        os.chdir(original_dir)

def calculate_std_dev(colvar_file: str) -> Dict[str, float]:
    """Calculate standard deviation for each column in COLVAR file."""
    expected_columns = ["cyl.z", "L1", "L2", "L3", "L4", "V1", "V2", "V3", "V4", "V5", 
                       "V6", "V7", "V8", "V9", "V10", "V11", "cosang"]
    
    # Read header to get column indices
    with open(colvar_file, 'r') as f:
        header = next(f).replace('#! FIELDS ', '').strip().split()
        indices = {name: i for i, name in enumerate(header)}
    
    # Read data and calculate standard deviations
    data = np.loadtxt(colvar_file, comments='#')
    std_devs = {}
    
    missing_columns = []
    for name in expected_columns:
        if name in indices:
            col_data = data[:, indices[name]]
            std_dev = np.std(col_data)
            std_devs[name] = round(float(std_dev), 2)
        else:
            missing_columns.append(name)
            # Provide a default value for missing columns
            std_devs[name] = 0.5  # Default value
            print(f"Warning: Column {name} not found in COLVAR file, using default value of 0.5")
    
    if missing_columns:
        print(f"Warning: The following columns were missing from COLVAR: {', '.join(missing_columns)}")
    
    return std_devs

def write_std_results(std_devs: Dict[str, float], output_file: str) -> None:
    """Write standard deviation results to file."""
    with open(output_file, 'w') as f:
        f.write("Column\tStandard Deviation\n")
        for col, std in std_devs.items():
            f.write(f"{col}\t{std}\n")

def get_reference_pdb(plumed_content: str) -> Optional[str]:
    """Extract reference PDB filename from PLUMED content."""
    match = re.search(r'REFERENCE=(.*?\.pdb)', plumed_content)
    return match.group(1) if match else None

def copy_additional_files(source_folder: str, destination_folder: str, files_to_copy: List[str]) -> None:
    """Copy required files to the replica folders."""
    # Look for files in the subfolder itself, not in the plumed directory
    for file in files_to_copy:
        source_file = os.path.join(source_folder, file)
        if os.path.exists(source_file):
            shutil.copy2(source_file, destination_folder)
        else:
            print(f"Warning: {file} not found in {source_folder}")
def generate_sbatch_script(folder_path: str, folder_name: str) -> None:
    """Generate the SLURM batch script for the replica system."""
    sbatch_content = f"""#!/bin/sh
#SBATCH --partition=private-gervasio-gpu
#SBATCH --time 48:00:00
#SBATCH --gpus=1
#SBATCH --constraint="V7|V8"
#SBATCH --job-name {folder_name.split('_', 1)[1] if '_' in folder_name else folder_name}
#SBATCH --error jobname-error.e%j
#SBATCH --output jobname-out.o%j
#SBATCH --ntasks 8
#SBATCH --cpus-per-task 2
#SBATCH --nodes 1
export OMP_NUM_THREADS=2
alias splumedhome='module load GCC/12.2.0 && module load OpenMPI && module load CUDA/12.0.0 && module load Python && source /home/users/f/febrerma/plumed2-2.9/sourceme.sh; source /home/users/f/febrerma/gromacs-2023/install_single_cuda/bin/GMXRC'
splumedhome
srun gmx_mpi mdrun -deffnm topol -multidir 0 1 2 3 4 5 6 7 -replex 1000 -hrex -plumed plumed.dat -nsteps 100000000 -notunepme
"""
    
    with open(os.path.join(folder_path, "run_simulation.sh"), 'w') as f:
        f.write(sbatch_content)

def create_base_plumed_content(setup: ReplicaSetup) -> str:
    """Create the base PLUMED content with initial modifications."""
    content = setup.original_content
    
    # Add base OPES_METAD_EXPLORE
    base_content = insert_content(content, f"""OPES_METAD_EXPLORE ...
   LABEL=opes
   ARG=cyl.z,cosang
   SIGMA={setup.sigma_values['cyl.z']:.3f},{setup.sigma_values['cosang']:.3f}
   FILE=Kernels.data
   STATE_RFILE=compressed.Kernels.data
   STATE_WFILE=compressed.Kernels.data
   PACE=10000
   BARRIER=100
... OPES_METAD_EXPLORE""")
    
    # Update print line
    return update_print_line(base_content, 
        "opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11")

def insert_content(content: str, new_content: str) -> str:
    """Insert new content into PLUMED file at appropriate location."""
    lines = content.split('\n')
    print_line_index = next(i for i, line in enumerate(lines) if line.startswith('PRINT'))
    endplumed_line_index = next(i for i, line in enumerate(lines) if line.startswith('ENDPLUMED'))
    
    return '\n'.join(lines[:print_line_index] + [new_content] + 
                    lines[print_line_index:endplumed_line_index] + 
                    lines[endplumed_line_index:])

def update_print_line(content: str, new_args: str) -> str:
    """Update the PRINT line with new arguments."""
    lines = content.split('\n')
    print_line_index = next(i for i, line in enumerate(lines) if line.startswith('PRINT'))
    lines[print_line_index] = re.sub(r'ARG=.*?(?=\s)', f'ARG={new_args}', lines[print_line_index])
    return '\n'.join(lines)

def insert_ecv_opes(content: str, temp: int) -> str:
    """Insert ECV_MULTITHERMAL and OPES_EXPANDED sections."""
    lines = content.split('\n')
    opes_index = next(i for i, line in enumerate(lines) if line.startswith('# --- (4) OPES  ---'))
    
    ecv_opes_lines = f"\necv: ECV_MULTITHERMAL ARG=ene TEMP_MAX={temp}\nopesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=100\n"
    
    lines.insert(opes_index + 1, ecv_opes_lines)
    return '\n'.join(lines)

def update_temp_max(content: str, temp: int) -> str:
    """Update the temperature maximum in the ECV_MULTITHERMAL section."""
    return re.sub(r'TEMP_MAX=\d+', f'TEMP_MAX={temp}', content)

def setup_replica(setup: ReplicaSetup, replica_num: int, content: str, 
                 temperature: Optional[int] = None) -> str:
    if replica_num == 0:
        return content
    
    # Maintain a list of all previous modifications
    modifications = []
    for i in range(1, replica_num + 1):
        cv_mappings = {
            1: ('L4', 'L4'),
            2: ('V6', 'V6'),
            3: ('L1', 'L1'),
            4: ('V8', 'V8'),
            5: ('V4', 'V4'),
            6: ('V10', 'V10'),
            7: ('V2', 'V2')
        }
        
        if i in cv_mappings:
            cv, label = cv_mappings[i]
            new_content = f"""OPES_METAD_EXPLORE ...
   LABEL=opese{i}
   ARG={cv}
   SIGMA={setup.sigma_values[cv]:.3f}
   FILE=Kernels{i}.data
   STATE_RFILE=compressed_Kernels{i}.data
   STATE_WFILE=compressed.Kernels{i}.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE"""
            modifications.append(new_content)
            
            # Build bias list for PRINT line
            bias_list = [f"opese{j}.bias" for j in range(1, i + 1)]
            if i >= 4:
                bias_list.insert(0, "opesX.bias")
            
            base_args = "opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias"
            if i >= 3:
                base_args += ",ene"
            base_args += ",cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11"
            
            content = update_print_line(content, f"{base_args},{','.join(bias_list)}")
    
    # Apply all modifications
    content = content.replace("STATE_RFILE=compressed.Kernels.data", 
                            "STATE_RFILE=compressed_Kernels.data")
    for mod in modifications:
        content = insert_content(content, mod)
    
    # Add temperature-dependent modifications
    if temperature and replica_num >= 4:
        if "ecv: ECV_MULTITHERMAL" not in content:
            content = insert_ecv_opes(content, temperature)
        else:
            content = update_temp_max(content, temperature)
    
    return content

def setup_replica_folder(setup: ReplicaSetup, replica_num: int, 
                        temperature: Optional[int] = None,
                        source_folder: str = None) -> None:
    """Set up complete replica folder with all necessary files."""
    replica_path = os.path.join(setup.base_path, str(replica_num))
    os.makedirs(replica_path, exist_ok=True)
    
    # Generate modified PLUMED content for this replica
    content = setup_replica(setup, replica_num, setup.original_content, temperature)
    
    # Write PLUMED file
    with open(os.path.join(replica_path, 'plumed.dat'), 'w') as f:
        f.write(content)
    
    # Copy additional required files from the original subfolder
    if source_folder:  # Use the source_folder for copying files
        files_to_copy = ['npt.gro', 'topol.top', 'topol.tpr']
        if setup.reference_pdb:
            files_to_copy.append(setup.reference_pdb)
        
        copy_additional_files(source_folder, replica_path, files_to_copy)


def process_replica_setup(folder_path: str, sigma_values: Dict[str, float], 
                         plumed_base_path: str) -> None:
    """Set up the complete replica system."""
    folder_name = os.path.basename(folder_path)
    base_path = os.path.join(plumed_base_path, folder_name)
    os.makedirs(base_path, exist_ok=True)
    
    # Read original PLUMED content
    original_plumed_path = os.path.join(folder_path, "plumed.dat")
    with open(original_plumed_path, 'r') as f:
        original_content = f.read()
    
    # Get reference PDB name
    reference_pdb = get_reference_pdb(original_content)
    
    # Create setup object
    setup = ReplicaSetup(
        original_content=original_content,
        sigma_values=sigma_values,
        reference_pdb=reference_pdb,
        base_path=base_path,
        folder_name=folder_name
    )
    
    # Temperature settings for replicas 4-7
    temperatures = {
        4: 310,
        5: 330,
        6: 350,
        7: 370
    }
    
    # Set up each replica
    for i in range(8):
        setup_replica_folder(setup, i, temperatures.get(i), folder_path)  # Pass folder_path here
    # Generate SLURM batch script
    generate_sbatch_script(base_path, folder_name)
    
    print(f"Replica setup completed for {folder_name}")

def main():
    parser = argparse.ArgumentParser(description="Complete replica setup pipeline")
    parser.add_argument("--nsteps", type=int, default=500000,
                       help="Number of MD steps for production (default: 500000)")
    parser.add_argument("--ncores", type=int, default=1,
                       help="Number of CPU cores to use (default: 1)")
    parser.add_argument("--skip_production", action="store_true",
                       help="Skip the production run step")
    args = parser.parse_args()
    
    try:
        # Step 1: Select input folder
        input_folder = select_input_folder()
        prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
        base_dir = os.path.join(prep_dir, input_folder)
        
        print(f"\nProcessing system_preparation folder: {input_folder}")
        
        # Step 2: Process all subfolders
        subfolders = [d for d in os.listdir(base_dir) 
                     if os.path.isdir(os.path.join(base_dir, d)) and not d == "plumed"]
        
        if not subfolders:
            raise RuntimeError("No subfolders found")
        
        # Create plumed base directory
        plumed_base_path = os.path.join(base_dir, "plumed")
        os.makedirs(plumed_base_path, exist_ok=True)
        
        successful = 0
        failed = 0
        
        for subfolder in sorted(subfolders):
            subfolder_path = os.path.join(base_dir, subfolder)
            print(f"\nProcessing subfolder: {subfolder}")
            
            try:
                # Step 3: Run production simulation (if not skipped)
                if not args.skip_production:
                    if not run_production(subfolder_path, args.nsteps, args.ncores):
                        print(f"Skipping further processing for {subfolder} due to production failure")
                        failed += 1
                        continue
                
                # Step 4: Calculate standard deviations
                colvar_file = os.path.join(subfolder_path, "COLVAR")
                if not os.path.exists(colvar_file):
                    print(f"COLVAR file not found in {subfolder}")
                    failed += 1
                    continue
                
                std_devs = calculate_std_dev(colvar_file)
                std_output = os.path.join(subfolder_path, "colvar_results.txt")
                write_std_results(std_devs, std_output)
                
                # Step 5: Set up replica system
                process_replica_setup(subfolder_path, std_devs, plumed_base_path)
                
                successful += 1
                
            except Exception as e:
                print(f"Error processing {subfolder}: {str(e)}")
                import traceback
                traceback.print_exc()
                failed += 1
                continue
        
        print("\nProcessing Summary:")
        print(f"Total subfolders processed: {len(subfolders)}")
        print(f"Successful: {successful}")
        print(f"Failed: {failed}")
        
        if successful > 0:
            print(f"\nReplica systems have been set up in: {plumed_base_path}")
        
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()