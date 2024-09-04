import os
import subprocess
import glob
import shutil
import parmed as pmd

# Function to perform operations inside each folder
def process_folder(folder_name, base_dir):
    # Skip if not a folder
    if not os.path.isdir(os.path.join(base_dir, folder_name)):
        return
    
    # Navigate to the GAFF2 folder
    gaff2_folder_path = os.path.join(base_dir, folder_name, 'parameters', 'GAFF2')
    os.chdir(gaff2_folder_path)
    
    mol_top_path = os.path.join(gaff2_folder_path, 'MOL.top')
    
    # If MOL.top doesn't exist, create it
    if not os.path.exists(mol_top_path):
        # Execute tleap command
        subprocess.run(['tleap', '-f', 'tleap.in'])
        
        # Execute Python code for parmed
        parm = pmd.load_file('MOL.prmtop')
        parm.save('MOL.top', format='gromacs')
    
    # Navigate back to the base directory
    os.chdir(base_dir)
    
    return mol_top_path

def main():
    # Get the current working directory
    current_dir = os.getcwd()
    
    # Set the base directory to ../../playmolecule_results
    base_dir = os.path.abspath(os.path.join(current_dir, '..', '..', 'playmolecule_results'))
    
    # Ask user for the output folder name
    output_folder_name = input("Enter a name for the output folder (it will be created in ../../system_parameters/topologies/): ")
    
    # Set the output directory
    output_dir = os.path.abspath(os.path.join(current_dir, '..', '..', 'system_parameters', 'topologies', output_folder_name))
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Working in directory: {base_dir}")
    print(f"Results will be saved in: {output_dir}")
    
    # Change to the base directory
    os.chdir(base_dir)
    
    # Process all folders in the playmolecule_results directory
    for folder in os.listdir(base_dir):
        try:
            mol_top_path = process_folder(folder, base_dir)
            
            if mol_top_path and os.path.exists(mol_top_path):
                new_mol_top_path = os.path.join(output_dir, f"{folder}.top")
                shutil.copy(mol_top_path, new_mol_top_path)
                print(f"Copied {mol_top_path} to {new_mol_top_path}")
            else:
                print(f"Warning: MOL.top not found for folder {folder}")
        except Exception as e:
            print(f"Error processing folder: {folder}. Error: {str(e)}")

if __name__ == "__main__":
    main()
