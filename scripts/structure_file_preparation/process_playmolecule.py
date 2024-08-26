import os
import subprocess
import glob
import shutil
import parmed as pmd

# Function to perform operations inside each folder
def process_folder(folder_name):
    # Skip if not a folder
    if not os.path.isdir(folder_name):
        return
    
    # Navigate to the GAFF2 folder
    gaff2_folder_path = os.path.join(folder_name, 'parameters', 'GAFF2')
    os.chdir(gaff2_folder_path)
    
    # Execute tleap command
    subprocess.run(['tleap', '-f', 'tleap.in'])
    
    # Execute Python code for parmed
    parm = pmd.load_file('MOL.prmtop')
    parm.save('MOL.top', format='gromacs')
    
    # Navigate back to the main directory
    os.chdir('../../../')

# Main function to loop through all folders in the current directory
if __name__ == "__main__":
    # Create a directory one level up
    results_folder = '../playmolecule_results'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    
    for folder in os.listdir('.'):
        try:
            process_folder(folder)
        
            # Copy MOL.top to the results folder, renaming it to {name}.top
            if os.path.isdir(folder):
                mol_top_path = os.path.join(folder, 'parameters', 'GAFF2', 'MOL.top')
                if os.path.exists(mol_top_path):
                    new_mol_top_path = os.path.join(results_folder, f"{folder}.top")
                    shutil.copy(mol_top_path, new_mol_top_path)
        except:
            print("Error processing folder: {}".format(folder))
