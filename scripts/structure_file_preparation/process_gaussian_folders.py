import os
import subprocess
import glob
import shutil

def process_folder(folder_name, base_dir):
    # Extract the 'name' from the folder_name
    name = folder_name.replace('_gaussian', '')
    
    # Navigate into the folder
    os.chdir(os.path.join(base_dir, folder_name))
    
    # Execute shell commands
    subprocess.run(['antechamber', '-i', f'{name}.esp', '-fi', 'gesp', '-o', f'{name}.prep', '-fo', 'prepc', '-c', 'resp', '-ge', f'{name}.esp', '-at', 'gaff2', '-pl', '15'])
    subprocess.run(['parmchk2', '-i', f'{name}.prep', '-f', 'prepc', '-o', f'frcmod.{name}'])
    subprocess.run(['antechamber', '-i', f'{name}.prep', '-fi', 'prepc', '-o', f'{name}.pdb', '-fo', 'pdb'])
    
    # Navigate back to the base directory
    os.chdir(base_dir)
    
    return name

def copy_resp_files(output_folder, processed_folders, base_dir):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    for folder in processed_folders:
        resp_file = os.path.join(base_dir, folder, 'ANTECHAMBER_RESP2.OUT')
        if os.path.exists(resp_file):
            new_name = f"{folder.replace('_gaussian', '')}_RESP2.OUT"
            shutil.copy(resp_file, os.path.join(output_folder, new_name))
            print(f"Copied {resp_file} to {os.path.join(output_folder, new_name)}")
        else:
            print(f"Warning: ANTECHAMBER_RESP2.OUT not found in {folder}")

def main():
    # Get the current working directory
    current_dir = os.getcwd()
    
    # Set the base directory to ../../gaussian_results
    base_dir = os.path.abspath(os.path.join(current_dir, '..', '..', 'gaussian_results'))
    
    # Set the charges directory to ../../system_parameters/charges
    charges_dir = os.path.abspath(os.path.join(current_dir, '..', '..', 'system_parameters', 'charges'))
    
    # Change to the base directory
    os.chdir(base_dir)
    print(f"Working in directory: {base_dir}")
    
    processed_folders = []
    
    # Process all *_gaussian folders in the gaussian_results directory
    for folder in glob.glob('*_gaussian'):
        print(f"Processing {folder}")
        name = process_folder(folder, base_dir)
        processed_folders.append(folder)
    
    # Ask user for the output folder name
    output_folder_name = input("Enter a name for the output folder (it will be created in ../../system_parameters/charges/): ")
    output_folder = os.path.join(charges_dir, output_folder_name)
    
    # Copy RESP files to the new folder
    copy_resp_files(output_folder, processed_folders, base_dir)

if __name__ == "__main__":
    main()
