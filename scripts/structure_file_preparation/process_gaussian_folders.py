import os
import subprocess
import glob

# Function to perform operations inside each {name}_gaussian folder
def process_folder(folder_name):
    # Extract the 'name' from the folder_name
    name = folder_name.replace('_gaussian', '')
    
    # Navigate into the folder
    os.chdir(folder_name)
    
    # Execute shell commands
    subprocess.run(['antechamber', '-i', f'{name}.esp', '-fi', 'gesp', '-o', f'{name}.prep', '-fo', 'prepc', '-c', 'resp', '-ge', f'{name}.esp', '-at', 'gaff2', '-pl', '15'])
    subprocess.run(['parmchk2', '-i', f'{name}.prep', '-f', 'prepc', '-o', f'frcmod.{name}'])
    subprocess.run(['antechamber', '-i', f'{name}.prep', '-fi', 'prepc', '-o', f'{name}.pdb', '-fo', 'pdb'])
    
    # Navigate back to the parent directory
    os.chdir('..')

# Main function to loop through all {name}_gaussian folders in the current directory
if __name__ == "__main__":
    for folder in glob.glob('*_gaussian'):
        process_folder(folder)
