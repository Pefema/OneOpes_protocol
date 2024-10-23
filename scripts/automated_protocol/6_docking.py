def prepare_receptor(pdb_file):
    """Prepare receptor using MGLTools prepare_receptor4 script with minimal processing."""
    # Save original file with _not_docked suffix
    not_docked_file = pdb_file.replace('.pdb', '_not_docked.pdb')
    if not os.path.exists(not_docked_file):
        subprocess.run(['cp', pdb_file, not_docked_file], check=True)

    # Get paths
    abs_pdb_file = os.path.abspath(pdb_file)
    output_pdbqt = abs_pdb_file.replace('.pdb', '.pdbqt')
    
    # Get directory and filename
    working_dir = os.path.dirname(abs_pdb_file)
    pdb_filename = os.path.basename(abs_pdb_file)
    output_filename = os.path.basename(output_pdbqt)
    
    # Store current directory
    original_dir = os.getcwd()
    
    try:
        # Change to working directory
        os.chdir(working_dir)
        
        # Run prepare_receptor4 with local filenames
        cmd = f"prepare_receptor4.py -r {pdb_filename} -o {output_filename}"
        print(f"Running command in {working_dir}: {cmd}")
        
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
            
        # Verify the output file exists
        if not os.path.exists(output_filename):
            raise Exception(f"PDBQT file was not created: {output_filename}")
            
        return output_pdbqt
        
    finally:
        # Always return to original directory
        os.chdir(original_dir)

def prepare_ligand(pdb_file):
    """Prepare ligand using MGLTools prepare_ligand4 script while strictly preserving structure."""
    # Save original file with _not_docked suffix
    not_docked_file = pdb_file.replace('.pdb', '_not_docked.pdb')
    if not os.path.exists(not_docked_file):
        subprocess.run(['cp', pdb_file, not_docked_file], check=True)

    # Get paths
    abs_pdb_file = os.path.abspath(pdb_file)
    output_pdbqt = abs_pdb_file.replace('.pdb', '.pdbqt')
    
    # Get directory and filename
    working_dir = os.path.dirname(abs_pdb_file)
    pdb_filename = os.path.basename(abs_pdb_file)
    output_filename = os.path.basename(output_pdbqt)
    
    # Store current directory
    original_dir = os.getcwd()
    
    try:
        # Change to working directory
        os.chdir(working_dir)
        
        # Run prepare_ligand4 with local filenames
        cmd = f"prepare_ligand4.py -l {pdb_filename} -o {output_filename} -U -C -p -k -Z"
        print(f"Running command in {working_dir}: {cmd}")
        
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
            
        # Verify the output file exists
        if not os.path.exists(output_filename):
            raise Exception(f"PDBQT file was not created: {output_filename}")
            
        return output_pdbqt
        
    finally:
        # Always return to original directory
        os.chdir(original_dir)

def dock_molecules_in_folder(input_folder):
    """Process all subfolders using the selected host against each guest."""
    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation'))
    base_dir = os.path.join(prep_dir, input_folder)

    # Get all subfolders
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not subfolders:
        print("Error: No subfolders found")
        return

    # Select host PDB from the first subfolder
    first_subfolder = os.path.join(base_dir, subfolders[0])
    host_basename = select_host_pdb(first_subfolder)
    if not host_basename:
        print("Error: No host PDB file selected")
        return

    host_file = os.path.join(first_subfolder, host_basename)
    print(f"\nSelected host file: {host_basename}")

    # Prepare the receptor once since it will be used for all docking runs
    try:
        print(f"\nPreparing receptor: {host_file}")
        receptor_pdbqt = prepare_receptor(host_file)
        print(f"Receptor preparation successful: {receptor_pdbqt}")
    except Exception as e:
        print(f"Error preparing receptor: {str(e)}")
        return

    # Process each subfolder
    for subfolder in subfolders:
        subfolder_path = os.path.join(base_dir, subfolder)
        print(f"\nProcessing subfolder: {subfolder}")

        # Find guest PDB file
        guest_basename = find_guest_pdb(subfolder_path, host_basename)
        if not guest_basename:
            print(f"Skipping {subfolder}: No guest PDB file found")
            continue

        guest_file = os.path.join(subfolder_path, guest_basename)
        print(f"Found guest file: {guest_basename}")

        try:
            # Copy host files to current subfolder if not already there
            if not os.path.exists(os.path.join(subfolder_path, host_basename)):
                subprocess.run(['cp', host_file, subfolder_path], check=True)
            if not os.path.exists(os.path.join(subfolder_path, os.path.basename(receptor_pdbqt))):
                subprocess.run(['cp', receptor_pdbqt, subfolder_path], check=True)

            print(f"Preparing ligand: {guest_file}")
            ligand_pdbqt = prepare_ligand(guest_file)

            print("Calculating search box...")
            center, size = calculate_box(os.path.join(subfolder_path, os.path.basename(receptor_pdbqt)))

            print("Creating configuration file...")
            config_path = create_config_file(center, size, 
                                          os.path.join(subfolder_path, os.path.basename(receptor_pdbqt)), 
                                          os.path.join(subfolder_path, os.path.basename(ligand_pdbqt)), 
                                          subfolder_path)

            print("Running AutoDock Vina...")
            docked_pdbqt = run_vina(config_path, subfolder_path)

            # Save docked complex with original filename
            output_pdb = os.path.join(subfolder_path, guest_basename)

            print("Converting docked complex to PDB...")
            convert_pdbqt_to_pdb_mda(docked_pdbqt, output_pdb)

            print(f"Docking completed successfully for {subfolder}")
            print(f"Complex saved as: {guest_basename}")

        except Exception as e:
            print(f"Error processing {subfolder}: {str(e)}")
            continue
