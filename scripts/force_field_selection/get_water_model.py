import os
import sys

def get_water_model_content(water_model):
    water_model_path = os.path.abspath(os.path.join(
        os.getcwd(), '..', '..', 'system_parameters', 'water_force_fields',
        water_model, f'{water_model}_atomtypes.txt'
    ))
    
    if not os.path.exists(water_model_path):
        print(f"Error: Water model file not found at {water_model_path}")
        sys.exit(1)
    
    with open(water_model_path, 'r') as f:
        return f.read().strip()

def process_topol_file(file_path, water_model_content, water_model):
    with open(file_path, 'r') as f:
        content = f.readlines()

    atomtypes_section = False
    insert_index_atomtypes = -1
    insert_index_system = -1

    for i, line in enumerate(content):
        if '[ atomtypes ]' in line:
            atomtypes_section = True
        elif atomtypes_section and line.strip() == '':
            insert_index_atomtypes = i
            atomtypes_section = False
        elif '[ system ]' in line:
            insert_index_system = i
            break

    if insert_index_atomtypes == -1 or insert_index_system == -1:
        print(f"Error: Could not find appropriate insertion points in {file_path}")
        return

    # Insert water model atomtypes
    content.insert(insert_index_atomtypes, f"{water_model_content}\n")

    # Insert water model include statement
    content.insert(insert_index_system, f'\n#include "{water_model}.itp"\n\n')

    output_file = os.path.join(os.path.dirname(file_path), 'topol_water.top')
    with open(output_file, 'w') as f:
        f.writelines(content)

    print(f"Created {output_file}")

def select_water_model():
    water_models_dir = os.path.abspath(os.path.join(
        os.getcwd(), '..', '..', 'system_parameters', 'water_force_fields'
    ))
    
    water_models = [d for d in os.listdir(water_models_dir) 
                    if os.path.isdir(os.path.join(water_models_dir, d))]
    
    if not water_models:
        print("Error: No water models found.")
        sys.exit(1)
    
    print("Available water models:")
    for i, model in enumerate(water_models, 1):
        print(f"{i}. {model}")
    
    while True:
        try:
            choice = int(input("Enter the number of the water model you want to use: "))
            if 1 <= choice <= len(water_models):
                return water_models[choice - 1]
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

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
            choice = int(input("Enter the number of the input folder you want to process: "))
            if 1 <= choice <= len(folders):
                return folders[choice - 1]
            else:
                print("Invalid choice. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def main():
    water_model = select_water_model()
    input_folder = select_input_folder()

    water_model_content = get_water_model_content(water_model)

    prep_dir = os.path.abspath(os.path.join(os.getcwd(), '..', '..', 'system_preparation', input_folder))
    
    if not os.path.exists(prep_dir):
        print(f"Error: Input folder not found at {prep_dir}")
        sys.exit(1)

    for root, dirs, files in os.walk(prep_dir):
        for file in files:
            if file == 'topol_0.top':
                file_path = os.path.join(root, file)
                process_topol_file(file_path, water_model_content, water_model)

    print(f"\nProcessing complete. Water model {water_model} has been added to all topol_0.top files in {input_folder}.")

if __name__ == "__main__":
    main()
