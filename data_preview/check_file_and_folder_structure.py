import os
import argparse


def check_readme_exists(folder_path, verbose=True):
    readme_path = os.path.join(folder_path, 'README')
    if os.path.isfile(readme_path):
        if verbose:
            print(f"File 'README' exists.")
        return True
    else:
        print(f"File 'README' does not exist in the folder '{folder_path}'.")
        return False


def check_subfolders_exist(folder_path, check_optional_subf=True, verbose=True):
    mandatory_subfolders = ['GS', 'GS/x', 'GS/z', 'Jij', 'MC']
    non_mandatory_subfolders = ['GS/y']
    
    mandatory_exist = []
    non_mandatory_exist = []

    for subfolder in mandatory_subfolders:
        subfolder_path = os.path.join(folder_path, subfolder)
        if os.path.isdir(subfolder_path):
            if verbose:
                print(
                    f"Mandatory subfolder '{subfolder}' exists.")
            mandatory_exist.append(True)
        else:
            if verbose:
                print(
                    f"Mandatory subfolder '{subfolder}' does not exist in the folder '{folder_path}'.")
            mandatory_exist.append(False)
    if check_optional_subf:
        for subfolder in non_mandatory_subfolders:
            subfolder_path = os.path.join(folder_path, subfolder)
            if os.path.isdir(subfolder_path):
                if verbose:
                    print(
                        f"Non-mandatory subfolder '{subfolder}'.")
                non_mandatory_exist.append(True)
            else:
                if verbose:
                    print(
                        f"Non-mandatory subfolder '{subfolder}' does not exist in the folder '{folder_path}'.")
                non_mandatory_exist.append(False)

    return all(mandatory_exist), all(non_mandatory_exist)


def check_structure_cif_file_exists(folder_path, verbose=True):
    structure_file_path = os.path.join(folder_path, 'structure.cif')
    if os.path.isfile(structure_file_path):
        if verbose:
            print(
                f"File 'structure.cif' exists.")
        return True
    else:
        if verbose:
            print(
                f"File 'structure.cif' does not exist.")
        return False


def check_out_last_file_exists(folder_path, verbose=True):
    out_last_file_path = os.path.join(folder_path, 'out_last')
    if os.path.isfile(out_last_file_path):
        if verbose:
            print(
                f"File 'out_last' exists in the folder '{folder_path}'.")
        return True
    else:
        if verbose:
            print(
                f"File 'out_last' does not exist in the folder '{folder_path}'.")
        return False


def check_structure(folder_path, check_subfolders=True, check_optional_subf=True, check_README=False, check_structure_cif=True, check_out_last_files=True, verbose=True):
    """
    folder_path (str): The path to the folder to check.
    check_README (bool): Whether to check for the existence of a README file. Default is False.
    check_subfolders (bool): Whether to check for the existence of mandatory and optional subfolders. Default is True.
    check_optional_subf (bool): Whether to check for the existence of optional subfolders. Default is True.
    check_structure_cif (bool): Whether to check for the existence of a structure.cif file. Default is True.
    check_out_last_files (bool): Whether to check for the existence of out_last_x, out_last_y, and out_last_z files. Default is True.
    verbose (bool): Whether to print detailed messages during the check. Default is True.
    Returns:
    bool: True if all specified checks pass, False otherwise.
    """

    main_folder_exists = False
    readme_exists = False
    mandatory_subfolders_exist = False
    optional_subfolders_exist = False
    structure_cif_exists = False
    out_last_x_exists = False
    out_last_y_exists = False
    out_last_z_exists = False
    
    if not os.path.isdir(folder_path):
        if verbose:
            print(
                f"The main folder '{folder_path}' for the dataset does not exist.")
        main_folder_exists = False
    else:
        if verbose:
            print(f"The main folder '{folder_path}' for the dataset exists.")
        main_folder_exists = True

        if check_README:
            readme_exists = check_readme_exists(folder_path, verbose)
        if check_subfolders:
            mandatory_subfolders_exist, optional_subfolders_exist = check_subfolders_exist(folder_path, check_optional_subf, verbose)
            if mandatory_subfolders_exist:
                if check_out_last_files:
                    out_last_x_exists = check_out_last_file_exists(
                        folder_path+'/GS/x', verbose)
                    out_last_z_exists = check_out_last_file_exists(
                        folder_path+'/GS/z', verbose)
            if optional_subfolders_exist and check_optional_subf:
                if check_out_last_files:
                    out_last_y_exists = check_out_last_file_exists(
                        folder_path+'/GS/y', verbose)
        if check_structure_cif:
            structure_cif_exists = check_structure_cif_file_exists(
                folder_path, verbose)

        return_values = [main_folder_exists, readme_exists,
                         mandatory_subfolders_exist, structure_cif_exists, out_last_x_exists, out_last_z_exists]

        if verbose:
            print(f"All return values: {return_values}")
        
        if all(return_values):
            if optional_subfolders_exist == True:
                if verbose:
                    print("## Optional subfolders exist.")
                    if out_last_y_exists:
                        print("## Optional out_last_y file exists.")
            else:
                if verbose:
                    print("## Optional subfolders do not exist.") 
            if verbose:
                print("# All checks passed.")
            return True
        else:
            if verbose:
                print("# Some checks failed.")
            return False


def main(args):

    dataset_name = 'Co2Fe2H4'  # Default value

    check_optional_subf = True  # Default value
    
    verbose = True
    
    if args.dataset_name:
        if verbose:
            dataset_name = args.dataset_name

    if args.check_optional_subf is None:
        print(" args.check_optional_subf is None")
    
    if args.check_optional_subf is not None:        
        if args.check_optional_subf == True:
            if verbose:
                print("Check optional subfolders : True")
        else:
            if verbose:
                print("Check optional subfolders : False")
        check_optional_subf = args.check_optional_subf


    check_structure(dataset_name, check_optional_subf=check_optional_subf, check_README=True, check_subfolders=True, check_structure_cif=True, check_out_last_files=True, verbose=verbose)


if __name__ == "__main__":
    # Initialize parser
    parser = argparse.ArgumentParser(
        description="Check folder structure of a simulation data set for the 'Electronic' and 'Atomistic' model.")

    # Adding optional argument
    parser.add_argument('-dataset_name', type=str, default='Co2Fe2H4',
                        help='Name of the folder or dataset to check (default: Co2Fe2H4)')

    parser.add_argument('-check_optional_subf', action='store_true', dest='check_optional_subf',
                        help='Check optional subfolders (default: True)')
    parser.add_argument('-no_check_optional_subf', action='store_false', dest='check_optional_subf',
                        help='Do not check optional subfolders')

    # Read arguments from command line
    args = parser.parse_args()

    print(args)

    main(args)
