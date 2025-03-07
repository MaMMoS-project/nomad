import os
import argparse


def check_readme_exists(folder_path, verbose=True):
    readme_path = os.path.join(folder_path, 'README')
    if os.path.isfile(readme_path):
        if verbose:
            print(f"File 'README' exists in the folder '{folder_path}'.")
        return True
    else:
        print(f"File 'README' does not exist in the folder '{folder_path}'.")
        return False


def check_subfolders_exist(folder_path, verbose=True):
    subfolders = ['GS', 'GS/x', 'GS/y', 'GS/z', 'Jij', 'MC']
    all_exist = []

    for subfolder in subfolders:
        subfolder_path = os.path.join(folder_path, subfolder)
        if os.path.isdir(subfolder_path):
            if verbose:
                print(
                    f"Subfolder '{subfolder}' exists in the folder '{folder_path}'.")
            all_exist.append(True)
        else:
            if verbose:
                print(
                    f"Subfolder '{subfolder}' does not exist in the folder '{folder_path}'.")
            all_exist.append(False)
    return all(all_exist)


def check_structure_cif_file_exists(folder_path, verbose=True):
    structure_file_path = os.path.join(folder_path, 'structure.cif')
    if os.path.isfile(structure_file_path):
        if verbose:
            print(
                f"File 'structure.cif' exists in the folder '{folder_path}'.")
        return True
    else:
        if verbose:
            print(
                f"File 'structure.cif' does not exist in the folder '{folder_path}'.")
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


def check_structure(folder_path, check_README=False, check_subfolders=True, check_structure_cif=True, check_out_last_files=True, verbose=True):
    
    main_folder_exists = False
    readme_exists = False
    subfolders_exist = False
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
            subfolders_exist = check_subfolders_exist(folder_path, verbose)
            if subfolders_exist:
                if check_out_last_files:
                    out_last_x_exists = check_out_last_file_exists(
                        folder_path+'/GS/x', verbose)
                    out_last_y_exists = check_out_last_file_exists(
                        folder_path+'/GS/y', verbose)
                    out_last_z_exists = check_out_last_file_exists(
                        folder_path+'/GS/z', verbose)

        if check_structure_cif:
            structure_cif_exists = check_structure_cif_file_exists(
                folder_path, verbose)

        return_values = [main_folder_exists, readme_exists,
                         subfolders_exist, structure_cif_exists, out_last_x_exists, out_last_y_exists, out_last_z_exists]

        if all(return_values):
            if verbose:
                print("# All checks passed.")
            return True
        else:
            if verbose:
                print("# Some checks failed.")
            return False


def main(args):

    dataset_name = 'Co2Fe2H4'  # Default value

    if args.dataset_name:
        # print("The dataset is given: " +
        #       args.dataset_name)
        dataset_name = args.dataset_name

    check_structure(args.dataset_name, check_README=True)


if __name__ == "__main__":
    # Initialize parser
    parser = argparse.ArgumentParser(
        description="Check folder structure of a simulation data set for the 'Electronic' and 'Atomistic' model.")

    # Adding optional argument
    parser.add_argument('-dataset_name', type=str, default='Co2Fe2H4',
                        help='Name of the folder or dataset to check (default: Co2Fe2H4)')

    # Read arguments from command line
    args = parser.parse_args()

    print(args)

    main(args)
