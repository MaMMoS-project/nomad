import os
import argparse

# TODO make checks not optional


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
    subfolders = ['GS', 'Jij', 'MC']

    for subfolder in subfolders:
        subfolder_path = os.path.join(folder_path, subfolder)
        if os.path.isdir(subfolder_path):
            if verbose:
                print(
                    f"Subfolder '{subfolder}' exists in the folder '{folder_path}'.")
            all_exist = True
        else:
            if verbose:
                print(
                    f"Subfolder '{subfolder}' does not exist in the folder '{folder_path}'.")
            all_exist = False
    return all_exist


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


def check_structure(folder_path, check_README=False, check_subfolders=True, check_structure_cif=True, verbose=True):
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
        if check_structure_cif:
            structure_cif_exists = check_structure_cif_file_exists(
                folder_path, verbose)

        if readme_exists and subfolders_exist and structure_cif_exists:
            if verbose:
                print("# All checks passed.")
            return True
        else:
            if verbose:
                print("# Some checks failed.")
            return False

# TODO: Add a check if the main folder exists
# Optional-TODO: Add a check if dataset is .zip file and extract it


def main(args):

    dataset_name = 'Co2Fe2H4'  # Default value

    if args.dataset_name:
        # print("The dataset is given: " +
        #       args.dataset_name)
        dataset_name = args.dataset_name

    check_structure(args.dataset_name)


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
