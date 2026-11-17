import os
import argparse


def check_readme_exists(folder_path, verbose='on'):
    readme_path = os.path.join(folder_path, 'README')
    if os.path.isfile(readme_path):
        if verbose in ['on', 'debug']:
            print(f"File 'README' exists.")
        return True
    else:
        if verbose in ['on', 'debug']:
            print(
                f"File 'README' does not exist in the folder '{folder_path}'.")
        return False


def check_subfolders_exist(folder_path, check_optional_subf=True, verbose='on'):
    mandatory_subfolders = ['GS', 'GS/x', 'GS/z', 'Jij', 'MC']
    non_mandatory_subfolders = ['GS/y']

    mandatory_exist = []
    non_mandatory_exist = []

    if verbose in ['debug']:
        print("# In check_subfolders_exist function.")

    for subfolder in mandatory_subfolders:
        subfolder_path = os.path.join(folder_path, subfolder)
        if os.path.isdir(subfolder_path):
            if verbose in ['on', 'debug']:
                print(f"Mandatory subfolder '{subfolder}' exists.")
            mandatory_exist.append(True)
        else:
            if verbose in ['on', 'debug']:
                print(
                    f"Mandatory subfolder '{subfolder}' does not exist in the folder '{folder_path}'.")
            mandatory_exist.append(False)
    if check_optional_subf:
        for subfolder in non_mandatory_subfolders:
            subfolder_path = os.path.join(folder_path, subfolder)
            if os.path.isdir(subfolder_path):
                if verbose in ['on', 'debug']:
                    print(f"Optional subfolder '{subfolder}' exists.")
                non_mandatory_exist.append(True)
            else:
                if verbose in ['on', 'debug']:
                    print(
                        f"Optional subfolder '{subfolder}' does not exist in the folder '{folder_path}'.")
                non_mandatory_exist.append(False)
    else:
        non_mandatory_exist = [None] * len(non_mandatory_subfolders)
        if verbose in ['on', 'debug']:
            print("Optional subfolders are not checked.")
        if verbose in ['debug']:
            print("non_mandatory_exist: ", non_mandatory_exist)

    return all(mandatory_exist), non_mandatory_exist


def check_structure_cif_file_exists(folder_path, verbose='on'):
    structure_file_path = os.path.join(folder_path, 'structure.cif')
    if os.path.isfile(structure_file_path):
        if verbose in ['on', 'debug']:
            print(f"File 'structure.cif' exists.")
        return True
    else:
        if verbose in ['on', 'debug']:
            print(f"File 'structure.cif' does not exist.")
        return False


def check_out_last_file_exists(folder_path, verbose='on'):
    out_last_file_path = os.path.join(folder_path, 'out_last')
    if os.path.isfile(out_last_file_path):
        if verbose in ['on', 'debug']:
            print(f"File 'out_last' exists in the folder '{folder_path}'.")
        return True
    else:
        if verbose in ['on', 'debug']:
            print(
                f"File 'out_last' does not exist in the folder '{folder_path}'.")
        return False


def check_posfile_file_exists(folder_path, verbose='on'):
    posfile_file_path = os.path.join(folder_path, 'posfile')
    if os.path.isfile(posfile_file_path):
        if verbose in ['on', 'debug']:
            print(
                f"File 'posfile' exists in the folder '{folder_path}'.")
        return True
    else:
        if verbose in ['on', 'debug']:
            print(
                f"File 'posfile' does not exist in the folder '{folder_path}'.")
        return False

# TODO: check is jfile needed anyhow
def check_jfile_file_exists(folder_path, verbose='on'):
    jfile_file_path = os.path.join(folder_path, 'jfile')
    if os.path.isfile(jfile_file_path):
        if verbose in ['on', 'debug']:
            print(
                f"File 'jfile' exists in the folder '{folder_path}'.")
        return True
    else:
        if verbose in ['on', 'debug']:
            print(
                f"File 'jfile' does not exist in the folder '{folder_path}'.")
        return False


def check_momfile_file_exists(folder_path, verbose='on'):
    momfile_file_path = os.path.join(folder_path, 'momfile')
    if os.path.isfile(momfile_file_path):
        if verbose in ['on', 'debug']:
            print(
                f"File 'momfile' exists in the folder '{folder_path}'.")
        return True
    else:
        if verbose in ['on', 'debug']:
            print(
                f"File 'momfile' does not exist in the folder '{folder_path}'.")
        return False


def check_structure(
    folder_path,
    check_subfolders=True,
    check_optional_subf=True,
    check_README=False,
    check_structure_cif=True,
    check_out_last_files=True,
    verbose='on',
):

    if verbose in ['debug']:
        print("# In check_structure function.")
    main_folder_exists = False
    readme_exists = False
    mandatory_subfolders_exist = False
    optional_subfolders_exist = False
    structure_cif_exists = False
    out_last_x_exists = False
    out_last_y_exists = False
    out_last_z_exists = False

    if not os.path.isdir(folder_path):
        if verbose in ['on', 'debug']:
            print(
                f"The main folder '{folder_path}' for the dataset does not exist.")
        main_folder_exists = False
    else:
        if verbose in ['on', 'debug']:
            print(f"The main folder '{folder_path}' for the dataset exists.")
        main_folder_exists = True

        if check_README:
            readme_exists = check_readme_exists(folder_path, verbose)
        if check_subfolders:
            mandatory_subfolders_exist, optional_subfolders_exist = check_subfolders_exist(
                folder_path, check_optional_subf, verbose)
            if verbose in ['debug']:
                print(
                    f"Mandatory subfolders exist: {mandatory_subfolders_exist}")
                print(
                    f"Optional subfolders exist: {optional_subfolders_exist}")
            if mandatory_subfolders_exist:
                if check_out_last_files:
                    out_last_x_exists = check_out_last_file_exists(
                        folder_path+'/GS/x', verbose)
                    out_last_z_exists = check_out_last_file_exists(
                        folder_path+'/GS/z', verbose)

                posfile_exists = check_posfile_file_exists(
                    folder_path+'/MC', verbose)
                jfile_exists = check_jfile_file_exists(
                    folder_path+'/MC', verbose)
                momfile_exists = check_momfile_file_exists(
                    folder_path+'/MC', verbose)

                if (posfile_exists or momfile_exists):
                    if verbose in ['on', 'debug']:
                        print(
                            f"Files 'posfile', 'jfile', or 'momfile' exist in the folder '{folder_path}/MC'.")
                    mom_pos_file_exists = True
                else:
                    if verbose in ['on', 'debug']:
                        print(
                            f"None of the files 'posfile', 'jfile', or 'momfile' exist in the folder '{folder_path}/MC'.")
                    mom_pos_file_exists = False

            if optional_subfolders_exist and check_optional_subf:
                if check_out_last_files:
                    out_last_y_exists = check_out_last_file_exists(
                        folder_path+'/GS/y', verbose)
        if check_structure_cif:
            structure_cif_exists = check_structure_cif_file_exists(
                folder_path, verbose)

        return_values = [main_folder_exists, readme_exists, mandatory_subfolders_exist,
                         structure_cif_exists, out_last_x_exists, out_last_z_exists,
                         mom_pos_file_exists]

        if verbose == 'debug':
            print(f"All return values: {return_values}")

        if all(return_values):
            if optional_subfolders_exist == [True] * len(optional_subfolders_exist):
                if verbose in ['on', 'debug']:
                    print("## Optional subfolders exist.")
                    if out_last_y_exists:
                        print("## Optional out_last_y file exists.")
            else:
                if optional_subfolders_exist == [False] * len(optional_subfolders_exist):
                    if verbose in ['on', 'debug']:
                        print("## Optional subfolders do not exist.")
                else:
                    if optional_subfolders_exist == [None] * len(optional_subfolders_exist):
                        if verbose in ['on', 'debug']:
                            print("## Optional subfolders are not checked.")
            if verbose in ['off', 'on', 'debug']:
                print("# All checks passed.")
            return True
        else:
            if verbose in ['off', 'on', 'debug']:
                print("# Some checks failed.")
            return False


# --- New uppsala-data structure (https://github.com/MaMMoS-project/uppsala-data) ---
# Dataset root: metadata.yaml, intrinsic_properties.yaml, structure.cif (required),
#               DOSCAR and README.md (optional)
# RSPt/:        common_input/{atomdens,kmap,spts,symcof,symt.inp} (all required),
#               gs_x/ and gs_z/ (required), gs_y/ (optional); each gs_* contains
#               data + out_last + (hist XOR out_MF); Jij/ with data + out_last
# UppASD/:      MC_1/ (required) with jfile, momfile, posfile, inpsd.dat,
#               thermal.csv, thermal.dat; MC_2/ and MC_3/ optional (same content)

UPPSALA_ROOT_FILES = ['metadata.yaml', 'intrinsic_properties.yaml', 'structure.cif']
UPPSALA_COMMON_INPUT_FILES = ['atomdens', 'kmap', 'spts', 'symcof', 'symt.inp']
UPPSALA_GS_DIRS_REQUIRED = ['gs_x', 'gs_z']
UPPSALA_GS_DIRS_OPTIONAL = ['gs_y']
UPPSALA_GS_FILES = ['data', 'out_last']
UPPSALA_MC_DIRS_REQUIRED = ['MC_1']
UPPSALA_MC_DIRS_OPTIONAL = ['MC_2', 'MC_3']
UPPSALA_MC_FILES = ['jfile', 'momfile', 'posfile', 'inpsd.dat', 'thermal.csv', 'thermal.dat']


def is_uppsala_structure(folder_path):
    """Detect the new uppsala-data layout by the presence of RSPt/ and UppASD/."""
    return os.path.isdir(os.path.join(folder_path, 'RSPt')) and os.path.isdir(
        os.path.join(folder_path, 'UppASD'))


def check_uppsala_structure(folder_path, verbose='on'):
    """Check a dataset against the new uppsala-data layout (see dataset-schema.yaml
    in https://github.com/MaMMoS-project/uppsala-data). Returns True/False."""
    results = []

    def record(ok, msg):
        if verbose in ['on', 'debug']:
            print(('OK   ' if ok else 'FAIL ') + msg)
        results.append(bool(ok))

    if not os.path.isdir(folder_path):
        print(f"The main folder '{folder_path}' for the dataset does not exist.")
        return False

    for f in UPPSALA_ROOT_FILES:
        record(os.path.isfile(os.path.join(folder_path, f)), f"file '{f}'")

    rspt = os.path.join(folder_path, 'RSPt')
    common_input = os.path.join(rspt, 'common_input')
    record(os.path.isdir(common_input), "directory 'RSPt/common_input'")
    if os.path.isdir(common_input):
        for f in UPPSALA_COMMON_INPUT_FILES:
            record(os.path.isfile(os.path.join(common_input, f)),
                   f"file 'RSPt/common_input/{f}'")

    for gs in UPPSALA_GS_DIRS_REQUIRED:
        gs_dir = os.path.join(rspt, gs)
        record(os.path.isdir(gs_dir), f"directory 'RSPt/{gs}'")
        if os.path.isdir(gs_dir):
            for f in UPPSALA_GS_FILES:
                record(os.path.isfile(os.path.join(gs_dir, f)),
                       f"file 'RSPt/{gs}/{f}'")
            has_hist = os.path.isfile(os.path.join(gs_dir, 'hist'))
            has_out_mf = os.path.isfile(os.path.join(gs_dir, 'out_MF'))
            record(has_hist != has_out_mf,
                   f"'RSPt/{gs}' has exactly one of 'hist'/'out_MF'")

    record(os.path.isdir(os.path.join(rspt, 'Jij')), "directory 'RSPt/Jij'")
    jij = os.path.join(rspt, 'Jij')
    if os.path.isdir(jij):
        for f in ['data', 'out_last']:
            record(os.path.isfile(os.path.join(jij, f)), f"file 'RSPt/Jij/{f}'")

    uppasd = os.path.join(folder_path, 'UppASD')
    for mc in UPPSALA_MC_DIRS_REQUIRED + UPPSALA_MC_DIRS_OPTIONAL:
        mc_dir = os.path.join(uppasd, mc)
        if mc in UPPSALA_MC_DIRS_OPTIONAL and not os.path.isdir(mc_dir):
            continue
        record(os.path.isdir(mc_dir), f"directory 'UppASD/{mc}'")
        if os.path.isdir(mc_dir):
            for f in UPPSALA_MC_FILES:
                record(os.path.isfile(os.path.join(mc_dir, f)),
                       f"file 'UppASD/{mc}/{f}'")

    if all(results):
        if verbose in ['off', 'on', 'debug']:
            print("# All uppsala-data structure checks passed.")
        return True
    if verbose in ['off', 'on', 'debug']:
        print("# Some uppsala-data structure checks failed.")
    return False


def validate_with_mammos_parser(folder_path, verbose='on'):
    """Run `mammos-parser uppsala-data validate-dataset` if the CLI is installed.
    Looks on PATH and next to the current Python executable (pixi env bin).
    Returns True/False, or None if mammos-parser is not available."""
    import shutil
    import subprocess
    import sys

    candidates = [shutil.which('mammos-parser')]
    env_bin_exe = os.path.join(os.path.dirname(sys.executable), 'mammos-parser')
    if os.path.isfile(env_bin_exe):
        candidates.append(env_bin_exe)
    exe = next((c for c in candidates if c), None)
    if exe is None:
        if verbose in ['on', 'debug']:
            print("mammos-parser not found; skipping external validation.")
        return None
    proc = subprocess.run(
        [exe, 'uppsala-data', 'validate-dataset', folder_path],
        capture_output=True, text=True)
    if verbose in ['on', 'debug']:
        if proc.stdout:
            print(proc.stdout)
        if proc.returncode != 0 and proc.stderr:
            print(proc.stderr)
    return proc.returncode == 0


def check_dataset_structure(folder_path, verbose='on'):
    """Dispatch on the detected layout: new uppsala-data structure or old
    GS/Jij/MC structure. Returns True/False, or None if only the external
    mammos-parser check ran and it is unavailable."""
    if is_uppsala_structure(folder_path):
        internal_ok = check_uppsala_structure(folder_path, verbose=verbose)
        external = validate_with_mammos_parser(folder_path, verbose=verbose)
        if external is False:
            return False
        return internal_ok
    return check_structure(folder_path, check_README=True, verbose=verbose)


def main(args):
    dataset_name = 'Co2Fe2H4'  # Default value
    check_optional_subf = True  # Default value
    verbose = 'on'

    if args.dataset_name:
        dataset_name = args.dataset_name

    if args.check_optional_subf is not None:
        check_optional_subf = args.check_optional_subf

    if args.verbose:
        verbose = args.verbose

    check_structure(dataset_name, check_optional_subf=check_optional_subf, check_README=True,
                    check_subfolders=True, check_structure_cif=True, check_out_last_files=True, verbose=verbose)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check folder structure of a simulation data set for the 'Electronic' and 'Atomistic' model.")
    parser.add_argument('-dataset_name', type=str, default='Co2Fe2H4',
                        help='Name of the folder or dataset to check (default: Co2Fe2H4)')
    parser.add_argument('-check_optional_subf', action='store_true',
                        dest='check_optional_subf', help='Check optional subfolders (default: True)')
    parser.add_argument('-no_check_optional_subf', action='store_false',
                        dest='check_optional_subf', help='Do not check optional subfolders')
    parser.add_argument('-verbose', type=str, choices=[
                        'off', 'on', 'debug'], default='on', help='Set verbosity level (default: on)')

    args = parser.parse_args()
    main(args)
