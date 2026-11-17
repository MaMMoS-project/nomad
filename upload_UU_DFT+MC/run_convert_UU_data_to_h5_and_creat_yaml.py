import subprocess
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# DONE-TODO: allow to rund in nested folders and create a .yaml and .h5 file in the current folder
# TODO: check if the folder-name is a valid chemical formula (Co2Fe16Y6 vs. Fe3-xZnxY(x=0.22)) and make sure the formula is provided elsewhere


def resolve_datasets_dir(datasets_dir):
    """Default input is a clone of https://github.com/MaMMoS-project/uppsala-data
    (its data/ folder) next to this script; fall back to ./datasets (old layout)."""
    if datasets_dir:
        return datasets_dir
    candidates = [
        os.path.join(SCRIPT_DIR, "uppsala-data", "data"),
        "./uppsala-data/data",
        "./datasets",
    ]
    for candidate in candidates:
        if os.path.isdir(candidate):
            return candidate
    return "./datasets"


# Usage: python run_convert_UU_data_to_h5_and_creat_yaml.py [subfolder_name]
def main():
    import argparse
    import zipfile

    parser = argparse.ArgumentParser(
        description="Run convert_UU_data_to_h5.py for a subfolder and create a matching archive.yaml file."
    )
    parser.add_argument(
        "--source_dir",
        nargs="?",
        default="Fe12Co4N2",
        help="Subfolder name under the datasets directory (default: Fe12Co4N2)",
    )
    parser.add_argument(
        "--datasets_dir",
        nargs="?",
        default=None,
        help="Directory containing the datasets (default: ./uppsala-data/data if present, else ./datasets)",
    )
    args = parser.parse_args()
    subfolder = args.source_dir
    datasets_dir = resolve_datasets_dir(args.datasets_dir)
    print(f"Datasets directory: {os.path.abspath(datasets_dir)}")
    subfolder_path = os.path.join(datasets_dir, subfolder)
    if not os.path.isdir(subfolder_path):
        print(f"Error: {subfolder_path} is not a valid directory.")
        sys.exit(1)

    from packages import check_file_and_folder_structure as cffs

    structure_ok = cffs.check_dataset_structure(subfolder_path)
    if structure_ok is False:
        print(f"Error: dataset '{subfolder_path}' failed the structure validation.")
        sys.exit(1)

    script = os.path.join(SCRIPT_DIR, "convert_UU_data_to_h5.py")
    h5_file = os.path.join(datasets_dir, f"{subfolder}.h5")
    yaml_file = os.path.join(datasets_dir, f"{subfolder}.archive.yaml")

    # Create uploads directory at the same level as datasets
    uploads_dir = os.path.join(os.path.dirname(datasets_dir) or ".", "uploads")
    os.makedirs(uploads_dir, exist_ok=True)
    zip_filename = os.path.join(uploads_dir, f"{subfolder}_archive.zip")

    # Run the conversion script and create the .h5 file in datasets
    cmd = [
        sys.executable,
        script,
        "--source_dir",
        subfolder_path,
        "--output_file",
        h5_file,
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Create a new archive.yaml file from the template
    template = os.path.join(SCRIPT_DIR, "UU_template.archive.yaml")
    with open(template, "r") as f:
        content = f.read()

    # Compute the values for the template placeholders
    from packages import data_tests_functions as dtf

    intrinsic_file = os.path.join(subfolder_path, "intrinsic_properties.yaml")
    if cffs.is_uppsala_structure(subfolder_path) and os.path.isfile(intrinsic_file):
        # New uppsala-data layout: take the derived intrinsic properties
        # (Js, Ms, Ku, Tc) directly from intrinsic_properties.yaml.
        props = dtf.read_intrinsic_properties(intrinsic_file)
        js_0 = props.get("Js", {}).get("value")
        ku_mj = props.get("Ku", {}).get("value")  # MJ/m^3
        tc = props.get("Tc", {}).get("value")
        ms = props.get("Ms", {}).get("value")
        # Js at 300 K from the MC thermal curve, if available
        thermal_csv = os.path.join(subfolder_path, "UppASD", "MC_1", "thermal.csv")
        js_300 = js_0
        if os.path.isfile(thermal_csv):
            thermal_df, _ = dtf.read_thermal_csv(thermal_csv)
            if "Js" in thermal_df.columns:
                js_300 = dtf.interpolate_thermal_property(thermal_df, "Js", 300.0)
        k1_300 = ku_mj * 1e6 if ku_mj is not None else None  # MJ/m^3 -> J/m^3
        print(
            f"From intrinsic_properties.yaml: Js={js_0} T, Ms={ms} kA/m, "
            f"Ku={ku_mj} MJ/m^3, Tc={tc} K"
        )
        print(f"Js at 300 K (interpolated from thermal.csv): {js_300} T")
    else:
        # Old GS/Jij/MC layout: compute the values from the raw output files.
        A_0, A_300, K_300, Js_300, Js_0 = dtf.compute_exchange_and_anisotropy_constants(
            os.path.join(subfolder_path, "MC"),
            dtf.find_line_val_dict(
                os.path.join(subfolder_path, "GS/x/out_last"),
                "Total moment [J=L+S] (mu_B):",
            ),
            dtf.compute_anisotropy_constant(
                os.path.join(subfolder_path, "GS"),
                [d for d in os.listdir(os.path.join(subfolder_path, "GS")) if len(d) == 1],
                dtf.get_unit_cell_volume(os.path.join(subfolder_path, "GS/x/out_last")),
            ),
            dtf.get_unit_cell_volume(os.path.join(subfolder_path, "GS/x/out_last")),
            plot_Js=False,
        )
        k1_300 = K_300
        tc = None
        ms = None

    def fmt(value):
        if value is None:
            return "null"
        if value == 0:
            return "0"
        return str(dtf.round_to_significant_digits(value, 4))

    # Replace placeholders in template
    content = content.replace("$chemical_formula$", subfolder)
    content = content.replace("$filename.h5$", f"{subfolder}.h5")
    content = content.replace("$Js_300$", fmt(js_300))
    content = content.replace("$K1_300$", fmt(k1_300))
    content = content.replace("$Js_0$", fmt(js_0))
    content = content.replace("$Tc$", fmt(tc))
    content = content.replace("$Ms$", fmt(ms))

    with open(yaml_file, "w") as f:
        f.write(content)
    print(f"Created {yaml_file}")
    print(
        f"Js_0: {fmt(js_0)} T, Js_300: {fmt(js_300)} T, K1_300: {fmt(k1_300)} J/m^3, "
        f"Tc: {fmt(tc)} K, Ms: {fmt(ms)} kA/m"
    )

    # Create a .zip containing the .h5, .yaml files, all content of source_dir in a 'datasets' subfolder, and the packages folder
    with zipfile.ZipFile(zip_filename, "w") as zipf:
        # Add .h5 and .yaml files to the root of the zip
        if os.path.exists(h5_file):
            zipf.write(h5_file, arcname=os.path.basename(h5_file))
        if os.path.exists(yaml_file):
            zipf.write(yaml_file, arcname=os.path.basename(yaml_file))

        # Add the Jupyter notebook to the root of the zip
        notebook_file = os.path.join(SCRIPT_DIR, "Evaluation_dataset_UU.ipynb")
        if os.path.exists(notebook_file):
            zipf.write(notebook_file, arcname=os.path.basename(notebook_file))

        # Add all content from the original source_dir to a 'datasets/[source_dir_name]' subfolder in the zip
        for root, dirs, files in os.walk(subfolder_path):
            for file in files:
                file_path = os.path.join(root, file)
                # Calculate the relative path from the subfolder_path
                relative_path = os.path.relpath(file_path, subfolder_path)
                # Add the file to the 'datasets/[source_dir_name]' subfolder in the zip
                arcname = os.path.join("datasets", subfolder, relative_path)
                zipf.write(file_path, arcname=arcname)

        # Add only .py files from the packages subfolder to the zip
        packages_dir = os.path.join(SCRIPT_DIR, "packages")
        if os.path.exists(packages_dir):
            for root, dirs, files in os.walk(packages_dir):
                for file in files:
                    if file.endswith(".py"):
                        file_path = os.path.join(root, file)
                        # Use the relative path from the script directory to maintain folder structure
                        arcname = os.path.relpath(file_path, SCRIPT_DIR)
                        zipf.write(file_path, arcname=arcname)
    print(
        f"Created {zip_filename} containing .h5, .yaml files, Jupyter notebook, all content from {subfolder} in 'datasets/{subfolder}' subfolder, and .py files from packages folder."
    )


if __name__ == "__main__":
    main()
