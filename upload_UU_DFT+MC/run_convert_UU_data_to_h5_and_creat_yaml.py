import subprocess
import sys
import os


# DONE-TODO: allow to rund in nested folders and create a .yaml and .h5 file in the current folder
# TODO: check if the folder-name is a valid chemical formula (Co2Fe16Y6 vs. Fe3-xZnxY(x=0.22)) and make sure the formula is provided elsewhere


# Usage: python run_convert_UU_data_to_h5.py <subfolder_name>
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
        help="Subfolder name under 'datasets' directory (default: Fe12Co4N2)",
    )
    parser.add_argument(
        "--datasets_dir",
        nargs="?",
        default="./datasets",
        help="Directory containing the data sets (default: ./datasets)",
    )
    args = parser.parse_args()
    subfolder = args.source_dir
    datasets_dir = args.datasets_dir
    subfolder_path = os.path.join(datasets_dir, subfolder)
    if not os.path.isdir(subfolder_path):
        print(f"Error: {subfolder_path} is not a valid directory.")
        sys.exit(1)
    script = "convert_UU_data_to_h5.py"
    h5_file = os.path.join(datasets_dir, f"{subfolder}.h5")
    yaml_file = os.path.join(datasets_dir, f"{subfolder}.archive.yaml")

    # Create uploads directory at the same level as datasets
    uploads_dir = os.path.join(os.path.dirname(datasets_dir), "uploads")
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

    # Create a new archive.yaml file from the template in datasets
    template = "UU_template.archive.yaml"
    with open(template, "r") as f:
        content = f.read()

    # Compute Js_300 and K1_300 values using the same logic as in Evaluation_dataset_UU.py
    from packages import data_tests_functions as dtf

    # Get the computed values from the data processing
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

    # Replace placeholders in template
    content = content.replace("$chemical_formula$", subfolder)
    content = content.replace("$filename.h5$", f"{subfolder}.h5")
    content = content.replace(
        "$Js_300$", str(dtf.round_to_significant_digits(Js_300, 4))
    )
    content = content.replace(
        "$K1_300$", str(dtf.round_to_significant_digits(K_300, 4))
    )
    content = content.replace("$Js_0$", str(dtf.round_to_significant_digits(Js_0, 4)))

    with open(yaml_file, "w") as f:
        f.write(content)
    print(f"Created {yaml_file}")
    print(f"Computed Js_0: {dtf.round_to_significant_digits(Js_0, 4)} T")
    print(f"Computed Js_300: {dtf.round_to_significant_digits(Js_300, 4)} T")
    print(f"Computed K1_300: {dtf.round_to_significant_digits(K_300, 4)} J/m³")

    # Create a .zip containing the .h5, .yaml files, all content of source_dir in a 'datasets' subfolder, and the packages folder
    with zipfile.ZipFile(zip_filename, "w") as zipf:
        # Add .h5 and .yaml files to the root of the zip
        if os.path.exists(h5_file):
            zipf.write(h5_file, arcname=os.path.basename(h5_file))
        if os.path.exists(yaml_file):
            zipf.write(yaml_file, arcname=os.path.basename(yaml_file))

        # Add the Jupyter notebook to the root of the zip
        notebook_file = "Evaluation_dataset_UU.ipynb"
        if os.path.exists(notebook_file):
            zipf.write(notebook_file, arcname=notebook_file)

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
        packages_dir = "packages"
        if os.path.exists(packages_dir):
            for root, dirs, files in os.walk(packages_dir):
                for file in files:
                    if file.endswith(".py"):
                        file_path = os.path.join(root, file)
                        # Use the relative path from current directory to maintain folder structure
                        arcname = file_path
                        zipf.write(file_path, arcname=arcname)
    print(
        f"Created {zip_filename} containing .h5, .yaml files, Jupyter notebook, all content from {subfolder} in 'datasets/{subfolder}' subfolder, and .py files from packages folder."
    )


if __name__ == "__main__":
    main()
