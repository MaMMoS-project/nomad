import subprocess
import sys
import os


# TODO: allow to rund in nested folders and create a .yaml and .h5 file in the current folder
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
        default="Co2Fe16Y6",
        help="Subfolder name under 'dataSets' directory (default: Co2Fe16Y6)",
    )
    parser.add_argument(
        "--dataSets_dir",
        nargs="?",
        default="./dataSets",
        help="Directory containing the data sets (default: ./dataSets)",
    )
    args = parser.parse_args()
    subfolder = args.source_dir
    dataSets_dir = args.dataSets_dir
    subfolder_path = os.path.join(dataSets_dir, subfolder)
    if not os.path.isdir(subfolder_path):
        print(f"Error: {subfolder} is not a valid directory.")
        sys.exit(1)
    script = "convert_UU_data_to_h5.py"
    h5_file = os.path.join(dataSets_dir, f"{subfolder}.h5")
    yaml_file = os.path.join(dataSets_dir, f"{subfolder}.archive.yaml")
    zip_filename = os.path.join(dataSets_dir, f"{subfolder}_archive.zip")

    # Run the conversion script and create the .h5 file in dataSets
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

    # Create a new archive.yaml file from the template in dataSets
    template = "template.archive.yaml"
    with open(template, "r") as f:
        content = f.read()
    content = content.replace("$chemical_formula$", subfolder)
    content = content.replace("$filename.h5$", f"{subfolder}.h5")
    with open(yaml_file, "w") as f:
        f.write(content)
    print(f"Created {yaml_file}")

    # Create a .zip containing the .h5, .yaml, and .cif files in dataSets
    cif_files = [f for f in os.listdir(subfolder_path) if f.endswith(".cif")]
    with zipfile.ZipFile(zip_filename, "w") as zipf:
        if os.path.exists(h5_file):
            zipf.write(h5_file, arcname=os.path.basename(h5_file))
        if os.path.exists(yaml_file):
            zipf.write(yaml_file, arcname=os.path.basename(yaml_file))
        for cif in cif_files:
            cif_path = os.path.join(subfolder_path, cif)
            zipf.write(cif_path, arcname=cif)
    print(f"Created {zip_filename} containing .h5, .yaml, and .cif files.")


if __name__ == "__main__":
    main()
