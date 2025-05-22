import subprocess
import sys
import os
import shutil


# Usage: python run_convert_UU_data_to_h5.py <subfolder_name>
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Run convert_UU_data_to_h5.py for a subfolder and create a matching archive.yaml file."
    )
    parser.add_argument(
        "--source_dir",
        nargs="?",
        default="Co2Fe16Y6",
        help="Subfolder name (default: Co2Fe16Y6)",
    )
    args = parser.parse_args()
    subfolder = args.source_dir
    if not os.path.isdir(subfolder):
        print(f"Error: {subfolder} is not a valid directory.")
        sys.exit(1)
    script = "convert_UU_data_to_h5.py"
    cmd = [sys.executable, script, "--source_dir", subfolder]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Create a new archive.yaml file from the template
    template = "template.archive.yaml"
    new_archive = f"{subfolder}.archive.yaml"
    with open(template, "r") as f:
        content = f.read()
    # Replace all occurrences of e.g. 'Co2Fe16Y6' with the subfolder name
    content = content.replace("$chemical_formula$", subfolder)
    # Replace all occurrences of e.g. 'Co2Fe16Y6.h5' with '<subfolder>.h5'
    content = content.replace(f"$filename.h5$", f"{subfolder}.h5")
    # Save the new archive.yaml
    with open(new_archive, "w") as f:
        f.write(content)
    print(f"Created {new_archive}")


if __name__ == "__main__":
    main()
