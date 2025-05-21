"""
convert_UU_data_to_h5.py
Version: 0.0.1
Author: [Your Name]
Date: 2025-05-14
Description: Converts folder structure and files to HDF5, including text and CSV files.
"""

__version__ = "0.0.1"

import os
import h5py
import numpy as np
import argparse
import datetime
import data_tests_functions as dtf

parser = argparse.ArgumentParser(
    description="Convert folder structure and files to HDF5."
)
parser.add_argument(
    "--source_dir",
    type=str,
    default="Co2Fe16Y6",
    help="Source directory to convert (default: Co2Fe16Y6)",
)
parser.add_argument(
    "--exclude_DOSCAR",
    action="store_true",
    help="Exclude DOSCAR files from being added to the HDF5 file (default: include DOSCAR)",
)
parser.add_argument(
    "--js_minimal",
    action="store_true",
    help="Export only Js to the HDF5 file (default: export full structure)",
)
args = parser.parse_args()

source_dir = args.source_dir
exclude_DOSCAR = args.exclude_DOSCAR
js_minimal = args.js_minimal
hdf5_file = "output_" + os.path.basename(source_dir.rstrip("/"))
if exclude_DOSCAR:
    hdf5_file += "_DOSCAR_excluded"
if js_minimal:
    hdf5_file += "_Js_minimal"
hdf5_file += ".h5"

# Remove the HDF5 file if it already exists
if os.path.exists(hdf5_file):
    os.remove(hdf5_file)

with h5py.File(hdf5_file, "w") as h5f:
    # Save script version in HDF5 file attributes
    h5f.attrs["conversion_script_version"] = __version__
    h5f.attrs["conversion_script_filename"] = os.path.basename(__file__)
    # H5MD root attribute
    h5f.attrs["H5MD_version"] = "1.0.0"
    # Creator group for H5MD compliance
    creator = h5f.require_group("creator")
    creator.attrs["name"] = "convert_UU_data_to_h5.py"
    creator.attrs["version"] = __version__
    creator.attrs["author"] = "[Wilfried Hortschitz]"
    creator.attrs["date"] = datetime.datetime.now().isoformat()

    # Only export raw_data if not js_minimal
    if not js_minimal:
        raw_data_group = h5f.require_group("raw_data")
        for root, dirs, files in os.walk(source_dir):
            rel_path = os.path.relpath(root, source_dir)
            group = raw_data_group.require_group(rel_path)
            for file in files:
                file_path = os.path.join(root, file)
                if exclude_DOSCAR and file == "DOSCAR":
                    print(f"Excluding DOSCAR file: {file_path}")
                    continue
                if file.endswith(".csv"):
                    data = np.loadtxt(file_path, delimiter=",")
                    group.create_dataset(file, data=data)
                elif "." not in file:
                    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                        text = f.read()
                    group.create_dataset(file, data=text)

    # Always create the 'results' group and store Js as a dataset, with unit as a subgroup (H5MD style) and link as its attribute
    results_group = h5f.require_group("results")

    composition_dataset = results_group.create_dataset(
        "composition", data=["Co2Fe16Y6"]
    )

    mammos_data = dtf.get_mammos_data(source_dir)

    js_dataset = results_group.create_dataset("A_0", data=mammos_data[0])
    unit_group = results_group.require_group("A_unit")
    unit_group.attrs["unit"] = "T"
    js_dataset.attrs["link_in_ontology"] = "ExchangeStiffnessConstantA"

    js_dataset = results_group.create_dataset("A_300", data=mammos_data[1])
    unit_group = results_group.require_group("A_unit")
    unit_group.attrs["unit"] = "T"
    js_dataset.attrs["link_in_ontology"] = "ExchangeStiffnessConstantA"

    js_dataset = results_group.create_dataset("K1_300", data=mammos_data[2])
    unit_group = results_group.require_group("Js_unit")
    unit_group.attrs["unit"] = "T"
    js_dataset.attrs["link_in_ontology"] = "SpontaneousMagneticPolarisation"

    js_dataset = results_group.create_dataset("Js_300", data=mammos_data[3])
    unit_group = results_group.require_group("K1_unit")
    unit_group.attrs["unit"] = "J/m^3"
    unit_group.attrs["unit_link_in_ontology"] = "JoulePerCubicMetre"

    js_dataset.attrs["link_in_ontology"] = "MagnetocrystallineAnisotropyConstantK1"
    js_dataset.attrs["ontology_iri"] = (
        "https://mammos-project.github.io/MagneticMaterialsOntology/doc/magnetic_material_mammos.html#EMMO_2bb87117-30f9-5b3a-b406-731836a3902f"
    )


if __name__ == "__main__":
    print(f"Convert_UU_data_to_h5.py version {__version__}")
