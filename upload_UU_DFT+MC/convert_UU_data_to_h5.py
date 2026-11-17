"""
convert_UU_data_to_h5.py
Version: 0.0.3
Author: Wilfried Hortschitz
Date: 2025-05-22
Description: Converts folder structure and files to HDF5, including text and CSV files.
Supports both the old layout (GS/Jij/MC) and the new uppsala-data layout
(https://github.com/MaMMoS-project/uppsala-data: RSPt/ + UppASD/ with
intrinsic_properties.yaml and UppASD/MC_*/thermal.csv).
"""

__version__ = "0.0.3"

import os
import h5py
import numpy as np
import argparse
import datetime
import json
from packages import data_tests_functions as dtf
from packages import check_file_and_folder_structure as cffs

parser = argparse.ArgumentParser(
    description="Convert folder structure and files to HDF5."
)
parser.add_argument(
    "--source_dir",
    type=str,
    default="./dataSets/Mn2CrB4",
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
parser.add_argument(
    "--structure",
    type=str,
    choices=["auto", "old", "new"],
    default="auto",
    help="Dataset layout: 'new' = uppsala-data RSPt/+UppASD/, 'old' = GS/Jij/MC, 'auto' = detect (default: auto)",
)
parser.add_argument(
    "--output_file",
    type=str,
    default=None,
    help="Path and name for the output HDF5 file (default: <source_dir_basename>.h5 in current directory)",
)
args = parser.parse_args()

source_dir = args.source_dir
exclude_DOSCAR = args.exclude_DOSCAR
js_minimal = args.js_minimal

if args.structure == "auto":
    structure = "new" if cffs.is_uppsala_structure(source_dir) else "old"
else:
    structure = args.structure
print(f"Dataset structure detected: {structure}")

# Determine output file path
if args.output_file:
    hdf5_file = args.output_file
else:
    hdf5_file = os.path.basename(source_dir.rstrip("/")) + ".h5"

if exclude_DOSCAR:
    hdf5_file = hdf5_file.replace(".h5", "_DOSCAR_excluded.h5")
if js_minimal:
    hdf5_file = hdf5_file.replace(".h5", "_Js_minimal.h5")

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
    creator.attrs["author"] = "Wilfried Hortschitz"
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
                    try:
                        data = np.loadtxt(file_path, delimiter=",")
                        group.create_dataset(file, data=data)
                    except ValueError:
                        # mammos csv v3 (uppsala-data thermal.csv): multi-row header
                        df, meta = dtf.read_thermal_csv(file_path)
                        ds = group.create_dataset(file, data=df.to_numpy())
                        ds.attrs["columns"] = ",".join(meta["columns"])
                        ds.attrs["units"] = json.dumps(meta["units"])
                        ds.attrs["ontology_iri"] = json.dumps(meta["ontology_iris"])
                elif "." not in file:
                    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                        text = f.read()
                    group.create_dataset(file, data=text)

    # Always create the 'results' group and store Js as a dataset, with unit as a subgroup (H5MD style) and link as its attribute
    results_group = h5f.require_group("results")

    # composition_dataset = results_group.create_dataset("composition", data="Co2Fe16Y6")
    # composition_dataset_array = results_group.create_dataset(
    #     "composition_in_array", data=["Co2Fe16Y6"]
    # )

    intrinsic_file = os.path.join(source_dir, "intrinsic_properties.yaml")

    if structure == "new" and os.path.isfile(intrinsic_file):
        # New uppsala-data layout: store the derived intrinsic properties
        # (Js, Ms, Ku, Tc) directly, with their ontology metadata as attrs.
        props = dtf.read_intrinsic_properties(intrinsic_file)
        for name, entry in props.items():
            ds = results_group.create_dataset(name, data=entry["value"])
            if entry.get("unit"):
                ds.attrs["unit"] = entry["unit"]
            if entry.get("ontology_label"):
                ds.attrs["link_in_ontology"] = entry["ontology_label"]
            if entry.get("ontology_iri"):
                ds.attrs["ontology_iri"] = entry["ontology_iri"]
            if entry.get("description"):
                ds.attrs["description"] = entry["description"]
    else:
        mammos_data = dtf.get_mammos_data(source_dir)

        # TODO: add BHmax (what's that?), Lex, Tc (why? it will not be precise),
        # Total_magnetic_moment (we have it as a polarization),
        # and Spin_wave_stiffness_constant to h5 file

        A0_dataset = results_group.create_dataset("A_0", data=mammos_data[0])
        A0_dataset.attrs["link_in_ontology"] = "ExchangeStiffnessConstant"
        A0_dataset.attrs["description"] = "Exchange stiffness constant at 0 K"
        unit_group = results_group.require_group("A_unit")
        # TODO: why A_unit as a separate data entry?
        unit_group.attrs["unit"] = "J/m"
        unit_group.attrs["unit_link_in_ontology"] = "JoulePerMetre"

        A300_dataset = results_group.create_dataset("A_300", data=mammos_data[1])
        A300_dataset.attrs["link_in_ontology"] = "ExchangeStiffnessConstant"
        A300_dataset.attrs["description"] = "Exchange stiffness constant at 300 K"
        unit_group = results_group.require_group("A_unit")
        unit_group.attrs["unit"] = "J/m"
        unit_group.attrs["unit_link_in_ontology"] = "JoulePerMetre"

        K1_300_dataset = results_group.create_dataset("K1_300", data=[mammos_data[2]])
        # TODO: why vector? meaning [mammos_data[2]] vs mammos_data[2]
        K1_300_dataset.attrs["description"] = (
            "Magnetocrystalline anisotropy constant at 300 K"
        )
        K1_300_dataset.attrs["link_in_ontology"] = "MagnetocrystallineAnisotropyConstantK1"
        K1_300_dataset.attrs["unit"] = "J/m^3"
        K1_300_dataset.attrs["unit_link_in_ontology"] = "JoulePerCubicMetre"

        Js_300_dataset = results_group.create_dataset("Js_300", data=[mammos_data[3]])
        # TODO: why vector? meaning [mammos_data[2]] vs mammos_data[2]
        Js_300_dataset.attrs["description"] = "Spontaneous magnetic polarization at 300 K"
        Js_300_dataset.attrs["link_in_ontology"] = "SpontaneousMagneticPolarisation"
        Js_300_dataset.attrs["ontology_iri"] = (
            "https://mammos-project.github.io/MagneticMaterialsOntology/doc/magnetic_material_mammos.html#EMMO_2bb87117-30f9-5b3a-b406-731836a3902f"
        )
        Js_300_dataset.attrs["unit"] = "T"

        Js_300_dataset_as_scalar = results_group.create_dataset(
            "Js_300_scalar", data=mammos_data[3]
        )
        # TODO: why scalar as a separate entry?

        Js_0_dataset = results_group.create_dataset("Js_0", data=[mammos_data[4]])
        # TODO: why vector? meaning [mammos_data[2]] vs mammos_data[2]
        Js_0_dataset.attrs["description"] = "Spontaneous magnetic polarization at 0 K"
        Js_0_dataset.attrs["link_in_ontology"] = "SpontaneousMagneticPolarisation"
        Js_0_dataset.attrs["ontology_iri"] = (
            "https://mammos-project.github.io/MagneticMaterialsOntology/doc/magnetic_material_mammos.html#EMMO_2bb87117-30f9-5b3a-b406-731836a3902f"
        )
        Js_0_dataset.attrs["unit"] = "T"

        Js_0_dataset_as_scalar = results_group.create_dataset(
            "Js_0_scalar", data=mammos_data[4]
        )

if __name__ == "__main__":
    print(f"Convert_UU_data_to_h5.py version {__version__}")
