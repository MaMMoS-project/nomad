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
    parser.add_argument(
        "--cif",
        choices=["auto", "on", "off"],
        default="off",
        help="Embed structure.cif as inline NOMAD run section (3-D structure preview): "
        "'off' (default) = never embed it, chemical_formula from the folder name is used; "
        "'auto' = embed if structure.cif exists; 'on' = require it",
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

    # Whether to embed the cif as an inline run section (structure preview)
    cif_path = os.path.join(subfolder_path, "structure.cif")
    use_cif = args.cif == "on" or (args.cif == "auto" and os.path.isfile(cif_path))
    if args.cif == "on" and not os.path.isfile(cif_path):
        print(f"Error: --cif on but no structure.cif in {subfolder_path}")
        sys.exit(1)
    if not use_cif:
        print(
            "cif not embedded (--cif off or no structure.cif): "
            "chemical_formula from the folder name will be used instead."
        )

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

    # Compute the values for the archive fields
    from packages import data_tests_functions as dtf
    from packages import nomad_archive_extras as nae

    intrinsic_file = os.path.join(subfolder_path, "intrinsic_properties.yaml")
    props = None
    if cffs.is_uppsala_structure(subfolder_path) and os.path.isfile(intrinsic_file):
        # New uppsala-data layout: take the derived intrinsic properties
        # (Js, Ms, Ku, Tc) directly from intrinsic_properties.yaml, using the
        # mammos-entity names/labels/units as-is (Ku stays in MJ/m^3).
        props = dtf.read_intrinsic_properties(intrinsic_file)
        js = props.get("Js", {}).get("value")
        ku = props.get("Ku", {}).get("value")  # MJ/m^3
        tc = props.get("Tc", {}).get("value")
        ms = props.get("Ms", {}).get("value")
        # Js at 300 K from the MC thermal curve, if available
        thermal_csv = os.path.join(subfolder_path, "UppASD", "MC_1", "thermal.csv")
        js_300 = js
        if os.path.isfile(thermal_csv):
            thermal_df, _ = dtf.read_thermal_csv(thermal_csv)
            if "Js" in thermal_df.columns:
                js_300 = dtf.interpolate_thermal_property(thermal_df, "Js", 300.0)
        print(
            f"From intrinsic_properties.yaml: Js={js} T, Ms={ms} kA/m, "
            f"Ku={ku} MJ/m^3, Tc={tc} K"
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
        # Old layout computes K_300 in J/m^3; the template quantity Ku uses MJ/m^3
        js = Js_0
        ku = K_300 / 1e6 if K_300 is not None else None
        tc = None
        ms = None

    def rnd(value):
        if value is None or value == 0:
            return value
        return dtf.round_to_significant_digits(value, 4)

    # Create the archive.yaml from the template (structured fill)
    import yaml

    template = os.path.join(SCRIPT_DIR, "UU_template.archive.yaml")
    with open(template, "r") as f:
        doc = yaml.safe_load(f)
    data = doc["data"]
    if not use_cif:
        # No embedded structure: provide the chemical formula so the NOMAD
        # Chemical base section populates formula/elements (Material section).
        formula = nae.chemical_formula_from_name(subfolder)
        data["chemical_formula"] = formula
        import re

        if not re.fullmatch(r"([A-Z][a-z]?[0-9.]*)+", formula):
            print(
                f"Warning: folder name '{subfolder}' -> formula '{formula}' does not "
                "look like a plain chemical formula; NOMAD may not parse it."
            )
    data["name"] = f"UU-DFT+MC-{subfolder}"
    data["short_name"] = f"UU-{subfolder}"
    data["Js"] = rnd(js)
    data["Js_300"] = rnd(js_300)
    data["Ku"] = rnd(ku)
    data["Tc"] = rnd(tc)
    data["Ms"] = rnd(ms)
    if props is not None:
        data["mammos_entities"] = nae.mammos_entity_records(props)
    else:
        data["mammos_entities"] = nae.records_for_computed(js, js_300, ku)

    # Fill quantity descriptions and the entry description from the mammos
    # files / ontology (not hardcoded in the template).
    quantities = doc["definitions"]["sections"]["Simulation_data_from_UU"]["quantities"]
    if props is not None:
        for q in ["Js", "Ku", "Tc", "Ms"]:
            entry = props.get(q, {})
            quantities[q]["description"] = nae.ontology_description(
                entry.get("ontology_label") or "",
                entry.get("ontology_iri") or "",
                entry.get("description") or "",
            )
        js_entry = props.get("Js", {})
        quantities["Js_300"]["description"] = nae.ontology_description(
            js_entry.get("ontology_label") or "",
            js_entry.get("ontology_iri") or "",
            "Spontaneous magnetic polarization at 300 K "
            "(interpolated from UppASD/MC_1/thermal.csv)",
        )
        yaml_desc = nae.mammos_yaml_description(intrinsic_file)
        if yaml_desc:
            data["description"] = f"{data.get('description') or ''}\n{yaml_desc}".strip()
    else:
        for r in data["mammos_entities"]:
            if r["value_key"] in quantities:
                quantities[r["value_key"]]["description"] = nae.ontology_description(
                    r["ontology_label"], r["ontology_iri"], r["description"]
                )

    # Embed the crystal structure as an inline top-level run section (parsed
    # locally with ASE), so NOMAD shows a 3-D structure preview in the entry
    # overview (server-side OptimadeNormalizer via metadata.domain=dft).
    # The .cif itself is shipped in the zip as raw data (see below).
    if use_cif:
        doc["run"] = [nae.structure_run_section(cif_path)]
        doc["metadata"] = {"domain": "dft"}
        # Do NOT keep chemical_formula together with the inline run section:
        # the Chemical base section would derive `elements` from the formula,
        # colliding with the system-derived composition of the run section
        # (server-side ResultsNormalizer error: "Could not populate
        # compositional data as 'elements' is already defined"). The
        # composition is derived from the embedded structure instead.
        data.pop("chemical_formula", None)

    with open(yaml_file, "w") as f:
        yaml.safe_dump(doc, f, sort_keys=False, allow_unicode=True)
    print(f"Created {yaml_file}")
    print(
        f"Js: {rnd(js)} T, Js_300: {rnd(js_300)} T, Ku: {rnd(ku)} MJ/m^3, "
        f"Tc: {rnd(tc)} K, Ms: {rnd(ms)} kA/m"
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

        # Add the crystal structure file to the root of the zip as raw data
        # (the parsed structure is embedded in the archive.yaml run section)
        cif_path = os.path.join(subfolder_path, "structure.cif")
        if os.path.isfile(cif_path):
            zipf.write(cif_path, arcname=f"cif_structure_{subfolder}.cif")

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
