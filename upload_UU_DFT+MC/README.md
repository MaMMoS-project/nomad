# Code used to create a upload on Nomad for data from UU

## Data source (new structure)

UU simulation data is available from
[github.com/MaMMoS-project/uppsala-data](https://github.com/MaMMoS-project/uppsala-data).
Clone it into this project directory (the clone is gitignored):

```bash
git clone https://github.com/MaMMoS-project/uppsala-data.git
```

The scripts default to `./uppsala-data/data` as the datasets directory. Each
dataset subfolder there (named by a reduced chemical formula, e.g. `Fe12Co4N2`)
follows the schema in `uppsala-data/dataset-schema.yaml`:

- `metadata.yaml`, `intrinsic_properties.yaml` (derived `Js`, `Ms`, `Ku`, `Tc`
  with ontology IRIs), `structure.cif` (required); `DOSCAR`, `README.md` (optional)
- `RSPt/` with `common_input/`, ground-state dirs `gs_x/`, `gs_z/` (optional
  `gs_y/`; each with `data`, `out_last`, and `hist` XOR `out_MF`), and `Jij/`
- `UppASD/MC_1/` (optional `MC_2/`, `MC_3/` for Binder-cumulant Tc) with `jfile`,
  `momfile`, `posfile`, `inpsd.dat`, `thermal.csv` (mammos csv v3 with ontology
  header; columns T, Ms, Js, E, Cv, chi, U_L), `thermal.dat`

Datasets in this layout are validated with
`mammos-parser uppsala-data validate-dataset <dataset_dir>` (installed via
pixi; see `pixi.toml`) plus an internal structural check
(`packages/check_file_and_folder_structure.py: check_uppsala_structure`). The
`.archive.yaml` placeholders are filled directly from
`intrinsic_properties.yaml` (`Js` -> `Js_0`, `Ku` -> `K1_300`, `Tc`, `Ms`), and
`Js_300` is interpolated from `UppASD/MC_1/thermal.csv` at 300 K.

### Change from the prior static solution

Earlier versions mapped the mammos-entity values onto fixed template
quantities: `Ku` was written to `K1_300` (labelled "Magnetocrystalline
anisotropy constant at 300 K") with a MJ/m^3 -> J/m^3 conversion, and `Js`
to `Js_0`. This was wrong: the ME `Ku` (ontology label
`UniaxialAnisotropyConstant`) is not "K1 at 300 K". Now the mammos-entity
**variable names, ontology labels and units are used as-is** — the archive
quantities are `Js`, `Ku` (MJ/m^3), `Tc`, `Ms` (plus `Js_300`, the Js curve
interpolated at 300 K from `thermal.csv`), and quantity/entry descriptions
are filled at runtime from `intrinsic_properties.yaml` / the EMMO ontology
instead of being hardcoded in the template. In the `.h5` file the new-layout
path already stored the ME names and the per-entry `description` from the
yaml (e.g. `Tc` -> "Tc from peak in specific heat"); only the old-layout path
keeps hardcoded descriptions.

The generated archive additionally carries ontology metadata, following
MaMMoS-project/to-nomad (branch `feature/subcollection-subsections`):

- a repeating `mammos_entities` subsection lists each quantity with its
  `name`, `value_key`, `ontology_label`, `ontology_iri` (MaMMoS/EMMO
  magnetic-materials ontology), `unit` and `description` — e.g. `Tc` is
  stored as `CurieTemperature` with IRI
  `https://w3id.org/emmo#EMMO_6b5af5a8_a2d8_4353_a1d6_54c9f778343d`;
- With `--cif auto/on`, `structure.cif` is parsed locally with ASE and
  embedded as an inline top-level `run` section (plus
  `metadata: {domain: dft}`), so NOMAD shows a 3-D structure preview in the
  entry overview; the `.cif` itself is also copied to the zip root as
  `cif_structure_<dataset>.cif` (raw data / provenance — NOMAD never reads
  the cif at upload time). When the run section is embedded,
  `chemical_formula` is **omitted** from the archive: keeping both makes the
  `Chemical` base section derive `elements` from the formula, which collides
  with the system-derived composition of the run section (server-side
  `ResultsNormalizer` error: "Could not populate compositional data as
  'elements' is already defined", observed on a Fe12Co4N2 upload). With the
  default `--cif off`, no run section is embedded and `chemical_formula`
  (folder name without `_comment` suffix) is used instead.
- entry `name` / `short_name` are filled (`UU-DFT+MC-<dataset>` /
  `UU-<dataset>`) so they no longer show "no value" in Data -> Quantities.
- The schema uses plain `EntryData` + `HDF5Normalizer` base sections (no ELN
  `Chemical` base), so the new NOMAD UI shows entry fields (datetime,
  chemical_formula, ...) read-only instead of as input fields, and
  `chemical_formula` is declared as our own quantity (filled from the folder
  name when `--cif off`).
- The generated `.h5` is linked via `data_file` (a `str` quantity that MUST be
  declared in the section — NOMAD's `HDF5Normalizer` opens any string
  quantity whose value ends in `.h5` as `<upload>/raw/<value>`; without the
  declared quantity the value is dropped at parse time and nothing is
  mapped, which showed up as "no value" for all h5-linked quantities).
  Quantities present in the h5 `results` group (`Js`, `Ku`, `Tc`, `Ms`; old
  layout: `Js_0`, `K1_300`, `Js_300`) get `m_annotations.hdf5.path:
  /results/<name>` and no inline value — the server-side `HDF5Normalizer`
  fills them from the uploaded h5 (same mechanism as to-nomad's
  `hdf_reference_mode`). Quantities not in the h5 (e.g. `Js_300` for new
  datasets) are written inline.

## Data source (old structure, legacy)

Alternatively, **create a folder named `datasets`** in your project directory.
Into the folder called 'datasets' one should put the simulation dataset to be
processed. Each dataset subfolder should have the following structure:

- A `.cif` file containing the crystal structure information.
- A `.README` file with details about the dataset and the expected results.
- A `GS` folder, which contains files related to ground state calculations (e.g., output files).
- A `JiJ` folder, which includes files for exchange interaction calculations (e.g., Jij calculation outputs).
- An `MC` folder, which holds files for the Monte Carlo simulations.

This organization ensures that all relevant data and documentation for each chemical composition are kept together and easy to locate.

## Usage

### Processing Single Datasets

To process a single dataset, use the `run_convert_UU_data_to_h5_and_creat_yaml.py` script:

```bash
# Process a specific dataset (e.g., Fe12Co4N2); datasets dir defaults to
# ./uppsala-data/data if present, else ./datasets
python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir Fe12Co4N2

# Process with custom datasets directory
python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir Fe12Co4N2 --datasets_dir ./uppsala-data/data
```

Note: `--datasets_dir` only locates the parent folder containing the dataset
subfolders. Each run of `run_convert_UU_data_to_h5_and_creat_yaml.py` processes
exactly ONE dataset — the one named by `--source_dir`. Giving `--datasets_dir`
alone does NOT loop over all datasets in that folder. For that, use the batch
script below.

```bash
# Control the cif embedding in the .archive.yaml (3-D structure preview).
# Default is OFF: no run section is embedded and chemical_formula (from the
# folder name, '_comment' suffix stripped) is used for the Material section.
# Use 'auto' to embed structure.cif when present, 'on' to require it.
# Note: with the run section embedded, chemical_formula must be omitted
# (server-side collision). Independently of --cif, the .cif is ALWAYS
# included in the upload zip as raw data when present
# (cif_structure_<dataset>.cif at the zip root and under datasets/<dataset>/).
python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir Fe12Co4N2 --cif auto
```

### Processing Multiple Datasets

To process multiple datasets at once, use the `batch_process_datasets.py` script. It iterates over all subfolders of the datasets directory and calls `run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir <name> --datasets_dir <dir>` once per dataset, continuing past failures and printing a summary at the end:

```bash
# Process ALL subfolders of the datasets directory
# (default: ./uppsala-data/data if present, else ./datasets)
python batch_process_datasets.py

# Batch with a custom datasets directory
python batch_process_datasets.py --datasets_dir ./uppsala-data/data

# Process only selected datasets (comma-separated names)
python batch_process_datasets.py --only Fe12Co4N2,Co2Fe2H4
```

Notes:
- The batch script runs each dataset with the default `--cif off` (it does not forward a `--cif` flag); run the single-dataset command with `--cif auto|on` for datasets that need the embedded structure preview.
- Equivalent manual alternative (shell loop):
  `for d in uppsala-data/data/*/; do python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir "$(basename "$d")"; done`

## Output Files

After processing, the created output files will be stored in the `uploads/` directory as `.zip` files ready for uploading to NOMAD:

```
uploads/
├── Fe12Co4N2_archive.zip
├── Co2Fe2H4_archive.zip
└── [dataset_name]_archive.zip
```

Each zip file contains:
- `.h5` file with processed data
- `.archive.yaml` file with metadata and computed magnetic properties
- `Evaluation_dataset_UU.ipynb` Jupyter notebook for data analysis
- `datasets/[dataset_name]/` folder with original raw data
- `packages/` folder with Python analysis tools

These zip files are self-contained and ready for upload to the NOMAD repository.

