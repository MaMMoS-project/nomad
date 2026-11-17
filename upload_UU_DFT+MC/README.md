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

The generated archive additionally carries ontology metadata, following
MaMMoS-project/to-nomad (branch `feature/subcollection-subsections`):

- a repeating `mammos_entities` subsection lists each quantity with its
  `name`, `value_key`, `ontology_label`, `ontology_iri` (MaMMoS/EMMO
  magnetic-materials ontology), `unit` and `description` — e.g. `Tc` is
  stored as `CurieTemperature` with IRI
  `https://w3id.org/emmo#EMMO_6b5af5a8_a2d8_4353_a1d6_54c9f778343d`;
- `structure.cif` is parsed locally with ASE and embedded as an inline
  top-level `run` section (plus `metadata: {domain: dft}`), so NOMAD shows a
  3-D structure preview in the entry overview; the `.cif` itself is also
  copied to the zip root as `cif_structure_<dataset>.cif` (raw data /
  provenance — NOMAD never reads the cif at upload time).

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

### Processing Multiple Datasets

To process multiple datasets at once, use the `batch_process_datasets.py` script:

```bash
# Process all subfolders of the datasets directory
python batch_process_datasets.py

# Process only selected datasets
python batch_process_datasets.py --only Fe12Co4N2,Co2Fe2H4
```

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

