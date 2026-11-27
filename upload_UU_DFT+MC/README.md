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

The derived files (`intrinsic_properties.yaml`, `UppASD/MC_*/thermal.csv`) in
each dataset are **generated and validated with
[mammos-parser](https://github.com/MaMMoS-project/mammos-parser)**
(`mammos-parser uppsala-data generate-derived-files <dir>`), not by hand. The
uploaded `<dataset>.h5` is a mammos-entity `EntityCollection` — the shipped
`Evaluation_dataset_UU.ipynb` shows how to read it (`me.from_hdf5`,
`me.from_yaml`, `me.from_csv`).

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
  `<dataset>_structure.cif` (raw data / provenance — NOMAD never reads
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
- The schema uses the ELN `Chemical` base section (as in the original/initial
  UU template) + `EntryData` + `HDF5Normalizer`. The `Chemical` base provides
  the `chemical_formula` quantity (filled from the folder name when
  `--cif off`); its `normalize()` populates `results.material` server-side
  (elements + all formula formats), which fills the Overview Material section
  and makes element search work. We deliberately do NOT write `results`
  directly into the `.archive.yaml`. The section additionally declares an
  `elements` **string** quantity (comma-separated, e.g. `Al, C, Fe`), filled
  from the same folder-name formula via `elements_from_formula()` — no
  conflict with `chemical_formula` (the ELN base only reads
  `chemical_formula`; `elements` is a plain data-section quantity). It must
  NOT be a list (`shape: ['*']`): nomad-lab 2.0.1 `HDF5Normalizer.normalize()`
  runs `re.match(h5_re, value)` on *every* str-typed quantity of the section
  and crashes on a list value (`TypeError: expected string or bytes-like
  object, got 'list'`), which aborts the whole section normalize — no h5
  fill, and the `chemical_formula` → `results.material` population never
  runs (found by local validation on 2.0.1, would also fail server-side).
  Both `elements` and `chemical_formula` are omitted when `--cif on`. Trade-off: ELN
  base fields (datetime, chemical_formula, ...) show as editable input fields
  in the NOMAD 2.0 UI — accepted because this matches the initial datasets.
- The full temperature-dependent curves are exported as a repeating
  `thermal_properties` subsection — one entry per `UppASD/MC_*/thermal.csv`
  (mammos csv v3), named `MC_1`, `MC_2`, …, described by the csv header
  comment ("Temperature-dependent quantities computed with UppASD"), with
  array quantities `T, Ms, Js, E, Cv, chi, U_L` (ME names/units as-is, extra
  csv columns are added dynamically) and a per-column `mammos_entities`
  record (ontology label/IRI/description from the csv header). The arrays are
  written inline — an `hdf5.path` annotation on the shared quantity
  definition would fill only the *last* repeating instance server-side
  (`create_custom_mapping` uses `m_get_sub_section(sub_section, -1)`), which
  would map every MC run to the same curve; the `.h5` still carries the
  curves under `/thermal_MC_<n>/<column>` as raw data. The intrinsic
  properties stay grouped under the entry description ("Intrinsic magnetic
  properties computed with RSPt and UppASD." from `intrinsic_properties.yaml`).
- With `--cif on` (run section embedded), `chemical_formula` is omitted and
  the composition is derived from the embedded structure server-side; keeping
  both caused the "elements is already defined" error (see above).
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
# (<dataset>_structure.cif at the zip root and under raw_datasets/<dataset>/).
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

# QUICK TEST rebuild: only ~3 datasets (QUICK_DATASETS in the script:
# Al2Fe6C, Fe16N2_exp_str, Co3Mn2Ge — plain name, '_suffix' name, multi-MC),
# ~40 s instead of several minutes for all 23. Run without flags to rebuild ALL.
python batch_process_datasets.py --quick
```

Notes:
- **Always regenerate ALL datasets before an upload** (`python batch_process_datasets.py` without flags). Zips already present in `uppsala-data/uploads/` can be stale — produced with an older archive template (e.g. `elements` as a list, which crashes the 2.0.1 `HDF5Normalizer`) or an older zip layout (cif named `cif_structure_<name>.cif`, raw tree under `datasets/` instead of `<name>_structure.cif` / `raw_datasets/`).
- The batch script runs each dataset with the default `--cif off` (it does not forward a `--cif` flag); run the single-dataset command with `--cif auto|on` for datasets that need the embedded structure preview.
- Equivalent manual alternative (shell loop):
  `for d in uppsala-data/data/*/; do python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir "$(basename "$d")"; done`

### Local validation (optional, Level 1)

Each generated `.archive.yaml` + `.h5` pair can be validated **locally** with
the real NOMAD parser + normalizer chain (no network, no credentials) before
uploading:

```bash
# one-time: create the opt-in nomad-test environment (nomad-lab pinned to the
# central instance version, currently 2.0.1)
pixi install -e nomad-test

# validate while processing ONE dataset (stages a temp dir with just this
# dataset's .archive.yaml + .h5 [+ cif], exits non-zero on failure,
# writes <dataset>.validated.json next to the archive):
python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir Fe12Co4N2 --validate

# or validate an existing archive directly (data file must sit next to it):
pixi run -e nomad-test validate uppsala-data/data/Fe12Co4N2.archive.yaml
```

`scripts/validate_archive_locally.py` stages the archive's parent directory
as a NOMAD staging upload, runs `parsers/archive` + all normalizers of the
pinned nomad-lab, and fails if any `hdf5`-annotated quantity is unfilled, any
figure trace is empty, or (UU-specific, `--cif off` archives) the ELN
`chemical_formula` did not populate `results.material.elements`. The full
normalized archive is dumped as JSON for inspection.

Caveat: local OK means "identical parse/normalize result expected on the
server" — it does NOT cover server-only steps (msgpack `write_archive`,
Mongo/Elasticsearch indexing, the server-only `SystemNormalizer`, GUI
rendering). The final check remains the server-side `process_status` after
upload. The nomad-test env deliberately does NOT inherit the runtime conda
deps (conda-forge `pandas` pins `python-dateutil==4.0.0`, conflicting with
nomad-lab's `s3fs`/`aiobotocore` requirement `python-dateutil<3`).

### Uploading to NOMAD via the API (Level 2)

Upload the zips' content into an existing NOMAD project via the REST API.
**Use the backend ingress for uploads**:
`https://nomad-lab.eu/backend/api/v1` — the public edge
`https://nomad-lab.eu/prod/v1/api/v1` caps request bodies at ~1 MB (HTTP 413
above that, verified 2026-11-25), while the backend applies NOMAD's documented
35 GB upload limit. Minimal example using the nomad-lab package from the
`nomad-test` env (token only from env `NOMAD_CLIENT_ACCESS_TOKEN` or
`~/.config/to-nomad/token`, never in repo files/chat/CLI):

```bash
NOMAD_CLIENT_URL="https://nomad-lab.eu/backend/api" \
NOMAD_CLIENT_ACCESS_TOKEN="$(cat ~/.config/to-nomad/token)" \
./.pixi/envs/nomad-test/bin/python -c "
from nomad.client import api
files = [('file', open(p, 'rb')) for p in ['<name>.archive.yaml', '<name>.h5', 'README.md']]
print(api.put('uploads/<PROJECT_ID>/raw/', files=files, auth=api.Auth()).status_code)
"
```

See [doc/uppsala-data-nomad-upload.md](doc/uppsala-data-nomad-upload.md) section
"Uploading to NOMAD via the API" for the full route/parameter reference
(`wait_for_processing`, processing polling, rawdir pagination, 504-on-slow-link
caveat). Any real upload requires an explicit user request naming files +
project id with `[Approved]` in the same message (AGENTS.md).

## Output Files

After processing, the created output files will be stored in the `uploads/` directory as `.zip` files ready for uploading to NOMAD:

```
uploads/
├── Fe12Co4N2_archive.zip
├── Co2Fe2H4_archive.zip
└── [dataset_name]_archive.zip
```

Each zip file contains:
- `.h5` file with processed data (a mammos-entity `EntityCollection`)
- `.archive.yaml` file with metadata and computed magnetic properties
- `README.md` — general upload README (no dataset-specific names/values; points
  to the entry's Overview/Data tabs and the raw files, where they are found).
  The NOMAD GUI shows exactly this file in the upload overview (it fetches
  `raw/README.md`); it is byte-identical in every zip, so multi-zip uploads
  are safe
- `Evaluation_dataset_UU.ipynb` Jupyter notebook showing how to read the data
  with `mammos_entity`; it uses local paths only (the `.h5` / `raw_datasets/` next
  to the notebook), so the same notebook works in every unpacked upload
- `raw_datasets/[dataset_name]/` folder with original raw data
- `packages/` folder with Python analysis tools

These zip files are self-contained and ready for upload to the NOMAD repository.

