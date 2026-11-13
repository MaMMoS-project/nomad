# NEEL high-throughput (HT) HDF5 datasets: structure and versions

Findings from inspecting the NEEL datasets in `datasets/` with h5py, and from
comparing against the upstream MaMMoS tools that read/write this format
(`MaMMoS-project/Advanced_Data_Visualization`, `MaMMoS-project/ht-data-parser`,
`MaMMoS-project/DaHU` — all authored at Institut Néel, CNRS).

Datasets currently in `datasets/`:

| File | root attr | root writer attr | measurement group writers | inner layout |
|---|---|---|---|---|
| `NdCeFeB_2-5.hdf5` | `HT_class = 'HTroot'` (no root `HT_type`) | none | all `'0.1 beta'` (`edx_writer`, `moke_writer`, `esrf_writer`) | flat (oldest) |
| `SmFeV_3378.hdf5` | `HT_type = 'library'` | `library_writer = 0.1` | `edx_writer '0.1 beta'`, `moke_writer '0.2'`, `esrf_writer '0.2'`, `profil_writer '0.3'` | flat (v0.x) |
| `NdCeFeB_4_6.hdf5` | `HT_type = 'dataset'` | `dataset_writer = 0.1` | all `0.9` (`edx_writer`, `moke_writer`, `esrf_writer`, `profil_writer`) | nested `positions/` (v0.9) |
| `NdCeFeB_3_13.hdf5` | `HT_type = 'library'` | `library_writer = 0.1` | all `0.9` (incl. new `sem_writer '0.9'`) | nested `positions/` (v0.9) |

Key point: the root `HT_type` (`library` vs `dataset`) describes the file's role
(library = experimentalist-centered collection of scans; dataset = sample-based
reduced file, per the MaMMoS HT paper and ht-data-parser docs), but it does NOT
predict the inner layout — `NdCeFeB_3_13.hdf5` is a `library` file with the new
v0.9 nested layout. The reliable layout marker is the per-group writer version
(`0.9`) or the presence of a `positions/` sub-group.

## Root level

- File attributes: either `HT_type` (`'library'` or `'dataset'`) plus
  `library_writer` / `dataset_writer` (both `0.1`), or — in the oldest file
  `NdCeFeB_2-5.hdf5` — only `HT_class = 'HTroot'` with no writer version.
- Root groups: one group per measurement technique, one `sample` (and in
  `NdCeFeB_4_6.hdf5` an empty `samples` group with `HT_type='samples'`).

## Provenance: NOMAD uploads

- `NdCeFeB_2-5.hdf5` is the dataset that was used for the MaMMoS initial dataset
  published on NOMAD:
  - MaMMoS initial dataset: doi:10.17172/NOMAD/2025.06.29-4
  - `MaMMoS-NEEL_FeNdCeB_magnets`: doi:10.17172/NOMAD/2025.06.29-2
  (These are also listed under "Connected Datasets and their DOIS" in
  `Metadata for upload - comment_NEEL.md`.)
- `SmFeV_3378.hdf5` was also shared on NOMAD as a dataset:
  doi:10.17172/nomad.8wfe-0shj

## Measurement groups (first level)

Each measurement group carries attributes:

- `HT_type`: one of `edx`, `moke`, `xrd` (old), `xrd_wafer` (new ESRF),
  `profil`, `sem` (only in `NdCeFeB_3_13.hdf5`). Canonical enum in
  ht-data-parser `src/enums.py`: `HtType = xrd, edx, moke, sem, profil`.
- `<technique>_writer`: version string of the writer that produced the group
  (`'0.1 beta'` in `NdCeFeB_2-5.hdf5`; `'0.1 beta'`/`'0.2'`/`'0.3'` in
  `SmFeV_3378.hdf5`; `'0.9'` in the new files).
- `instrument`: instrument name (e.g. `Bruker Quantax Xflash-7`, `S-MOKE`,
  `bm02 esrf` / `bm02 - esrf`, `Bruker DektakXT`, `Zeiss Ultra Plus`).
- `profil` groups in the new files additionally have `source` (origin file).
- `NdCeFeB_2-5.hdf5` has two ESRF groups: `..._ESRF` (raw) and
  `..._ESRF_filtered`, the latter carrying a `note` attribute describing the
  filter (pyFAI `medfilt1d_ng`, percentile=(40,60)) and `hdf5_reader='keep'` —
  i.e. raw vs filtered variants can coexist as separate groups.

Group naming differs by file:

- `SmFeV_3378.hdf5`: plain names `edx`, `moke`, `xrd`, `dektak`.
- `NdCeFeB_2-5.hdf5`: dataset-prefixed names `NdCeFeB_1-1_EDX`,
  `NdCeFeB_2-5_10s650_MOKE`, `NdCeFeB_2-5_10s650_ESRF` (+ `_ESRF_filtered`)
  (prefix is NOT the file name).
- `NdCeFeB_4_6.hdf5`: dataset-prefixed names `NdCeFeB_3_9_EDX`, `NdCeFeB_4_6_MOKE`,
  `NdCeFeB_4_6_ESRF`, `NdCeFeB_4_8_PROFIL`.
- `NdCeFeB_3_13.hdf5`: mixed — `MOKE`, `ESRF`, `SEM` plain, EDX group named
  `NdCeFeB_3_13`, profil group `NdCeFeB_4_8_PROFIL`. Therefore group names must
  not be used to classify the technique; always use the `HT_type` attribute.

## Inner layout: two versions

### Old layout (writers <= 0.3, e.g. `NdCeFeB_2-5.hdf5`, `SmFeV_3378.hdf5`)

```
<measurement group>/
  (-10.0,-10.0)/            <- one group per wafer position, direct children
    instrument/
    measurement/
    results/
```

### New layout (writer 0.9, e.g. `NdCeFeB_4_6.hdf5`, `NdCeFeB_3_13.hdf5`)

```
<measurement group>/
  positions/
    (-10.0,-10.0)/          <- one group per wafer position
      instrument/
      measurement/
      results/
  experiment_info/          <- copy of the sample metadata
  additional_scans/
```

Detection rule used by the conversion script (and identical to the upstream
reader `Advanced_Data_Visualization/neel_data_vis/readers/read_hdf5.py`):
if a `positions` sub-group exists, coordinate groups live under
`<group>/positions/`; otherwise directly under `<group>`.

## Position groups

- Name format: `(<x>,<y>)` with one decimal, e.g. `(-10.0,-10.0)`, parsed with
  `\((-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)\)`.
- Attributes: `ignored` (bool, all `False` in current files) and `index`
  (scan index, int array or string).
- Children: `instrument/`, `measurement/`, `results/` (the `SEM` group in
  `NdCeFeB_3_13.hdf5` has only `instrument/` and `measurement/`, no `results/`).

### `instrument/`

- `x_pos`, `y_pos`: scalar float datasets with `units` attribute (`'mm'` in all
  current files; the conversion script converts `µm`/`um` to `mm` if needed).
- Technique-specific headers: `DetectorHeader`, `ESMAHeader`, `KnownHeader`,
  `Spectrum`, `SpectrumHardwareHeader`, `SpectrumHeader` (EDX); MOKE has
  `x_pos`, `y_pos`, `dataset_name`, `date`, `pulse_voltage`, `coil_to__sample_gap`,
  `scanning_dimension_x/y`, `number_of_points_x/y`, `all_points_measured?`, etc.

### `measurement/`

- EDX: `counts` (4096 int64, `cps`), `energy` (4096 float64, `keV`).
- MOKE: `shot_1..shot_4` and `shot_mean` groups with `pulse_*`,
  `magnetization_*`, `reflectivity_*`, `integrated_pulse_*` arrays plus `time`
  (`µs`).

### `results/`

- EDX: one group per element named `Element <symbol>` (e.g. `Element Nd`),
  containing `Element` (atomic number, float), `AtomPercent` (units `at.%`),
  `MassPercent` (units `mass.%`), `NetIntens`, `Sigma`, `DetectionLimit`,
  `XLine`, `Types`, ... plus a `TRTResult` group. Elements without
  `AtomPercent` (e.g. C, O, Si, Ta substrate signals) are skipped by the
  conversion script.
- MOKE: groups `coercivity_m0`, `coercivity_dmdh`, `intercept_field`,
  `remanent_kerr_signal`, each with `mean`/`negative`/`positive` datasets
  (units `T`); scalars `max_kerr_signal`, `reflectivity` (units `V`); and a
  `parameters` group. The conversion script uses `coercivity_m0/mean` and
  converts T to A/m via mu0.
- ESRF/XRD (`xrd_wafer`): phase results (lattice params, phase fractions) —
  not yet processed by the conversion script.
- Profil: thickness results — not yet processed.

## `sample` / `samples` / `experiment_info`

- `NdCeFeB_2-5.hdf5`: `sample` has attr `HT_class='sample'` with scalar datasets
  `fabrication_date` ('16.04.2025'), `operator` ('William'), `sample_name`
  ('NdCeFeB_2-5') plus layer groups `layer 0`, `layer 1`, `layer 2` (each with
  `Element`, `Thickness`).
- `SmFeV_3378.hdf5`: `sample` has attr `HT_class='sample'` and scalar datasets
  `fabrication_date`, `file_type='Library'`, `operator`, `sample_name`.
- `NdCeFeB_4_6.hdf5`: `sample` (no attrs) contains groups `NdCeFeB_4_6` and
  `layer_1` with deposition parameters (`element='NdCeFeB'`, `power=50 W`,
  `distance=50 mm`, `thickness=800 nm`, `temperature`, `time`, `angle`,
  `type='active'`); plus an empty `samples` group (`HT_type='samples'`).
- `NdCeFeB_3_13.hdf5`: `sample` has attr `HT_type='sample'` and mixes the old
  scalar fields (`fabrication_date`, `operator`, `sample_name`) with the new
  `NdCeFeB_3_13` / `layer_1` groups (`thickness=400 nm`, `time=1800 s`).
- New-layout measurement groups each contain an `experiment_info/sample` copy
  of this metadata.

## Compatibility with upstream MaMMoS repos

Field names used by the conversion script match the upstream tools:

| Concept | upstream source | conversion script |
|---|---|---|
| Technique detection | `read_hdf5.py`: `group.attrs["HT_type"] == data_type.lower()`; `enums.py HtType` | `get_ht_type()` / `find_groups_by_ht_type()` |
| Layout detection | `read_hdf5.py`: `if "positions" in group.keys()` | `get_position_group_paths()` |
| Position naming | `read_hdf5.py`: `f"({round(x,1)},{round(y,1)})"` | `extract_coordinates_from_name()` regex |
| Position coords | `read_hdf5.py`: `instrument["x_pos"][()]`, units `mm` | `get_instrument_positions()` |
| EDX composition | `read_edx.py`: `Element X` groups, `AtomPercent`/`MassPercent` + `units` | `get_element_data()` |
| MOKE coercivity | `read_moke.py`: results groups with `mean` + `units`; writer uses `coercivity_m0` (T) | `get_moke_data()` (`coercivity_m0/mean`, T→A/m) |

Upstream authorship (all Institut Néel): `Advanced_Data_Visualization`
(williamrigaut, latest commit "added compatibility with new hdf5 format"),
`DaHU` and `ht-data-parser` (Escouflenfer; ht-data-parser exports HDF5
"NeXus-inspired" with the MaMMoS ontology and defines the `HtType` enum).

## Conversion script support (branch `neel-ht-dataset-vs-library`)

- `find_groups_by_ht_type()` detects measurement groups via the `HT_type`
  attribute (case-insensitive), with a name-substring fallback for files
  without the attribute.
- `get_position_group_paths()` supports both layouts (`positions/` sub-group
  or direct children).
- Each sample is tagged `data_type='edx'/'moke'` during processing and
  `merge_edx_and_moke_data()` uses this tag (name match kept as fallback) —
  required because e.g. the EDX group of `NdCeFeB_3_13.hdf5` is named
  `NdCeFeB_3_13` and contains no `edx` substring.
- Running without arguments processes all `*.hdf5` in `datasets/`; CLI:
  `-f/--file`, `-m/--mass-fraction`, `-z/--create-zip`, `-q/--quiet`, `-h/--help`.

Verification (all four datasets): 4 × 289 EDX positions + 4 × 225 MOKE
positions, 2056/2056 coordinate matches and 3395/3395 element verifications
(100%), 1156 generated `*.archive.yaml` files in `generated_schemas/`
(225 EDX+MOKE merged + 64 EDX-only per dataset), no UNKNOWN or empty
compositions.

## Known quirks / open points

- Root attribute changed over time: `NdCeFeB_2-5.hdf5` uses `HT_class='HTroot'`
  (no root `HT_type`/writer), later files use `HT_type` (`library`/`dataset`)
  plus `library_writer`/`dataset_writer`.
- `xrd` was renamed to `xrd_wafer` in the new files; the conversion script
  currently processes only EDX and MOKE.
- `NdCeFeB_2-5.hdf5` ships raw and filtered ESRF data as separate groups
  (`..._ESRF` and `..._ESRF_filtered` with a `note` attribute); which variant
  to describe/upload is an open point (raw vs filtered).
- New `sem` HT_type in `NdCeFeB_3_13.hdf5` (no `results/` sub-group).
- `NdCeFeB_4_6.hdf5` contains an empty `samples` group; `sample` metadata moved
  into `experiment_info` inside each measurement group.
- Companion CSVs (`NdCeFeB_4_6.csv`, `datasets/NdCeFeB_3_13/.../NdCeFeB_3_13.csv`)
  hold wafer maps incl. ESRF phase data, coercivity, and process parameters.
