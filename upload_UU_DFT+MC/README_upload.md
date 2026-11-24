# UU (Uppsala) DFT+MC dataset upload — MaMMoS

This upload contains one atomistic spin-dynamics dataset (**RSPt DFT + UppASD
Monte Carlo**) computed at Uppsala University, Department of Physics and
Astronomy; Materials Theory, within the EU Horizon Europe project **MaMMoS**
(grant agreement No 101135546, https://mammos-project.github.io/).

The raw and derived data come from
[github.com/MaMMoS-project/uppsala-data](https://github.com/MaMMoS-project/uppsala-data)
(one folder per composition); the derived files (`intrinsic_properties.yaml`,
`UppASD/MC_*/thermal.csv`) are generated and validated with
[mammos-parser](https://github.com/MaMMoS-project/mammos-parser).

## Where to find the dataset-specific information in NOMAD

This README is identical for all UU uploads; the concrete dataset name, chemical
composition and values of **this** entry are in:

- the **entry page → Overview** tab: Material section (chemical formula and
  elements), Comment;
- the **entry page → Data** tab: the data section with the quantities
  (`Js`, `Js_300`, `Ku`, `Tc`, `Ms` — linked to the uploaded `.h5`, and the
  repeating `thermal_properties` subsection with the temperature-dependent
  curves `T, Ms, Js, E, Cv, chi, U_L`), plus `mammos_entities` with the EMMO
  ontology labels/IRIs;
- the **raw files** of this upload (see below), in particular
  `<dataset>.archive.yaml` and `<dataset>.h5`.

## Files in this upload

- `<dataset>.h5` — the data as a mammos-entity `EntityCollection`
  (`intrinsic_properties`, `thermal_MC_<n>`, `raw_data`); read with
  `import mammos_entity as me; me.from_hdf5('<dataset>.h5')`
- `<dataset>.archive.yaml` — the NOMAD entry for this dataset
- `README.md` — this file (general, identical for all UU uploads)
- `cif_structure_<dataset>.cif` and `datasets/<dataset>/` — crystal structure
  and the full raw data tree
- `Evaluation_dataset_UU.ipynb` — notebook showing how to read and analyse the
  data with `mammos_entity` (works with the local files of the unpacked upload)
- `packages/` — the conversion/validation code that built this upload
