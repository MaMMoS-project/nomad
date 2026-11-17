# Input datasets

Place the original NEEL HDF5 datasets (`.hdf5`) in this folder, e.g.
`NdCeFeB_2-5.hdf5`, `SmFeV_3378.hdf5`, `NdCeFeB_4_6.hdf5`, `NdCeFeB_3_13.hdf5`.

Then run from the parent folder:

```bash
python run_convert_NEEL_data_from_hdf5_and_create_schemas.py        # processes all *.hdf5 here
python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -h     # help
```

Generated NOMAD YAML schemas are written to `../generated_schemas/`.

Note: ALL `.hdf5` files in this folder are processed with ONE single template,
`../NEEL_template.archive.yaml` — there is no per-dataset template selection yet
(known limitation; the template's `description` text fits only one dataset).
See the README in the parent folder ("Known Limitations").

This README is only a placeholder so the folder exists at checkout; the `.hdf5`
files themselves are gitignored (they are several GB each) and must be obtained
from NEEL separately.
