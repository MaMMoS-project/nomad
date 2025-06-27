# Code used to create a upload on Nomad for data from UU

Before starting, **create a folder named `datasets`** in your project directory. This folder will store all simulation datasets to be processed.

Into the folder called 'datasets' one should put the simulation dataset to be processed. Each dataset subfolder should have the following structure:

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
# Process a specific dataset (e.g., Fe12Co4N2)
python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir Fe12Co4N2

# Process with custom datasets directory
python run_convert_UU_data_to_h5_and_creat_yaml.py --source_dir Fe12Co4N2 --datasets_dir ./datasets
```

### Processing Multiple Datasets

To process multiple datasets at once, use the `batch_process_datasets.py` script:

```bash
# Process all datasets listed in the batch script
python batch_process_datasets.py
```

Edit the `datasets_to_process` list in `batch_process_datasets.py` to specify which datasets to process.

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

