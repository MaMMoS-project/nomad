# NEEL Data Processing for NOMAD Upload

This repository contains code to process HDF5 datasets from NEEL and generate NOMAD-compatible YAML schemas for data publication.

## Overview

The script extracts metadata and measurement data from NEEL's HDF5 files, including:
- Sample positions on the wafer (x, y coordinates)
- Chemical/elemental composition from EDX analysis
- Coercivity measurements from MOKE analysis
- Atomic fractions with recalculated values including Boron

## Directory Structure

The script expects and creates the following directory structure:

```
project_root/
├── run_convert_NEEL_data_from_hdf5_and_create_schemas.py  # Main processing script
├── NEEL_template.archive.yaml                             # YAML template file
├── datasets/                                              # INPUT: Original HDF5 files
│   ├── NdCeFeB_2-5.hdf5                                   # Example dataset from NEEL
│   ├── NdFeB.hdf5                                         # Example dataset from NEEL
│   └── ...                                                # Other HDF5 files
├── generated_schemas/                                      # OUTPUT: Generated YAML schemas
│   ├── NdCeFeB_2-5_EDX_MOKE_xpos=0_0_ypos=0_0.archive.yaml
│   ├── NdCeFeB_2-5_EDX_MOKE_xpos=0_0_ypos=10_0.archive.yaml
│   └── ...                                                # One YAML per sample position
└── uploads/                                               # OUTPUT: Optional zip files
    ├── NdCeFeB_2-5_20250627_143022.zip                   # Contains HDF5 + YAML files
    └── ...                                                # Timestamped zip files
```

## Usage

### Basic Workflow

1. **Place original HDF5 datasets** from NEEL into the `datasets/` subfolder
2. **Run the processing script**:
   ```bash
   python run_convert_NEEL_data_from_hdf5_and_create_schemas.py
   ```
3. **Find generated YAML schemas** in the `generated_schemas/` folder
4. **Optionally create zip files** for NOMAD upload (when `create_zip=True`)

### Script Options

The script provides several functions with different options:

```python
# Process single file (default behavior)
process_single_file()                                        # Default file, no zip
process_single_file("filename.hdf5")                         # Specific file
process_single_file("filename.hdf5", create_zip=True)        # With zip file

# Process all files in datasets/
process_all_files()                                          # All files, no zip
process_all_files(create_zip=True)                           # All files with zip
process_all_files(include_mass_fraction=True)                # Include mass fractions

# Direct function calls
main(single_file="filename.hdf5")                            # Specific file
main(single_file="", create_zip=True)                        # All files with zip
```

## Output Files

### Generated YAML Schemas
- Location: `generated_schemas/` folder
- Format: One YAML file per sample position
- Naming: `{dataset}_{datatype}_xpos={x}_ypos={y}.archive.yaml`
- Content: NOMAD-compatible metadata and measurement data

### Optional Zip Files
- Location: `uploads/` folder
- Contains: Original HDF5 files + generated YAML schemas
- Created: Only when `create_zip=True` parameter is used
- Naming: Timestamped for uniqueness

## Important Notes

### Requirements
- **Original datasets from NEEL MUST be placed in the `datasets/` subfolder**
- The `NEEL_template.archive.yaml` template file must be present
- Python with h5py library is required

### File Size Considerations
- **ZIP files can become very large** (several GB) if original HDF5 datasets are large
- Zip file creation is **disabled by default** to prevent accidentally creating huge files
- Only enable zip creation (`create_zip=True`) when specifically needed for NOMAD upload

### Data Processing Features
- Automatic coordinate extraction and verification
- Chemical element identification and atomic fraction calculation
- Recalculation of atomic fractions including Boron for NdCeFeB compounds
- MOKE coercivity data integration when available
- Mass fraction support (optional, commented out by default)

## File Structure Details

Each generated YAML file contains:
- Sample identification and coordinates
- Elemental composition with atomic fractions
- Optional coercivity measurements from MOKE
- Metadata for NOMAD compatibility
- Recalculated atomic fractions including Boron (for applicable compounds)