import os
import re
import sys
import zipfile
import argparse
import pathlib
import numpy as np
from datetime import datetime

# =============================================================================
# Periodic Table Data - Single Source of Truth
# =============================================================================
# Symbol → Atomic Number mapping (first 118 elements)
PERIODIC_TABLE = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Nh": 113,
    "Fl": 114,
    "Mc": 115,
    "Lv": 116,
    "Ts": 117,
    "Og": 118,
}

# Derived: Atomic Number → Symbol (reverse mapping)
ATOMIC_NUMBER_TO_SYMBOL = {v: k for k, v in PERIODIC_TABLE.items()}

# Derived: Set of all known element symbols (immutable)
KNOWN_ELEMENTS = frozenset(PERIODIC_TABLE.keys())

# Import neel_data_vis for high-level xarray-based HDF5 access (REQUIRED)
# Try multiple locations: installed package first, then local subfolder
NEEL_DATA_VIS_AVAILABLE = False
read_hdf5 = None

# First, try to import from installed package
try:
    import neel_data_vis.readers.read_hdf5 as read_hdf5

    NEEL_DATA_VIS_AVAILABLE = True
    print(
        "neel_data_vis library found (installed package), using xarray-based data access"
    )

except ImportError:
    # If not installed, try to load from local subfolder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_neel_data_vis_path = os.path.join(
        script_dir, "Advanced_Data_Visualization_mod", "Advanced_Data_Visualization"
    )

    if os.path.exists(local_neel_data_vis_path):
        # Add the local path to sys.path temporarily
        if local_neel_data_vis_path not in sys.path:
            sys.path.insert(0, local_neel_data_vis_path)

        try:
            import neel_data_vis.readers.read_hdf5 as read_hdf5

            NEEL_DATA_VIS_AVAILABLE = True
            print(
                "neel_data_vis library found (local subfolder), using xarray-based data access"
            )
            print(f"  Location: {local_neel_data_vis_path}")
        except ImportError as e:
            print(f"ERROR: Failed to import neel_data_vis from local subfolder: {e}")
            print("This script requires neel_data_vis library for data extraction.")
            sys.exit(1)
    else:
        print(
            "ERROR: neel_data_vis not available (neither installed nor in local subfolder)"
        )
        print(f"  Checked local path: {local_neel_data_vis_path}")
        print("This script requires neel_data_vis library for data extraction.")
        sys.exit(1)

"""
NEEL Data Processing Script for NOMAD Schema Generation

DIRECTORY STRUCTURE:
This script expects and creates the following directory structure:

    project_root/
    ├── run_convert_NEEL_data_from_hdf5_and_create_schemas.py  (this script)
    ├── NEEL_template.archive.yaml                             (template file)
    ├── datasets/                                              (INPUT: place original .hdf5 files here)
    │   ├── SmFeV_3378.hdf5                                    (example dataset from NEEL)
    │   ├── NdCeFeB_2-5.hdf5                                   (example dataset from NEEL)
    │   └── ...                                                (other .hdf5 files)
    ├── generated_schemas/                                      (OUTPUT: generated YAML schemas)
    │   ├── SmFeV_3378_EDX_MOKE_xpos=0_0_ypos=0_0.archive.yaml
    │   ├── SmFeV_3378_EDX_MOKE_xpos=0_0_ypos=10_0.archive.yaml
    │   └── ...                                                (one YAML per sample position)
    └── uploads/                                               (OUTPUT: optional zip files for NOMAD upload)
        ├── SmFeV_3378_20260114_143022.zip                    (contains .hdf5 + .yaml files)
        └── ...                                                (timestamped zip files)

WORKFLOW:
1. Place original HDF5 datasets from NEEL into the 'datasets/' subfolder
2. Run this script: python run_convert_NEEL_data_from_hdf5_and_create_schemas_xarray_ADV.py
3. Generated YAML schemas for NOMAD entries will be created in 'generated_schemas/'
4. Optionally (if create_zip=True), zip files containing all data will be created in 'uploads/'

IMPORTANT NOTES:
- Original datasets from NEEL MUST be placed in a subfolder called 'datasets/'
- The generated YAML schemas can be found under 'generated_schemas/' after running the script
- ZIP files for NOMAD upload are created in 'uploads/' only when create_zip=True
- ZIP files can become very large, especially if the original HDF5 datasets are large
- Each YAML file corresponds to one sample position with EDX/MOKE measurement data

FILE SIZE WARNING:
The zip files in 'uploads/' can grow quite large (several GB) if the original HDF5 
datasets are large, as they contain both the original data and generated schemas.
"""

# Note: Chemical formula calculation is implemented but not exported to YAML
#       (atomic fractions are recalculated including B and exported instead)
#
# DEPENDENCIES:
#   - neel_data_vis: Required for xarray-based HDF5 data extraction
#   - numpy: For numerical operations
#   - Standard library: os, re, sys, zipfile, argparse, pathlib, datetime


def compute_stoichiometric_coefficients_from_fractions(nd_fraction, ce_fraction):
    """
    Compute stoichiometric coefficients for Nd_a Ce_b Fe_c B compound from atomic fractions.

    Based on the assumption that Fe only is in the main phase (Nd,Ce)2Fe14B: c = 14
    nd_fraction = a/(a+b+14)
    ce_fraction = b/(a+b+14)
    fe_fraction = 14/(a+b+14)

    Therefore:
    a = (14*nd_fraction)/fe_fraction
    b = (14*ce_fraction)/fe_fraction
    where fe_fraction = 1-nd_fraction-ce_fraction

    Parameters:
    nd_fraction (float): Atomic fraction of Nd from EDX analysis.
    ce_fraction (float): Atomic fraction of Ce from EDX analysis.

    Returns:
    tuple: (a, b, c) - stoichiometric coefficients for Nd_a Ce_b Fe_c B compound
    """
    # NdCeFeB Composition conversion from EDX (without B) to 2:14:1 (including B)
    fe_fraction = 1 - nd_fraction - ce_fraction
    if fe_fraction <= 0:
        # Invalid input - Fe fraction must be positive
        return None, None, None

    a = (14 * nd_fraction) / fe_fraction
    b = (14 * ce_fraction) / fe_fraction
    c = 14  # Assuming Fe is only in the main phase

    return a, b, c


def calculate_atomic_fractions_with_boron(a, b, c=14):
    """
    Calculate atomic fractions including Boron from stoichiometric coefficients.

    For compound Nd_a Ce_b Fe_c B, the atomic fractions are:
    - Nd: a / (a + b + c + 1)
    - Ce: b / (a + b + c + 1)
    - Fe: c / (a + b + c + 1)
    - B:  1 / (a + b + c + 1)

    Parameters:
    a (float): Stoichiometric coefficient for Nd
    b (float): Stoichiometric coefficient for Ce
    c (float): Stoichiometric coefficient for Fe (default 14)

    Returns:
    dict: Dictionary with atomic fractions for each element
    """
    total = a + b + c + 1  # +1 for the single B atom

    return {
        "Nd": a / total,
        "Ce": b / total,
        "Fe": c / total,
        "B": 1 / total,
    }


def update_elements_with_recalculated_fractions(elements, recalculated_fractions):
    """
    Update the elements dictionary with recalculated atomic fractions.

    Parameters:
    elements (dict): Original elements dictionary from EDX data
    recalculated_fractions (dict): Recalculated atomic fractions including B

    Returns:
    dict: Updated elements dictionary with recalculated fractions
    """
    updated_elements = elements.copy()

    # Update existing elements with recalculated fractions
    for element_symbol, new_fraction in recalculated_fractions.items():
        if element_symbol in updated_elements:
            # Update existing element
            updated_elements[element_symbol]["atom_percent"] = new_fraction * 100.0
            updated_elements[element_symbol]["recalculated"] = True
        elif element_symbol == "B":
            # Add Boron if it doesn't exist
            updated_elements["B"] = {
                "symbol_match": True,
                "atom_percent": new_fraction * 100.0,
                "mass_percent": None,  # We don't calculate mass percent for B
                "recalculated": True,
                "added_from_formula": True,
            }

    return updated_elements


def check_for_nd_ce_fe_only(elements):
    """
    Check if the elements contain only Nd, Ce, and Fe (and optionally B).

    Args:
        elements (dict): Dictionary of element data

    Returns:
        tuple: (is_nd_ce_fe_only, nd_fraction, ce_fraction, fe_fraction)
    """
    if not elements:
        return False, None, None, None

    # Get element symbols
    element_symbols = set(elements.keys())

    # Check if we have only Nd, Ce, Fe (B is assumed to be present but not measured)
    expected_elements = {"Nd", "Ce", "Fe"}

    # Must have exactly these elements or a subset including at least Fe
    if not element_symbols.issubset(expected_elements):
        return False, None, None, None

    # Must have Fe, and at least one of Nd or Ce
    if "Fe" not in element_symbols:
        return False, None, None, None

    if not ("Nd" in element_symbols or "Ce" in element_symbols):
        return False, None, None, None

    # Calculate atomic fractions (normalized)
    total_atom_percent = sum(
        elem_data["atom_percent"] for elem_data in elements.values()
    )

    if total_atom_percent <= 0:
        return False, None, None, None

    # Get normalized atomic fractions
    nd_fraction = elements.get("Nd", {}).get("atom_percent", 0) / total_atom_percent
    ce_fraction = elements.get("Ce", {}).get("atom_percent", 0) / total_atom_percent
    fe_fraction = elements.get("Fe", {}).get("atom_percent", 0) / total_atom_percent

    return True, nd_fraction, ce_fraction, fe_fraction


def get_element_symbol_from_atomic_number(atomic_number):
    """
    Get element symbol from atomic number.

    Args:
        atomic_number (int): Atomic number

    Returns:
        str: Element symbol or 'Unknown' if not found
    """
    return ATOMIC_NUMBER_TO_SYMBOL.get(atomic_number, "Unknown")


# =============================================================================
# Data extraction functions using neel_data_vis (xarray approach)
# =============================================================================


def load_dataset_with_neel_data_vis(file_path, verbose=True):
    """
    Load HDF5 dataset using neel_data_vis library which returns an xarray Dataset.

    Args:
        file_path (str or pathlib.Path): Path to the HDF5 file
        verbose (bool): If True, print detailed output

    Returns:
        xarray.Dataset or None: The loaded dataset, or None if loading failed
    """
    if not NEEL_DATA_VIS_AVAILABLE:
        if verbose:
            print("neel_data_vis not available, cannot load dataset")
        return None

    try:
        path = pathlib.Path(file_path)
        # Use exclude_wafer_edges=True to exclude edge positions
        data = read_hdf5.get_full_dataset(path, exclude_wafer_edges=False)
        if verbose:
            print(
                f"Successfully loaded dataset with {len(data.data_vars)} data variables"
            )
            print(f"Available data variables: {list(data.data_vars)}")
            print(f"Coordinate dimensions: {list(data.coords)}")
        return data
    except Exception as e:
        if verbose:
            print(f"Error loading dataset with neel_data_vis: {e}")
        return None


def extract_all_data_with_xarray(file_path, verbose=True):
    """
    Extract all coordinate data (EDX and MOKE) from an HDF5 file using neel_data_vis xarray approach.

    This function uses the neel_data_vis xarray interface for clean coordinate-based data access.
    It iterates over the 2D grid of positions and extracts data for each valid position.

    Args:
        file_path (str or pathlib.Path): Path to the HDF5 file
        verbose (bool): If True, print detailed output

    Returns:
        dict: Dictionary with coordinate data organized by sample position
    """
    all_coordinate_data = {}

    # Load dataset using neel_data_vis
    data = load_dataset_with_neel_data_vis(file_path, verbose)

    if data is None:
        if verbose:
            print("Failed to load dataset with neel_data_vis")
        return None  # No data available

    # Check what data variables are available
    data_vars = list(data.data_vars)
    if verbose:
        print("\nAvailable data variables in dataset:")
        for var in data_vars:
            print(f"  - {var}")

    # Check for coercivity data (MOKE)
    # Look for various possible naming patterns for coercivity
    coercivity_var_name = None
    coercivity_patterns = [
        "coercivity_m0",
        "coercivity",
        "Coercivity",
        "hc",
        "Hc",
        "coercive_field",
        "Coercivity M0",
        "coercivity m0",
    ]

    for pattern in coercivity_patterns:
        if pattern in data:
            coercivity_var_name = pattern
            if verbose:
                print(f"Found coercivity variable (exact match): '{pattern}'")
            break

    # Also check for partial matches in data_vars
    if coercivity_var_name is None:
        for var in data_vars:
            var_lower = var.lower()
            # Check for coercivity patterns or m0 patterns (MOKE specific)
            if (
                "coercivity" in var_lower
                or "coercive" in var_lower
                or "_m0" in var_lower
                or " m0" in var_lower
            ):
                coercivity_var_name = var
                if verbose:
                    print(f"Found coercivity variable (partial match): '{var}'")
                break

    coercivity_available = coercivity_var_name is not None
    coercivity = None
    coercivity_unit = "T"  # Default unit

    if coercivity_available:
        coercivity = data[coercivity_var_name]
        coercivity_unit = coercivity.attrs.get("units", "T")
        if verbose:
            print(
                f"\nCoercivity data available as '{coercivity_var_name}' with unit: {coercivity_unit}"
            )
    else:
        if verbose:
            print("\nNo coercivity data variable found in xarray dataset")

    # Check for EDX element data - look for composition/atomic fraction variables
    # Common naming patterns:
    #   - "Fe Composition", "Sm Composition" (element name + Composition)
    #   - element_<symbol>_atom_percent, <symbol>_atomic_fraction
    edx_element_vars = {}

    for var in data_vars:
        var_lower = var.lower()

        # Pattern 1: "<Element> Composition" (e.g., "Fe Composition", "Sm Composition")
        if "composition" in var_lower:
            # Extract the element symbol (first word before "Composition")
            parts = var.split()
            if len(parts) >= 2 and parts[0] in KNOWN_ELEMENTS:
                edx_element_vars[parts[0]] = var
                continue

        # Pattern 2: "atom" + "percent/fraction" patterns
        if "atom" in var_lower and ("percent" in var_lower or "fraction" in var_lower):
            # Extract element symbol from variable name
            parts = var.split("_")
            for part in parts:
                if len(part) <= 3 and part.isalpha() and part[0].isupper():
                    if part in KNOWN_ELEMENTS:
                        edx_element_vars[part] = var
                        break

    if edx_element_vars and verbose:
        print("\nDetected EDX element variables:")
        for elem, var in edx_element_vars.items():
            print(f"  - {elem}: {var}")
    elif verbose:
        print("\nNo EDX element variables detected in dataset")

    # Collect ALL unique positions from ANY data variable that has valid data
    # This ensures we don't miss positions that have EDX but no dektak/MOKE data
    all_valid_positions = set()  # Set of (n, m, x, y) tuples

    # Get grid shape from any available variable
    first_var = data[data_vars[0]]
    shape = first_var.shape
    if verbose:
        print(f"\nGrid shape: {shape}")

    # Scan all data variables to find ALL positions with valid data
    vars_to_scan = []

    # Add EDX element variables (these typically have the most positions - 289)
    for var_name in edx_element_vars.values():
        vars_to_scan.append(var_name)

    # Add coercivity variable if available
    if coercivity_var_name:
        vars_to_scan.append(coercivity_var_name)

    # Add measured_height if available
    if "measured_height" in data:
        vars_to_scan.append("measured_height")

    if verbose:
        print(f"\nScanning {len(vars_to_scan)} data variables for valid positions...")

    for var_name in vars_to_scan:
        var_data = data[var_name]
        for n in range(shape[0]):
            for m in range(shape[1]):
                value = var_data[n, m].values
                if not np.isnan(value):
                    x_pos = float(var_data[n, m].x.values)
                    y_pos = float(var_data[n, m].y.values)
                    all_valid_positions.add((n, m, x_pos, y_pos))

    if verbose:
        print(f"Found {len(all_valid_positions)} unique positions with valid data")

    # Now iterate over ALL valid positions and extract available data
    edx_group_data = {}
    moke_group_data = {}

    for n, m, x_pos, y_pos in all_valid_positions:
        # Create sample key based on coordinates
        sample_key = f"({int(x_pos) if x_pos == int(x_pos) else x_pos},{int(y_pos) if y_pos == int(y_pos) else y_pos})"

        # Build sample data structure
        sample_data = {
            "x_pos_instrument": x_pos,
            "y_pos_instrument": y_pos,
            "x_pos_unit": "mm",  # neel_data_vis typically uses mm
            "y_pos_unit": "mm",
            "x_match": True,  # Using xarray, coordinates are consistent
            "y_match": True,
        }

        has_any_data = False

        # Extract MOKE coercivity if available from xarray
        if coercivity_available:
            hc_value = float(coercivity[n, m].values)
            if not np.isnan(hc_value):
                sample_data["moke_data"] = {
                    "coercivity_mean": hc_value,
                    "coercivity_unit": coercivity_unit,
                    "data_available": True,
                }
                # Add to MOKE group data
                moke_sample = sample_data.copy()
                moke_sample["x_pos_MOKE"] = x_pos
                moke_sample["y_pos_MOKE"] = y_pos
                moke_group_data[sample_key] = moke_sample
                has_any_data = True

        # Extract EDX element data if available
        elements = {}
        for elem_symbol, var_name in edx_element_vars.items():
            elem_data = data[var_name]
            elem_value = float(elem_data[n, m].values)
            if not np.isnan(elem_value):
                # Get atomic number for element
                atomic_number = get_atomic_number_from_symbol(elem_symbol)
                elements[elem_symbol] = {
                    "element_symbol": elem_symbol,
                    "atomic_number": atomic_number,
                    "atom_percent": elem_value,
                    "mass_percent": None,  # May not be available in xarray format
                    "symbol_match": True,  # Assumed correct in neel_data_vis format
                    "expected_symbol": elem_symbol,
                    "group_name": f"Element {elem_symbol}",
                }

        if elements:
            sample_data["elements"] = elements
            sample_data["x_pos_EDX"] = x_pos
            sample_data["y_pos_EDX"] = y_pos
            edx_group_data[sample_key] = sample_data
            has_any_data = True
        elif has_any_data:
            # Has MOKE but no EDX - already stored in moke_group_data
            pass

    # Organize data into the expected format
    filename = os.path.basename(file_path)
    all_coordinate_data[filename] = {}

    if edx_group_data:
        all_coordinate_data[filename]["edx_xarray"] = edx_group_data
        if verbose:
            print(f"\nExtracted EDX data for {len(edx_group_data)} positions")

    if moke_group_data:
        all_coordinate_data[filename]["moke_xarray"] = moke_group_data
        if verbose:
            print(f"Extracted MOKE data for {len(moke_group_data)} positions")

    return all_coordinate_data


def get_atomic_number_from_symbol(symbol):
    """
    Get atomic number from element symbol.

    Args:
        symbol (str): Element symbol (e.g., 'Fe', 'Nd')

    Returns:
        int: Atomic number or 0 if not found
    """
    return PERIODIC_TABLE.get(symbol, 0)


# =============================================================================
# Template processing and YAML generation functions
# =============================================================================


def handle_mass_fraction_template(template_content, include_mass_fraction):
    """
    Handle mass_fraction related modifications in the template content.

    Args:
        template_content (str): The template content
        include_mass_fraction (bool): Whether to include mass_fraction in the output

    Returns:
        str: Modified template content
    """
    if include_mass_fraction:
        # Uncomment mass_fraction lines in the template
        # Uncomment in quantities section
        template_content = re.sub(
            r"(\s*)# mass_fraction:\s*\n(\s*)#\s*type:\s*np\.float64",
            r"\1mass_fraction:\n\2  type: np.float64",
            template_content,
            flags=re.MULTILINE,
        )
        # Uncomment in data section
        template_content = re.sub(
            r"(\s*)# mass_fraction:\s*\$\$mass_fraction\$\$",
            r"\1mass_fraction: $$mass_fraction$$",
            template_content,
        )
    else:
        # Ensure mass_fraction lines are commented out
        # Comment in quantities section if not already commented
        template_content = re.sub(
            r"(\s*)mass_fraction:\s*\n(\s*)type:\s*np\.float64",
            r"\1# mass_fraction:\n\2#   type: np.float64",
            template_content,
            flags=re.MULTILINE,
        )
        # Comment in data section if not already commented
        template_content = re.sub(
            r"(\s*)mass_fraction:\s*\$\$mass_fraction\$\$",
            r"\1# mass_fraction: $$mass_fraction$$",
            template_content,
        )

    return template_content


def create_yaml_from_template(
    template_path, output_path, sample_data, sample_key, include_mass_fraction=False
):
    """
    Create a YAML file from template by replacing placeholders with actual data.

    Args:
        template_path (str): Path to the template YAML file
        output_path (str): Path for the output YAML file
        sample_data (dict): Sample data containing coordinates and elements
        sample_key (str): Sample key (coordinate group name)
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
    """
    try:
        # Read template file
        with open(template_path, "r", encoding="utf-8") as f:
            template_content = f.read()

        # Handle mass_fraction template modifications
        template_content = handle_mass_fraction_template(
            template_content, include_mass_fraction
        )

        # Extract coordinate data
        x_pos = sample_data.get("x_pos_instrument") or sample_data.get("x_pos_EDX", 0.0)
        y_pos = sample_data.get("y_pos_instrument") or sample_data.get("y_pos_EDX", 0.0)

        # Convert units from μm to mm if needed
        x_unit = sample_data.get("x_pos_unit", "")
        y_unit = sample_data.get("y_pos_unit", "")

        if x_unit == "μm" or x_unit == "um":
            x_pos = x_pos / 1000.0  # Convert μm to mm
        if y_unit == "μm" or y_unit == "um":
            y_pos = y_pos / 1000.0  # Convert μm to mm

        # Replace coordinate placeholders
        template_content = template_content.replace("$$xpos$$", str(x_pos))
        template_content = template_content.replace("$$ypos$$", str(y_pos))

        # Handle MOKE coercivity data
        moke_data = sample_data.get("moke_data", {})
        if moke_data.get("data_available", False) and "coercivity_mean" in moke_data:
            coercivity_value = moke_data["coercivity_mean"]
            print(f"Using extracted MOKE coercivity: {coercivity_value}")
            if moke_data.get("coercivity_unit"):
                print(
                    f"Coercivity unit original in hdf5: {moke_data['coercivity_unit']}"
                )
            if moke_data["coercivity_unit"] == "T":
                print("Coercivity unit is Tesla, converting to A/m")
                mu0 = 4 * np.pi * 1e-7  # Vacuum permeability in T*m/A
                coercivity_value_in_A_per_m = (
                    float(coercivity_value) / mu0
                )  # Convert T to A/m

            # Replace coercivity placeholder with actual value
            template_content = template_content.replace(
                "$$coercivity$$", str(coercivity_value_in_A_per_m)
            )
        else:
            print("No MOKE data available, removing coercivity field from template")
            # Remove the coercivity field entirely if no MOKE data
            template_content = re.sub(
                r"\s*CoercivityHcExternal:\s*\$\$coercivity\$\$\s*\n",
                "\n",
                template_content,
            )
            # Also remove from schema definition if present
            template_content = re.sub(
                r"\s*CoercivityHcExternal:\s*\n\s*type:\s*np\.float64\s*\n\s*unit:\s*A/m\s*\n\s*description:\s*\'Coercivity\s*from\s*MOKE\s*measurements\s*\(optional\)\'\s*\n",
                "\n",
                template_content,
            )

        # Handle elemental composition
        elements = sample_data.get("elements", {})

        # Check if we can recalculate atomic fractions with Boron
        recalculated_elements = elements
        if elements:
            is_nd_ce_fe_only, nd_fraction, ce_fraction, fe_fraction = (
                check_for_nd_ce_fe_only(elements)
            )
            if is_nd_ce_fe_only and nd_fraction is not None and ce_fraction is not None:
                # Compute stoichiometric coefficients
                a, b, c = compute_stoichiometric_coefficients_from_fractions(
                    nd_fraction, ce_fraction
                )
                if a is not None and b is not None and c is not None:
                    # Calculate recalculated atomic fractions including Boron
                    recalculated_fractions = calculate_atomic_fractions_with_boron(
                        a, b, c
                    )
                    # Update elements with recalculated fractions
                    recalculated_elements = update_elements_with_recalculated_fractions(
                        elements, recalculated_fractions
                    )
                    print(
                        f"Using recalculated atomic fractions including B for {sample_key}"
                    )

        if recalculated_elements:
            # Create elemental composition entries
            elemental_entries = []
            for element_symbol, element_info in recalculated_elements.items():
                if (
                    element_info.get("symbol_match", False)
                    and element_info.get("atom_percent") is not None
                ):  # Only include verified elements with AtomPercent
                    atom_percent = element_info.get("atom_percent", 0.0)
                    mass_percent = element_info.get("mass_percent")

                    # Convert atom percent to atomic fraction (divide by 100)
                    atomic_fraction = atom_percent / 100.0

                    # Add comment if this is a recalculated value
                    comment = ""
                    if element_info.get("recalculated", False):
                        if element_info.get("added_from_formula", False):
                            comment = "  # Added from computed formula Nd_aCe_bFe_14B"
                        else:
                            comment = "  # Recalculated from formula Nd_aCe_bFe_14B"

                    # Convert mass percent to mass fraction (divide by 100) if available and requested
                    if include_mass_fraction and mass_percent is not None:
                        mass_fraction = mass_percent / 100.0
                        entry = f"""    - element: {element_symbol}
      atomic_fraction: {atomic_fraction:.6f}
      mass_fraction: {mass_fraction:.6f}{comment}"""
                    else:
                        # If mass_percent is not available or not requested, omit mass_fraction
                        entry = f"""    - element: {element_symbol}
      atomic_fraction: {atomic_fraction:.6f}{comment}"""

                    elemental_entries.append(entry)

            if elemental_entries:
                # Replace the template elemental composition with actual data
                elemental_composition = "\n".join(elemental_entries)

                # Try to match template with mass_fraction first (both commented and uncommented)
                mass_fraction_pattern = r"elemental_composition:\s*-\s*element:\s*\$\$element\$\$\s*atomic_fraction:\s*\$\$atomic_fraction\$\$\s*(#\s*)?mass_fraction:\s*\$\$mass_fraction\$\$"
                if re.search(
                    mass_fraction_pattern,
                    template_content,
                    flags=re.MULTILINE | re.DOTALL,
                ):
                    template_content = re.sub(
                        mass_fraction_pattern,
                        f"elemental_composition:\n{elemental_composition}",
                        template_content,
                        flags=re.MULTILINE | re.DOTALL,
                    )
                else:
                    # Try to match template without mass_fraction
                    atomic_fraction_pattern = r"elemental_composition:\s*-\s*element:\s*\$\$element\$\$\s*atomic_fraction:\s*\$\$atomic_fraction\$\$"
                    template_content = re.sub(
                        atomic_fraction_pattern,
                        f"elemental_composition:\n{elemental_composition}",
                        template_content,
                        flags=re.MULTILINE | re.DOTALL,
                    )
            else:
                # No verified elements found, remove elemental composition
                mass_fraction_pattern = r"elemental_composition:\s*-\s*element:\s*\$\$element\$\$\s*atomic_fraction:\s*\$\$atomic_fraction\$\$\s*(#\s*)?mass_fraction:\s*\$\$mass_fraction\$\$"
                if re.search(
                    mass_fraction_pattern,
                    template_content,
                    flags=re.MULTILINE | re.DOTALL,
                ):
                    template_content = re.sub(
                        mass_fraction_pattern,
                        "elemental_composition: []",
                        template_content,
                        flags=re.MULTILINE | re.DOTALL,
                    )
                else:
                    # Try to match template without mass_fraction
                    atomic_fraction_pattern = r"elemental_composition:\s*-\s*element:\s*\$\$element\$\$\s*atomic_fraction:\s*\$\$atomic_fraction\$\$"
                    template_content = re.sub(
                        atomic_fraction_pattern,
                        "elemental_composition: []",
                        template_content,
                        flags=re.MULTILINE | re.DOTALL,
                    )
        else:
            # No elements found, remove elemental composition
            mass_fraction_pattern = r"elemental_composition:\s*-\s*element:\s*\$\$element\$\$\s*atomic_fraction:\s*\$\$atomic_fraction\$\$\s*(#\s*)?mass_fraction:\s*\$\$mass_fraction\$\$"
            if re.search(
                mass_fraction_pattern, template_content, flags=re.MULTILINE | re.DOTALL
            ):
                template_content = re.sub(
                    mass_fraction_pattern,
                    "elemental_composition: []",
                    template_content,
                    flags=re.MULTILINE | re.DOTALL,
                )
            else:
                # Try to match template without mass_fraction
                atomic_fraction_pattern = r"elemental_composition:\s*-\s*element:\s*\$\$element\$\$\s*atomic_fraction:\s*\$\$atomic_fraction\$\$"
                template_content = re.sub(
                    atomic_fraction_pattern,
                    "elemental_composition: []",
                    template_content,
                    flags=re.MULTILINE | re.DOTALL,
                )

        # Update the short_name to include sample coordinates
        coords_str = sample_key.replace("(", "").replace(")", "").replace(",", "_")
        template_content = re.sub(
            r"short_name:\s*'\$\$NEEL-Sample-001\$\$'",
            f"short_name: 'NEEL-Sample-{coords_str}'",
            template_content,
        )

        # Update the sample_name to include sample coordinates
        template_content = re.sub(
            r"sample_name:\s*'\$\$NEEL-Sample-002\$\$'",
            f"sample_name: 'NEEL-Sample-{coords_str}'",
            template_content,
        )

        # Remove chemical formula field from template (not exported)
        template_content = re.sub(
            r"\s*chemical_formula:\s*\$\$chemical_formula\$\$\s*\n",
            "",
            template_content,
        )

        # Handle mass_fraction related modifications in the template
        if include_mass_fraction:
            # Uncomment mass_fraction lines in the template
            # Uncomment in quantities section
            template_content = re.sub(
                r"(\s*)# mass_fraction:\s*\n(\s*)#\s*type:\s*np\.float64",
                r"\1mass_fraction:\n\2  type: np.float64",
                template_content,
                flags=re.MULTILINE,
            )
            # Uncomment in data section
            template_content = re.sub(
                r"(\s*)# mass_fraction:\s*\$\$mass_fraction\$\$",
                r"\1mass_fraction: $$mass_fraction$$",
                template_content,
            )
        else:
            # Ensure mass_fraction lines are commented out
            # Comment in quantities section if not already commented
            template_content = re.sub(
                r"(\s*)mass_fraction:\s*\n(\s*)type:\s*np\.float64",
                r"\1# mass_fraction:\n\2#   type: np.float64",
                template_content,
                flags=re.MULTILINE,
            )
            # Comment in data section if not already commented
            template_content = re.sub(
                r"(\s*)mass_fraction:\s*\$\$mass_fraction\$\$",
                r"\1# mass_fraction: $$mass_fraction$$",
                template_content,
            )

        # Write the output file
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(template_content)

        print(f"Created YAML file: {output_path}")
        return True

    except Exception as e:
        print(f"Error creating YAML file {output_path}: {e}")
        return False


def process_all_samples_to_yaml(
    all_coordinate_data,
    template_path,
    output_dir,
    include_mass_fraction=False,
    verbose=True,
):
    """
    Process all samples and create individual YAML files.

    Args:
        all_coordinate_data (dict): All coordinate data from processed files
        template_path (str): Path to the template YAML file
        output_dir (str): Directory to save output YAML files
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        verbose (bool): If True, print detailed output. If False, print minimal output. (default: True)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        # Clean up existing YAML files in the output directory
        existing_yaml_files = [
            f
            for f in os.listdir(output_dir)
            if f.endswith(".yaml") or f.endswith(".yml")
        ]
        if existing_yaml_files:
            print(f"Cleaning up {len(existing_yaml_files)} existing YAML file(s)...")
            for yaml_file in existing_yaml_files:
                file_path = os.path.join(output_dir, yaml_file)
                try:
                    os.remove(file_path)
                    print(f"  Deleted: {yaml_file}")
                except Exception as e:
                    print(f"  Error deleting {yaml_file}: {e}")

    yaml_files_created = 0

    # Merge EDX and MOKE data by coordinates
    merged_data = merge_edx_and_moke_data(all_coordinate_data)

    for filename, file_data in merged_data.items():
        file_base = os.path.splitext(filename)[0]  # Remove .hdf5 extension

        for combined_group, samples in file_data.items():
            for sample_key, sample_data in samples.items():
                # Extract coordinates for the filename
                x_pos = sample_data.get("x_pos_instrument", 0.0)
                y_pos = sample_data.get("y_pos_instrument", 0.0)

                # Create clean coordinate string for filename
                x_str = str(x_pos).replace(".", "_").replace("-", "neg")
                y_str = str(y_pos).replace(".", "_").replace("-", "neg")

                # Create combined filename indicating both EDX and MOKE data types
                data_types = []
                if sample_data.get("has_edx", False):
                    data_types.append("EDX")
                if sample_data.get("has_moke", False):
                    data_types.append("MOKE")

                data_type_str = "_".join(data_types) if data_types else "UNKNOWN"
                output_filename = f"{file_base}_{data_type_str}_xpos={x_str}_ypos={y_str}.archive.yaml"
                output_path = os.path.join(output_dir, output_filename)

                # Create YAML file from template
                if create_yaml_from_template(
                    template_path,
                    output_path,
                    sample_data,
                    sample_key,
                    include_mass_fraction,
                ):
                    yaml_files_created += 1

    print(f"\nCreated {yaml_files_created} YAML files in {output_dir}")
    return yaml_files_created


def create_upload_zip(
    hdf5_files, script_dir, datasets_dir, output_dir, verbose=True, include_readme=False
):
    """
    Create a zip file containing the original HDF5 files and all generated YAML files.

    Args:
        hdf5_files (list): List of HDF5 filenames that were processed
        script_dir (str): Script directory path
        datasets_dir (str): Directory containing HDF5 files
        output_dir (str): Directory containing generated YAML files
        verbose (bool): If True, print detailed output. If False, print minimal output.
        include_readme (bool): If True, include a README.txt file in the zip (default: False)

    Returns:
        str: Path to the created zip file, or None if creation failed
    """
    try:
        # Create uploads directory if it doesn't exist
        uploads_dir = os.path.join(script_dir, "uploads")
        if not os.path.exists(uploads_dir):
            os.makedirs(uploads_dir)
            if verbose:
                print(f"Created uploads directory: {uploads_dir}")

        # Generate timestamp for unique zip filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Determine zip filename based on processed files
        if len(hdf5_files) == 1:
            base_name = os.path.splitext(hdf5_files[0])[0]
            zip_filename = f"{base_name}_{timestamp}.zip"
        else:
            zip_filename = f"NEEL_data_batch_{timestamp}.zip"

        zip_path = os.path.join(uploads_dir, zip_filename)

        # Create the zip file
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            files_added = 0

            # Add HDF5 files
            for hdf5_file in hdf5_files:
                hdf5_path = os.path.join(datasets_dir, hdf5_file)
                if os.path.exists(hdf5_path):
                    zipf.write(hdf5_path, f"datasets/{hdf5_file}")
                    files_added += 1
                    if verbose:
                        print(f"Added HDF5 file: {hdf5_file}")

            # Add YAML files if output directory exists
            if os.path.exists(output_dir):
                yaml_files = [
                    f for f in os.listdir(output_dir) if f.endswith((".yaml", ".yml"))
                ]
                for yaml_file in yaml_files:
                    yaml_path = os.path.join(output_dir, yaml_file)
                    zipf.write(yaml_path, f"{yaml_file}")
                    files_added += 1
                    if verbose:
                        print(f"Added YAML file: {yaml_file}")

            # Add README file with processing information (optional)
            if include_readme:
                readme_content = f"""# NEEL Data Processing Results

Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Contents:
- datasets/: Original HDF5 files
- generated_schemas/: Generated NOMAD schema YAML files

## Processed Files:
{chr(10).join([f"- {f}" for f in hdf5_files])}

## YAML Files Generated: {len([f for f in os.listdir(output_dir) if f.endswith((".yaml", ".yml"))]) if os.path.exists(output_dir) else 0}

This archive contains the original HDF5 data files and the corresponding NOMAD schema files
generated from EDX analysis data with atomic composition information.
"""
                zipf.writestr("README.txt", readme_content)
                files_added += 1
                if verbose:
                    print("Added README.txt")

        if verbose:
            print(f"\nCreated zip file: {zip_filename}")
            print(f"Location: {zip_path}")
            print(f"Total files added: {files_added}")
        else:
            print(f"Created zip file: {zip_filename} ({files_added} files)")

        return zip_path

    except Exception as e:
        if verbose:
            print(f"Error creating zip file: {e}")
        return None


def process_single_file(
    filename="NdCeFeB_2-5.hdf5",
    verbose=True,
    include_mass_fraction=False,
    create_zip=False,
    include_readme=False,
):
    """
    Process a single HDF5 file.

    Args:
        filename (str): Name of the HDF5 file to process (default: 'NdCeFeB_2-5.hdf5')
        verbose (bool): If True, print detailed output. If False, print only summary.
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        create_zip (bool): Whether to create a zip file with HDF5 and YAML files (default: False)
        include_readme (bool): Whether to include README.txt in the zip file (default: False)
    """
    main(
        single_file=filename,
        verbose=verbose,
        include_mass_fraction=include_mass_fraction,
        create_zip=create_zip,
        include_readme=include_readme,
    )


def process_all_files(
    verbose=True, include_mass_fraction=False, create_zip=False, include_readme=False
):
    """
    Process all HDF5 files in the datasets directory.

    Args:
        verbose (bool): If True, print detailed output. If False, print only summary.
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        create_zip (bool): Whether to create a zip file with HDF5 and YAML files (default: False)
        include_readme (bool): Whether to include README.txt in the zip file (default: False)
    """
    main(
        single_file="",
        verbose=verbose,
        include_mass_fraction=include_mass_fraction,
        create_zip=create_zip,
        include_readme=include_readme,
    )  # Empty string to force processing all files


def main(
    single_file=None,
    verbose=True,
    include_mass_fraction=False,
    create_zip=False,
    include_readme=False,
):
    """
    Main function to process HDF5 files in the datasets subfolder using neel_data_vis.

    Args:
        single_file (str, optional): Specific filename to process. If None, processes all HDF5 files.
                                   Default: 'NdCeFeB_2-5.hdf5'
        verbose (bool): If True, print detailed output. If False, print only summary.
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        create_zip (bool): Whether to create a zip file with HDF5 and YAML files (default: False)
        include_readme (bool): Whether to include README.txt in the zip file (default: False)
    """
    # Set default single file if none specified
    if single_file is None:
        # single_file = "NdCeFeB_2-5.hdf5"
        single_file = "SmFeV_3378.hdf5"

    # Define the datasets directory relative to the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    datasets_dir = os.path.join(script_dir, "datasets")

    if not os.path.exists(datasets_dir):
        print(f"Datasets directory not found: {datasets_dir}")
        return

    # Determine which files to process
    if single_file and single_file != "":
        # Process only the specified file
        single_file_path = os.path.join(datasets_dir, single_file)
        if os.path.exists(single_file_path) and single_file.endswith(".hdf5"):
            hdf5_files = [single_file]
            print(f"Processing single file: {single_file}")
        else:
            print(f"Single file not found or not an HDF5 file: {single_file}")
            print("Available HDF5 files in datasets directory:")
            available_files = [
                f for f in os.listdir(datasets_dir) if f.endswith(".hdf5")
            ]
            if available_files:
                for i, filename in enumerate(available_files, 1):
                    print(f"  {i}. {filename}")
            else:
                print("  No HDF5 files found.")
            return
    else:
        # Get all HDF5 files in the datasets directory
        hdf5_files = [f for f in os.listdir(datasets_dir) if f.endswith(".hdf5")]

        if not hdf5_files:
            print("No HDF5 files found in the datasets directory.")
            return

        print(f"Processing all files - Found {len(hdf5_files)} HDF5 file(s):")
        for i, filename in enumerate(hdf5_files, 1):
            print(f"  {i}. {filename}")

    # Process each HDF5 file
    # Dictionary to store all coordinate data
    all_coordinate_data = {}
    total_groups_processed = 0
    total_coordinates_found = 0
    total_matches = 0
    total_mismatches = 0
    total_elements_found = 0
    total_element_matches = 0
    total_element_mismatches = 0

    for filename in hdf5_files:
        file_path = os.path.join(datasets_dir, filename)
        if verbose:
            print(f"\n{'=' * 60}")
            print(f"Processing file: {filename}")
            print(f"{'=' * 60}")
        else:
            print(f"Processing file: {filename}")

        # Use neel_data_vis xarray-based extraction
        if verbose:
            print("\nExtracting data using neel_data_vis (xarray method)...")

        xarray_data = extract_all_data_with_xarray(file_path, verbose)

        if xarray_data is not None and filename in xarray_data:
            # Successfully extracted data using xarray method
            file_data = xarray_data[filename]

            # Count statistics from xarray extraction
            for group_name, samples in file_data.items():
                total_groups_processed += 1
                total_coordinates_found += len(samples)

                for sample_key, sample_data in samples.items():
                    # Count coordinate matches (always True for xarray)
                    if sample_data.get("x_match", False) and sample_data.get(
                        "y_match", False
                    ):
                        total_matches += 1
                    else:
                        total_mismatches += 1

                    # Count elements
                    elements = sample_data.get("elements", {})
                    total_elements_found += len(elements)
                    for elem_info in elements.values():
                        if elem_info.get("symbol_match", False):
                            total_element_matches += 1
                        else:
                            total_element_mismatches += 1

            all_coordinate_data[filename] = file_data

            if verbose:
                print("\n✓ Successfully extracted data using xarray method")
                print(f"  Groups processed: {len(file_data)}")
                print(f"  Sample positions: {sum(len(s) for s in file_data.values())}")
        else:
            if verbose:
                print("\nxarray extraction returned no data for this file.")

    # Print final summary
    if verbose:
        print(f"\n{'=' * 80}")
        print("FINAL SUMMARY")
        print(f"{'=' * 80}")
        print("Extraction method: neel_data_vis (xarray)")
        print(f"Files processed: {len(hdf5_files)}")
        print(f"Groups processed: {total_groups_processed}")
        print(f"Coordinate pairs found: {total_coordinates_found}")
        print(f"Coordinate matches: {total_matches}")
        print(f"Coordinate mismatches: {total_mismatches}")
        print(f"Elements found: {total_elements_found}")
        print(f"Element symbol matches: {total_element_matches}")
        print(f"Element symbol mismatches: {total_element_mismatches}")

        if total_coordinates_found > 0:
            coord_success_rate = (total_matches / total_coordinates_found) * 100
            print(f"Coordinate success rate: {coord_success_rate:.1f}%")

            if total_matches == total_coordinates_found:
                print("✓ ALL COORDINATE CHECKS PASSED!")
            else:
                print("⚠ Some coordinate mismatches found.")
        else:
            print("No coordinate data found to verify.")

        if total_elements_found > 0:
            element_success_rate = (total_element_matches / total_elements_found) * 100
            print(f"Element verification success rate: {element_success_rate:.1f}%")

            if total_element_matches == total_elements_found:
                print("✓ ALL ELEMENT VERIFICATIONS PASSED!")
            else:
                print("⚠ Some element symbol mismatches found.")
        else:
            print("No element data found to verify.")

        print(f"{'=' * 80}")
    else:
        # Short summary
        print(
            f"\nSUMMARY [xarray]: Processed {len(hdf5_files)} file(s), found {total_coordinates_found} coordinate pairs, {total_elements_found} elements"
        )
        if total_coordinates_found > 0:
            coord_success_rate = (total_matches / total_coordinates_found) * 100
            print(
                f"Coordinate success: {total_matches}/{total_coordinates_found} ({coord_success_rate:.1f}%)"
            )
        if total_elements_found > 0:
            element_success_rate = (total_element_matches / total_elements_found) * 100
            print(
                f"Element verification: {total_element_matches}/{total_elements_found} ({element_success_rate:.1f}%)"
            )

    # Generate YAML files from template if we have data
    if all_coordinate_data:
        if verbose:
            print("\n" + "=" * 80)
            print("GENERATING YAML FILES FROM TEMPLATE")
            print("=" * 80)
        else:
            print("\nGenerating YAML files...")

        # Define paths
        template_path = os.path.join(script_dir, "NEEL_template.archive.yaml")
        output_dir = os.path.join(script_dir, "generated_schemas")

        if os.path.exists(template_path):
            yaml_files_created = process_all_samples_to_yaml(
                all_coordinate_data,
                template_path,
                output_dir,
                include_mass_fraction,
                verbose,
            )
            print(f"Successfully generated {yaml_files_created} YAML schema files!")

            # Create upload zip file (optional)
            if create_zip:
                if verbose:
                    print("\n" + "=" * 40)
                    print("CREATING UPLOAD ZIP FILE")
                    print("=" * 40)
                else:
                    print("\nCreating upload zip file...")

                zip_file_path = create_upload_zip(
                    hdf5_files,
                    script_dir,
                    datasets_dir,
                    output_dir,
                    verbose,
                    include_readme,
                )
                if zip_file_path:
                    if verbose:
                        print(
                            f"✓ Upload zip file ready: {os.path.basename(zip_file_path)}"
                        )
                    else:
                        print(f"✓ Upload zip ready: {os.path.basename(zip_file_path)}")
                else:
                    print("✗ Failed to create upload zip file.")
            else:
                if verbose:
                    print("\nSkipping zip file creation (create_zip=False)")
        else:
            print(f"Template file not found: {template_path}")
    else:
        if verbose:
            print("\nNo data available to generate YAML files.")


def merge_edx_and_moke_data(all_coordinate_data):
    """
    Merge EDX and MOKE data for samples at the same coordinates.

    Args:
        all_coordinate_data (dict): Dictionary containing both EDX and MOKE data

    Returns:
        dict: Merged data organized by coordinates
    """
    merged_data = {}

    for filename, file_data in all_coordinate_data.items():
        file_base = os.path.splitext(filename)[0]
        merged_data[filename] = {}

        # Collect all coordinate data from EDX and MOKE groups
        coordinate_map = {}  # Maps (x, y) coordinates to combined data

        for group_name, samples in file_data.items():
            for sample_key, sample_data in samples.items():
                # Get coordinates (prefer instrument coordinates)
                x_pos = sample_data.get("x_pos_instrument")
                y_pos = sample_data.get("y_pos_instrument")

                # Fallback to group-specific coordinates if instrument coordinates not available
                if x_pos is None:
                    x_pos = sample_data.get("x_pos_EDX") or sample_data.get(
                        "x_pos_MOKE", 0.0
                    )
                if y_pos is None:
                    y_pos = sample_data.get("y_pos_EDX") or sample_data.get(
                        "y_pos_MOKE", 0.0
                    )

                coord_key = (x_pos, y_pos)

                # Initialize coordinate entry if not exists
                if coord_key not in coordinate_map:
                    coordinate_map[coord_key] = {
                        "x_pos_instrument": x_pos,
                        "y_pos_instrument": y_pos,
                        "x_pos_unit": sample_data.get("x_pos_unit", ""),
                        "y_pos_unit": sample_data.get("y_pos_unit", ""),
                        "has_edx": False,
                        "has_moke": False,
                        "sample_key": sample_key,
                        "group_names": [],
                    }

                # Merge data based on group type
                if "edx" in group_name:
                    coordinate_map[coord_key]["has_edx"] = True
                    coordinate_map[coord_key]["elements"] = sample_data.get(
                        "elements", {}
                    )
                    coordinate_map[coord_key]["x_pos_EDX"] = sample_data.get(
                        "x_pos_EDX"
                    )
                    coordinate_map[coord_key]["y_pos_EDX"] = sample_data.get(
                        "y_pos_EDX"
                    )
                    coordinate_map[coord_key]["x_match"] = sample_data.get(
                        "x_match", False
                    )
                    coordinate_map[coord_key]["y_match"] = sample_data.get(
                        "y_match", False
                    )

                elif "moke" in group_name:
                    coordinate_map[coord_key]["has_moke"] = True
                    coordinate_map[coord_key]["moke_data"] = sample_data.get(
                        "moke_data", {}
                    )
                    coordinate_map[coord_key]["x_pos_MOKE"] = sample_data.get(
                        "x_pos_MOKE"
                    )
                    coordinate_map[coord_key]["y_pos_MOKE"] = sample_data.get(
                        "y_pos_MOKE"
                    )
                    # Update match status if not already set by EDX
                    if not coordinate_map[coord_key].get("x_match", False):
                        coordinate_map[coord_key]["x_match"] = sample_data.get(
                            "x_match", False
                        )
                    if not coordinate_map[coord_key].get("y_match", False):
                        coordinate_map[coord_key]["y_match"] = sample_data.get(
                            "y_match", False
                        )

                coordinate_map[coord_key]["group_names"].append(group_name)

        # Convert coordinate map back to group structure for compatibility
        combined_group_name = f"{file_base}_Combined"
        merged_data[filename][combined_group_name] = {}

        for coord_key, merged_sample_data in coordinate_map.items():
            x_pos, y_pos = coord_key
            # Create a new sample key based on coordinates
            sample_key = f"({x_pos},{y_pos})"
            merged_data[filename][combined_group_name][sample_key] = merged_sample_data

    return merged_data


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Process NEEL HDF5 datasets and generate NOMAD schema YAML files"
    )
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        default=None,
        help="Process a specific HDF5 file from the datasets folder (e.g., 'SmFeV_3378.hdf5'). "
        "If not specified, all HDF5 files in the datasets folder will be processed.",
    )
    parser.add_argument(
        "-m",
        "--mass-fraction",
        action="store_true",
        help="Include mass_fraction in the output YAML files",
    )
    parser.add_argument(
        "-z",
        "--create-zip",
        action="store_true",
        help="Create a zip file containing all HDF5 and generated YAML files",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress verbose output"
    )
    parser.add_argument(
        "-r",
        "--include-readme",
        action="store_true",
        help="Include README.txt in the zip file (only applies when -z is used)",
    )

    args = parser.parse_args()

    # Call main with parsed arguments
    main(
        single_file=args.file,
        verbose=not args.quiet,
        include_mass_fraction=args.mass_fraction,
        create_zip=args.create_zip,
        include_readme=args.include_readme,
    )

    # USAGE EXAMPLES:
    # ================================================================================
    #
    # 1. Process all HDF5 files in datasets/ folder (default)
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py
    #
    # 2. Process a specific file
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -f SmFeV_3378.hdf5
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py --file SmFeV_3378.hdf5
    #
    # 3. Include mass fractions in output YAML files
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -m
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py --mass-fraction
    #
    # 4. Create a zip file for NOMAD upload
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -z
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py --create-zip
    #
    # 5. Suppress verbose output (quiet mode)
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -q
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py --quiet
    #
    # 6. Combine options:
    #    # Process specific file with mass fractions
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -f SmFeV_3378.hdf5 -m
    #
    #    # Process specific file and create zip
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -f SmFeV_3378.hdf5 -z
    #
    #    # Process all files with mass fractions and zip
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -m -z
    #
    #    # Process specific file with quiet output
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -f SmFeV_3378.hdf5 -q
    #
    #    # All options combined
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -f SmFeV_3378.hdf5 -m -z -q
    #
    # 7. Show help message with all available options
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py -h
    #    python run_convert_NEEL_data_from_hdf5_and_create_schemas.py --help
    #
    # ================================================================================
    #
    # EXTRACTION METHOD:
    # This script uses neel_data_vis (xarray-based) for data extraction, which
    # provides clean access to the HDF5 data through xarray's coordinate-based indexing.
    #
    # AVAILABLE OPTIONS:
    # -f, --file FILE       : Process a specific HDF5 file from datasets/ folder
    # -m, --mass-fraction   : Include mass_fraction in output YAML files
    # -z, --create-zip      : Create a zip file containing HDF5 and YAML files
    # -q, --quiet           : Suppress verbose output (minimal mode)
    # -h, --help            : Show help message
    #
    # ===============================================================================
