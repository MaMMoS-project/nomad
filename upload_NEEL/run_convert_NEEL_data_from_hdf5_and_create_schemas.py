import h5py
import os
import re
import zipfile
from datetime import datetime

# Note: Chemical formula calculation is implemented but not exported to YAML
#       (atomic fractions are recalculated including B and exported instead)
# TODO: after new dataset was provided 20250626: check "HT_type" which should be e.g. "edx" as an attribute of the group named with "EDX"
# TODO: after new description was provided 20250626: insert new description of the data set (raw / filtered data)
# TODO: after new .pdfs are accessible 20250626: reference to the chada-docu as .pdf on mammos-project on github


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


def find_edx_groups(file_path):
    """
    Find first-level groups containing 'EDX' in their name within an HDF5 file.

    Args:
        file_path (str): Path to the HDF5 file

    Returns:
        list: List of first-level group paths containing 'EDX' in their name
    """
    edx_groups = []

    try:
        with h5py.File(file_path, "r") as f:
            # Only check first-level groups (direct children of root)
            for key in f.keys():
                item = f[key]
                if isinstance(item, h5py.Group) and "EDX" in key:
                    edx_groups.append(key)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return []

    return edx_groups


def list_group_entries(file_path, group_path, verbose=True):
    """
    List all entries (groups and datasets) within a specific group.

    Args:
        file_path (str): Path to the HDF5 file
        group_path (str): Path to the group within the HDF5 file
        verbose (bool): If True, print detailed output. If False, print minimal output.
    """
    try:
        with h5py.File(file_path, "r") as f:
            if group_path in f:
                group = f[group_path]
                if verbose:
                    print(f"\nEntries in group '{group_path}':")
                    print("-" * 50)

                    for key in group.keys():
                        item = group[key]
                        if isinstance(item, h5py.Group):
                            print(f"  [GROUP]   {key}")
                        elif isinstance(item, h5py.Dataset):
                            shape = item.shape
                            dtype = item.dtype
                            print(f"  [DATASET] {key} - Shape: {shape}, Type: {dtype}")
                        else:
                            print(f"  [OTHER]   {key}")

                    print(f"\nTotal entries: {len(group.keys())}")
                return len(group.keys())
            else:
                if verbose:
                    print(f"Group '{group_path}' not found in file.")
                return 0
    except Exception as e:
        if verbose:
            print(f"Error accessing group {group_path} in file {file_path}: {e}")
        return 0


def extract_coordinates_from_name(name):
    """
    Extract x and y coordinates from a group name like "(-10,25)".

    Args:
        name (str): Group name containing coordinates in parentheses

    Returns:
        tuple: (x, y) as floats, or (None, None) if pattern not found
    """
    # Pattern to match coordinates like (-10,25), (5.5, -3.2), etc.
    pattern = r"\((-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)\)"
    match = re.search(pattern, name)

    if match:
        x = float(match.group(1))
        y = float(match.group(2))
        return x, y
    return None, None


def get_instrument_positions(file_path, group_path, verbose=True):
    """
    Get x_pos and y_pos values and their units from the 'instrument' subfolder.

    Args:
        file_path (str): Path to the HDF5 file
        group_path (str): Path to the parent group
        verbose (bool): If True, print detailed output. If False, print minimal output.

    Returns:
        tuple: (x_pos, y_pos, x_unit, y_unit) or (None, None, None, None) if not found
    """
    try:
        with h5py.File(file_path, "r") as f:
            instrument_path = f"{group_path}/instrument"

            if instrument_path in f:
                if verbose:
                    print(f"Reading instrument data from: {instrument_path}")
                instrument_group = f[instrument_path]

                x_pos = None
                y_pos = None
                x_unit = None
                y_unit = None

                if "x_pos" in instrument_group:
                    if verbose:
                        print("Found 'x_pos' in instrument group.")
                    x_pos_data = instrument_group["x_pos"]
                    # For scalar 64-bit integer dataset
                    x_pos = float(x_pos_data[()])
                    if verbose:
                        print(f"  x_pos value: {x_pos}")

                    # Extract unit attribute if it exists
                    if "units" in x_pos_data.attrs:
                        x_unit = x_pos_data.attrs["units"]
                        if isinstance(x_unit, bytes):
                            x_unit = x_unit.decode("utf-8")
                        if verbose:
                            print(f"  x_pos unit: {x_unit}")

                if "y_pos" in instrument_group:
                    if verbose:
                        print("Found 'y_pos' in instrument group.")
                    y_pos_data = instrument_group["y_pos"]
                    # For scalar 64-bit integer dataset
                    y_pos = float(y_pos_data[()])
                    if verbose:
                        print(f"  y_pos value: {y_pos}")

                    # Extract unit attribute if it exists
                    if "units" in y_pos_data.attrs:
                        y_unit = y_pos_data.attrs["units"]
                        if isinstance(y_unit, bytes):
                            y_unit = y_unit.decode("utf-8")
                        if verbose:
                            print(f"  y_pos unit: {y_unit}")

                return x_pos, y_pos, x_unit, y_unit

    except Exception as e:
        if verbose:
            print(f"Error reading instrument data from {group_path}: {e}")

    return None, None, None, None


def process_edx_coordinates(file_path, edx_group, verbose=True):
    """
    Process EDX group to extract coordinates from sub-group names and verify against instrument data.

    Args:
        file_path (str): Path to the HDF5 file
        edx_group (str): Path to the EDX group
        verbose (bool): If True, print detailed output. If False, print minimal output.

    Returns:
        dict: Dictionary with coordinate data and verification results
    """
    EDX_sample_coordinate_data = {}

    try:
        with h5py.File(file_path, "r") as f:
            if edx_group in f:
                group = f[edx_group]

                for key in group.keys():
                    item = group[key]
                    if isinstance(item, h5py.Group):
                        # Extract coordinates from group name
                        x_name, y_name = extract_coordinates_from_name(key)
                        if verbose:
                            print(
                                f"Processing group: {key} with coordinates ({x_name}, {y_name})"
                            )
                        if x_name is not None and y_name is not None:
                            # This is a coordinate group, rename variable to sample_key
                            sample_key = key
                            # Get instrument positions
                            if verbose:
                                print(f"  - Found coordinates: ({x_name}, {y_name})")
                            # Construct sample-group path
                            sample_group_path = f"{edx_group}/{sample_key}"
                            if verbose:
                                print(
                                    f"  - Constructed sample-group path: {sample_group_path}"
                                )
                            x_instrument, y_instrument, x_unit, y_unit = (
                                get_instrument_positions(
                                    file_path, sample_group_path, verbose
                                )
                            )

                            # Get element data from results sub-group
                            element_data = get_element_data(
                                file_path, sample_group_path, verbose
                            )

                            # Store data
                            EDX_sample_coordinate_data[sample_key] = {
                                "x_pos_EDX": x_name,
                                "y_pos_EDX": y_name,
                                "x_pos_instrument": x_instrument,
                                "y_pos_instrument": y_instrument,
                                "x_pos_unit": x_unit,
                                "y_pos_unit": y_unit,
                                "elements": element_data,
                                "x_match": x_instrument is not None
                                and abs(x_name - x_instrument) < 1e-6,
                                "y_match": y_instrument is not None
                                and abs(y_name - y_instrument) < 1e-6,
                            }

    except Exception as e:
        if verbose:
            print(f"Error processing EDX coordinates for {edx_group}: {e}")

    return EDX_sample_coordinate_data


def get_element_data(file_path, sample_group_path, verbose=True):
    """
    Extract chemical elements from the 'results' sub-group and verify element names against atomic numbers.
    Only includes elements that have an 'AtomPercent' field.

    Args:
        file_path (str): Path to the HDF5 file
        sample_group_path (str): Path to the sample group
        verbose (bool): If True, print detailed output. If False, print minimal output.

    Returns:
        dict: Dictionary with element data and verification results (only for elements with AtomPercent)
    """
    element_data = {}

    try:
        with h5py.File(file_path, "r") as f:
            results_path = f"{sample_group_path}/results"

            if results_path in f:
                if verbose:
                    print(f"Reading results data from: {results_path}")
                results_group = f[results_path]

                for key in results_group.keys():
                    item = results_group[key]
                    if isinstance(item, h5py.Group) and key.startswith("Element "):
                        # Extract element symbol from group name (e.g., 'Element C' -> 'C')
                        element_symbol = key.replace("Element ", "").strip()
                        if verbose:
                            print(
                                f"Found element group: {key} -> Symbol: {element_symbol}"
                            )

                        # Check if there's an 'Element' dataset with atomic number
                        element_group = item
                        if "Element" in element_group:
                            element_dataset = element_group["Element"]
                            atomic_number = float(element_dataset[()])
                            if verbose:
                                print(f"  Atomic number from dataset: {atomic_number}")

                            # Verify if element symbol matches atomic number
                            expected_symbol = get_element_symbol_from_atomic_number(
                                int(atomic_number)
                            )
                            symbol_match = expected_symbol == element_symbol

                            # Get AtomPercent value if available
                            atom_percent = get_atom_percent_for_element(element_group)

                            # Get MassPercent value if available
                            mass_percent = get_mass_percent_for_element(element_group)

                            # Only include elements that have AtomPercent field
                            if atom_percent is not None:
                                element_data[element_symbol] = {
                                    "group_name": key,
                                    "element_symbol": element_symbol,
                                    "atomic_number": int(atomic_number),
                                    "expected_symbol": expected_symbol,
                                    "symbol_match": symbol_match,
                                    "atom_percent": atom_percent,
                                    "mass_percent": mass_percent,  # Can be None if not available
                                }

                                if verbose:
                                    print(
                                        f"  Expected symbol for atomic number {int(atomic_number)}: {expected_symbol}"
                                    )
                                    print(
                                        f"  Symbol match: {'✓' if symbol_match else '✗'}"
                                    )
                                    print(f"  AtomPercent: {atom_percent}")
                                    if mass_percent is not None:
                                        print(f"  MassPercent: {mass_percent}")
                                    else:
                                        print("  MassPercent: Not available")
                            else:
                                if verbose:
                                    print(
                                        f"  Skipping element {element_symbol}: No AtomPercent field found"
                                    )
            else:
                if verbose:
                    print(f"No 'results' group found in {sample_group_path}")

    except Exception as e:
        if verbose:
            print(f"Error reading element data from {sample_group_path}: {e}")

    # Check for Nd-Ce-Fe-B compound formula if verbose
    if verbose and element_data:
        is_nd_ce_fe_only, nd_fraction, ce_fraction, fe_fraction = (
            check_for_nd_ce_fe_only(element_data)
        )
        if is_nd_ce_fe_only:
            print("  Detected Nd-Ce-Fe composition for compound formula calculation")
            a, b, c = compute_stoichiometric_coefficients_from_fractions(
                nd_fraction, ce_fraction
            )
            if a is not None and b is not None and c is not None:
                compound_formula = f"Nd_{a:.2f}Ce_{b:.2f}Fe_{c:.2f}B"
                print(f"  Computed compound formula: {compound_formula}")
        else:
            present_elements = list(element_data.keys())
            print(
                f"  Elements present: {present_elements} - not suitable for Nd-Ce-Fe-B formula"
            )

    return element_data


def get_atom_percent_for_element(element_group):
    """
    Get AtomPercent value from an element group if present.

    Args:
        element_group (h5py.Group): HDF5 group for the element

    Returns:
        float or None: AtomPercent value if found, None otherwise
    """
    try:
        if "AtomPercent" in element_group:
            atom_percent_dataset = element_group["AtomPercent"]
            atom_percent = float(atom_percent_dataset[()])
            return atom_percent
        else:
            return None
    except Exception as e:
        print(f"  Error reading AtomPercent: {e}")
        return None


def get_mass_percent_for_element(element_group):
    """
    Get MassPercent value from an element group if present.

    Args:
        element_group (h5py.Group): HDF5 group for the element

    Returns:
        float or None: MassPercent value if found, None otherwise
    """
    try:
        if "MassPercent" in element_group:
            mass_percent_dataset = element_group["MassPercent"]
            mass_percent = float(mass_percent_dataset[()])
            return mass_percent
        else:
            return None
    except Exception as e:
        print(f"  Error reading MassPercent: {e}")
        return None


def get_element_symbol_from_atomic_number(atomic_number):
    """
    Get element symbol from atomic number.

    Args:
        atomic_number (int): Atomic number

    Returns:
        str: Element symbol or 'Unknown' if not found
    """
    # Periodic table mapping (first 118 elements)
    periodic_table = {
        1: "H",
        2: "He",
        3: "Li",
        4: "Be",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "Ne",
        11: "Na",
        12: "Mg",
        13: "Al",
        14: "Si",
        15: "P",
        16: "S",
        17: "Cl",
        18: "Ar",
        19: "K",
        20: "Ca",
        21: "Sc",
        22: "Ti",
        23: "V",
        24: "Cr",
        25: "Mn",
        26: "Fe",
        27: "Co",
        28: "Ni",
        29: "Cu",
        30: "Zn",
        31: "Ga",
        32: "Ge",
        33: "As",
        34: "Se",
        35: "Br",
        36: "Kr",
        37: "Rb",
        38: "Sr",
        39: "Y",
        40: "Zr",
        41: "Nb",
        42: "Mo",
        43: "Tc",
        44: "Ru",
        45: "Rh",
        46: "Pd",
        47: "Ag",
        48: "Cd",
        49: "In",
        50: "Sn",
        51: "Sb",
        52: "Te",
        53: "I",
        54: "Xe",
        55: "Cs",
        56: "Ba",
        57: "La",
        58: "Ce",
        59: "Pr",
        60: "Nd",
        61: "Pm",
        62: "Sm",
        63: "Eu",
        64: "Gd",
        65: "Tb",
        66: "Dy",
        67: "Ho",
        68: "Er",
        69: "Tm",
        70: "Yb",
        71: "Lu",
        72: "Hf",
        73: "Ta",
        74: "W",
        75: "Re",
        76: "Os",
        77: "Ir",
        78: "Pt",
        79: "Au",
        80: "Hg",
        81: "Tl",
        82: "Pb",
        83: "Bi",
        84: "Po",
        85: "At",
        86: "Rn",
        87: "Fr",
        88: "Ra",
        89: "Ac",
        90: "Th",
        91: "Pa",
        92: "U",
        93: "Np",
        94: "Pu",
        95: "Am",
        96: "Cm",
        97: "Bk",
        98: "Cf",
        99: "Es",
        100: "Fm",
        101: "Md",
        102: "No",
        103: "Lr",
        104: "Rf",
        105: "Db",
        106: "Sg",
        107: "Bh",
        108: "Hs",
        109: "Mt",
        110: "Ds",
        111: "Rg",
        112: "Cn",
        113: "Nh",
        114: "Fl",
        115: "Mc",
        116: "Lv",
        117: "Ts",
        118: "Og",
    }

    return periodic_table.get(atomic_number, "Unknown")


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
                m0 = 4 * 3.141592653589793 * 1e-7  # Vacuum permeability in T*m/A
                coercivity_value_in_A_per_m = (
                    float(coercivity_value) / m0
                )  # Convert T to A/m

            # Replace coercivity placeholder with actual value
            template_content = template_content.replace(
                "$$coercivity$$", str(coercivity_value_in_A_per_m)
            )
        else:
            print("No MOKE data available, removing coercivity field from template")
            # Remove the coercivity field entirely if no MOKE data
            template_content = re.sub(
                r"\s*CoercivityBHcExternal:\s*\$\$coercivity\$\$\s*\n",
                "\n",
                template_content,
            )
            # Also remove from schema definition if present
            template_content = re.sub(
                r"\s*CoercivityBHcExternal:\s*\n\s*type:\s*np\.float64\s*\n\s*unit:\s*A/m\s*\n\s*description:\s*\'Coercivity\s*from\s*MOKE\s*measurements\s*\(optional\)\'\s*\n",
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
    all_coordinate_data, template_path, output_dir, include_mass_fraction=False
):
    """
    Process all samples and create individual YAML files.

    Args:
        all_coordinate_data (dict): All coordinate data from processed files
        template_path (str): Path to the template YAML file
        output_dir (str): Directory to save output YAML files
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
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


def create_upload_zip(hdf5_files, script_dir, datasets_dir, output_dir, verbose=True):
    """
    Create a zip file containing the original HDF5 files and all generated YAML files.

    Args:
        hdf5_files (list): List of HDF5 filenames that were processed
        script_dir (str): Script directory path
        datasets_dir (str): Directory containing HDF5 files
        output_dir (str): Directory containing generated YAML files
        verbose (bool): If True, print detailed output. If False, print minimal output.

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
                    zipf.write(yaml_path, f"generated_schemas/{yaml_file}")
                    files_added += 1
                    if verbose:
                        print(f"Added YAML file: {yaml_file}")

            # Add README file with processing information
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
):
    """
    Process a single HDF5 file.

    Args:
        filename (str): Name of the HDF5 file to process (default: 'NdCeFeB_2-5.hdf5')
        verbose (bool): If True, print detailed output. If False, print only summary.
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        create_zip (bool): Whether to create a zip file with HDF5 and YAML files (default: False)
    """
    main(
        single_file=filename,
        verbose=verbose,
        include_mass_fraction=include_mass_fraction,
        create_zip=create_zip,
    )


def process_all_files(verbose=True, include_mass_fraction=False, create_zip=False):
    """
    Process all HDF5 files in the datasets directory.

    Args:
        verbose (bool): If True, print detailed output. If False, print only summary.
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        create_zip (bool): Whether to create a zip file with HDF5 and YAML files (default: False)
    """
    main(
        single_file="",
        verbose=verbose,
        include_mass_fraction=include_mass_fraction,
        create_zip=create_zip,
    )  # Empty string to force processing all files


def main(single_file=None, verbose=True, include_mass_fraction=False, create_zip=False):
    """
    Main function to process HDF5 files in the datasets subfolder.

    Args:
        single_file (str, optional): Specific filename to process. If None, processes all HDF5 files.
                                   Default: 'NdCeFeB_2-5.hdf5'
        verbose (bool): If True, print detailed output. If False, print only summary.
        include_mass_fraction (bool): Whether to include mass_fraction in the output (default: False)
        create_zip (bool): Whether to create a zip file with HDF5 and YAML files (default: False)
    """
    # Set default single file if none specified
    if single_file is None:
        single_file = "NdCeFeB_2-5.hdf5"

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

        # Find groups containing 'EDX'
        edx_groups = find_edx_groups(file_path)

        if edx_groups:
            if verbose:
                print(f"\nFound {len(edx_groups)} group(s) containing 'EDX':")
                for group_path in edx_groups:
                    print(f"  - {group_path}")

            # Process coordinate data for each EDX group
            file_coordinate_data = {}
            for group_path in edx_groups:
                total_groups_processed += 1

                # List entries for each EDX group
                list_group_entries(file_path, group_path, verbose)

                # Process coordinates
                EDX_sample_coordinate_data = process_edx_coordinates(
                    file_path, group_path, verbose
                )
                if EDX_sample_coordinate_data:
                    file_coordinate_data[group_path] = EDX_sample_coordinate_data
                    total_coordinates_found += len(EDX_sample_coordinate_data)

                    # Print coordinate analysis
                    if verbose:
                        print(f"\n--- Coordinate Analysis for {group_path} ---")
                    for sample_group, data in EDX_sample_coordinate_data.items():
                        if verbose:
                            print(f"Sample-group: {sample_group}")
                            print(f"  x_pos_EDX: {data['x_pos_EDX']}")
                            print(f"  y_pos_EDX: {data['y_pos_EDX']}")
                            print(f"  x_pos_instrument: {data['x_pos_instrument']}")
                            print(f"  y_pos_instrument: {data['y_pos_instrument']}")
                            print(f"  x_pos_unit: {data['x_pos_unit']}")
                            print(f"  y_pos_unit: {data['y_pos_unit']}")

                        x_status = "✓" if data["x_match"] else "✗"
                        y_status = "✓" if data["y_match"] else "✗"
                        if verbose:
                            print(f"  X coordinate match: {x_status}")
                            print(f"  Y coordinate match: {y_status}")

                        # Display element data
                        if data.get("elements"):
                            if verbose:
                                print(f"  Elements found: {len(data['elements'])}")
                            total_elements_found += len(data["elements"])
                            for element_symbol, element_info in data[
                                "elements"
                            ].items():
                                element_status = (
                                    "✓" if element_info["symbol_match"] else "✗"
                                )
                                if verbose:
                                    print(
                                        f"    {element_symbol} (Z={element_info['atomic_number']}): {element_status}"
                                    )
                                if element_info["symbol_match"]:
                                    total_element_matches += 1
                                else:
                                    total_element_mismatches += 1
                        else:
                            if verbose:
                                print("  No elements found")

                        if data["x_match"] and data["y_match"]:
                            total_matches += 1
                        else:
                            total_mismatches += 1
                        if verbose:
                            print()

            if file_coordinate_data:
                all_coordinate_data[filename] = file_coordinate_data
        else:
            if verbose:
                print("\nNo groups containing 'EDX' found in this file.")

        # Find groups containing 'MOKE'
        moke_groups = find_moke_groups(file_path)

        if moke_groups:
            if verbose:
                print(f"\nFound {len(moke_groups)} group(s) containing 'MOKE':")
                for group_path in moke_groups:
                    print(f"  - {group_path}")

            # Process coordinate data for each MOKE group
            for group_path in moke_groups:
                total_groups_processed += 1

                # List entries for each MOKE group
                list_group_entries(file_path, group_path, verbose)

                # Process coordinates
                MOKE_sample_coordinate_data = process_moke_coordinates(
                    file_path, group_path, verbose
                )
                if MOKE_sample_coordinate_data:
                    # Merge MOKE data with existing coordinate data or create new entry
                    if filename not in all_coordinate_data:
                        all_coordinate_data[filename] = {}
                    all_coordinate_data[filename][group_path] = (
                        MOKE_sample_coordinate_data
                    )
                    total_coordinates_found += len(MOKE_sample_coordinate_data)

                    # Print coordinate analysis
                    if verbose:
                        print(f"\n--- MOKE Coordinate Analysis for {group_path} ---")
                    for sample_group, data in MOKE_sample_coordinate_data.items():
                        if verbose:
                            print(f"Sample-group: {sample_group}")
                            print(f"  x_pos_MOKE: {data['x_pos_MOKE']}")
                            print(f"  y_pos_MOKE: {data['y_pos_MOKE']}")
                            print(f"  x_pos_instrument: {data['x_pos_instrument']}")
                            print(f"  y_pos_instrument: {data['y_pos_instrument']}")
                            print(f"  x_pos_unit: {data['x_pos_unit']}")
                            print(f"  y_pos_unit: {data['y_pos_unit']}")

                        x_status = "✓" if data["x_match"] else "✗"
                        y_status = "✓" if data["y_match"] else "✗"
                        if verbose:
                            print(f"  X coordinate match: {x_status}")
                            print(f"  Y coordinate match: {y_status}")

                        # Display MOKE data
                        if data.get("moke_data", {}).get("data_available", False):
                            moke_info = data["moke_data"]
                            if verbose:
                                print(
                                    f"  Coercivity mean: {moke_info['coercivity_mean']}"
                                )
                                if moke_info.get("coercivity_unit"):
                                    print(
                                        f"  Coercivity unit: {moke_info['coercivity_unit']}"
                                    )
                        else:
                            if verbose:
                                error_msg = data.get("moke_data", {}).get(
                                    "error", "Unknown error"
                                )
                                print(f"  No MOKE data available: {error_msg}")

                        if data["x_match"] and data["y_match"]:
                            total_matches += 1
                        else:
                            total_mismatches += 1
                        if verbose:
                            print()
        else:
            if verbose:
                print("\nNo groups containing 'MOKE' found in this file.")

    # Print final summary
    if verbose:
        print(f"\n{'=' * 80}")
        print("FINAL SUMMARY")
        print(f"{'=' * 80}")
        print(f"Files processed: {len(hdf5_files)}")
        print(f"EDX groups processed: {total_groups_processed}")
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
            f"\nSUMMARY: Processed {len(hdf5_files)} file(s), found {total_coordinates_found} coordinate pairs, {total_elements_found} elements"
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
                all_coordinate_data, template_path, output_dir, include_mass_fraction
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
                    hdf5_files, script_dir, datasets_dir, output_dir, verbose
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


def find_moke_groups(file_path):
    """
    Find first-level groups containing 'MOKE' in their name within an HDF5 file.

    Args:
        file_path (str): Path to the HDF5 file

    Returns:
        list: List of first-level group paths containing 'MOKE' in their name
    """
    moke_groups = []

    try:
        with h5py.File(file_path, "r") as f:
            # Only check first-level groups (direct children of root)
            for key in f.keys():
                item = f[key]
                if isinstance(item, h5py.Group) and "MOKE" in key:
                    moke_groups.append(key)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return []

    return moke_groups


def process_moke_coordinates(file_path, moke_group, verbose=True):
    """
    Process MOKE group to extract coordinates from sub-group names and verify against instrument data.

    Args:
        file_path (str): Path to the HDF5 file
        moke_group (str): Path to the MOKE group
        verbose (bool): If True, print detailed output. If False, print minimal output.

    Returns:
        dict: Dictionary with coordinate data and verification results
    """
    MOKE_sample_coordinate_data = {}

    try:
        with h5py.File(file_path, "r") as f:
            if moke_group in f:
                group = f[moke_group]

                for key in group.keys():
                    item = group[key]
                    if isinstance(item, h5py.Group):
                        # Extract coordinates from group name
                        x_name, y_name = extract_coordinates_from_name(key)
                        if verbose:
                            print(
                                f"Processing MOKE group: {key} with coordinates ({x_name}, {y_name})"
                            )
                        if x_name is not None and y_name is not None:
                            # This is a coordinate group, rename variable to sample_key
                            sample_key = key
                            # Get instrument positions
                            if verbose:
                                print(f"  - Found coordinates: ({x_name}, {y_name})")
                            # Construct sample-group path
                            sample_group_path = f"{moke_group}/{sample_key}"
                            if verbose:
                                print(
                                    f"  - Constructed sample-group path: {sample_group_path}"
                                )
                            x_instrument, y_instrument, x_unit, y_unit = (
                                get_instrument_positions(
                                    file_path, sample_group_path, verbose
                                )
                            )

                            # Get MOKE data from results sub-group
                            moke_data = get_moke_data(
                                file_path, sample_group_path, verbose
                            )

                            # Store data
                            MOKE_sample_coordinate_data[sample_key] = {
                                "x_pos_MOKE": x_name,
                                "y_pos_MOKE": y_name,
                                "x_pos_instrument": x_instrument,
                                "y_pos_instrument": y_instrument,
                                "x_pos_unit": x_unit,
                                "y_pos_unit": y_unit,
                                "moke_data": moke_data,
                                "x_match": x_instrument is not None
                                and abs(x_name - x_instrument) < 1e-6,
                                "y_match": y_instrument is not None
                                and abs(y_name - y_instrument) < 1e-6,
                            }

    except Exception as e:
        if verbose:
            print(f"Error processing MOKE coordinates for {moke_group}: {e}")

    return MOKE_sample_coordinate_data


def get_moke_data(file_path, sample_group_path, verbose=True):
    """
    Extract MOKE data from the 'results' sub-group, specifically the mean value of coercivity_m0.

    Args:
        file_path (str): Path to the HDF5 file
        sample_group_path (str): Path to the sample group
        verbose (bool): If True, print detailed output. If False, print minimal output.

    Returns:
        dict: Dictionary with MOKE data and verification results
    """
    moke_data = {}

    try:
        with h5py.File(file_path, "r") as f:
            results_path = f"{sample_group_path}/results"

            if results_path in f:
                if verbose:
                    print(f"Reading MOKE results data from: {results_path}")
                results_group = f[results_path]

                # Look for coercivity_m0 dataset
                if "coercivity_m0" in results_group:
                    coercivity_m0_group = results_group["coercivity_m0"]
                    if verbose:
                        print("Found 'coercivity_m0' in results group.")

                    # Get mean value if available
                    if "mean" in coercivity_m0_group:
                        mean_dataset = coercivity_m0_group["mean"]
                        mean_value = float(mean_dataset[()])
                        if verbose:
                            print(f"  Coercivity mean value: {mean_value}")

                        # Extract unit attribute if it exists
                        unit = None
                        if "units" in mean_dataset.attrs:
                            unit = mean_dataset.attrs["units"]
                            if isinstance(unit, bytes):
                                unit = unit.decode("utf-8")
                            if verbose:
                                print(f"  Coercivity unit: {unit}")

                        moke_data = {
                            "coercivity_mean": mean_value,
                            "coercivity_unit": unit,
                            "data_available": True,
                        }
                    else:
                        if verbose:
                            print("  No 'mean' dataset found in coercivity_m0")
                        moke_data = {
                            "data_available": False,
                            "error": "No mean value found",
                        }
                else:
                    if verbose:
                        print("  No 'coercivity_m0' found in results")
                    moke_data = {
                        "data_available": False,
                        "error": "No coercivity_m0 found",
                    }
            else:
                if verbose:
                    print(f"No 'results' group found in {sample_group_path}")
                moke_data = {"data_available": False, "error": "No results group found"}

    except Exception as e:
        if verbose:
            print(f"Error reading MOKE data from {sample_group_path}: {e}")
        moke_data = {"data_available": False, "error": str(e)}

    return moke_data


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
                if "EDX" in group_name:
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

                elif "MOKE" in group_name:
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
    # Default: process only NdCeFeB_2-5.hdf5
    main()

    # Alternative usage examples:
    # process_single_file()                                        # Process default file (NdCeFeB_2-5.hdf5), no mass fractions, no zip
    # process_single_file("another_file.hdf5")                     # Process a specific file, no mass fractions, no zip
    # process_single_file("file.hdf5", include_mass_fraction=True) # Include mass fractions in output, no zip
    # process_single_file("file.hdf5", create_zip=True)            # Process file and create zip file
    # process_all_files()                                          # Process all files in datasets directory, no mass fractions, no zip
    # process_all_files(include_mass_fraction=True)                # Process all files with mass fractions, no zip
    # process_all_files(create_zip=True)                           # Process all files and create zip file
    # main(single_file="specific_file.hdf5")                       # Direct call with specific file, no mass fractions, no zip
    # main(single_file="", include_mass_fraction=True)             # Direct call to process all files with mass fractions, no zip
    # main(single_file="", create_zip=True)                        # Direct call to process all files and create zip file

    # Note: The NEEL_template.archive.yaml has mass_fraction commented out by default.
    # When include_mass_fraction=True, the script will automatically uncomment those lines.
    # When include_mass_fraction=False (default), mass_fraction remains commented out.
    #
    # Zip file creation is disabled by default (create_zip=False).
    # When create_zip=True, the script will create a zip file containing all HDF5 and generated YAML files.
