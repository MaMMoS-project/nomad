import h5py
import os


def find_edx_groups(file_path):
    """
    Find all groups containing 'EDX' in their name within an HDF5 file.

    Args:
        file_path (str): Path to the HDF5 file

    Returns:
        list: List of group paths containing 'EDX' in their name
    """
    edx_groups = []

    try:
        with h5py.File(file_path, "r") as f:

            def visit_item(name, obj):
                if isinstance(obj, h5py.Group) and "EDX" in name:
                    edx_groups.append(name)

            f.visititems(visit_item)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return []

    return edx_groups


def list_group_entries(file_path, group_path):
    """
    List all entries (groups and datasets) within a specific group.

    Args:
        file_path (str): Path to the HDF5 file
        group_path (str): Path to the group within the HDF5 file
    """
    try:
        with h5py.File(file_path, "r") as f:
            if group_path in f:
                group = f[group_path]
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
            else:
                print(f"Group '{group_path}' not found in file.")
    except Exception as e:
        print(f"Error accessing group {group_path} in file {file_path}: {e}")


def main():
    """
    Main function to process HDF5 files in the datasets subfolder.
    """
    # Define the datasets directory relative to the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    datasets_dir = os.path.join(script_dir, "datasets")

    if not os.path.exists(datasets_dir):
        print(f"Datasets directory not found: {datasets_dir}")
        return

    # Get all HDF5 files in the datasets directory
    hdf5_files = [f for f in os.listdir(datasets_dir) if f.endswith(".hdf5")]

    if not hdf5_files:
        print("No HDF5 files found in the datasets directory.")
        return

    print(f"Found {len(hdf5_files)} HDF5 file(s):")
    for i, filename in enumerate(hdf5_files, 1):
        print(f"  {i}. {filename}")

    # Process each HDF5 file
    for filename in hdf5_files:
        file_path = os.path.join(datasets_dir, filename)
        print(f"\n{'=' * 60}")
        print(f"Processing file: {filename}")
        print(f"{'=' * 60}")

        # Find groups containing 'EDX'
        edx_groups = find_edx_groups(file_path)

        if edx_groups:
            print(f"\nFound {len(edx_groups)} group(s) containing 'EDX':")
            for group_path in edx_groups:
                print(f"  - {group_path}")

            # List entries for each EDX group
            for group_path in edx_groups:
                list_group_entries(file_path, group_path)
        else:
            print("\nNo groups containing 'EDX' found in this file.")


if __name__ == "__main__":
    main()
