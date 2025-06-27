#!/usr/bin/env python3
"""
Batch processing script for UU datasets.
Runs run_convert_UU_data_to_h5_and_creat_yaml.py for multiple datasets.
"""

import subprocess
import sys
import os
from typing import List


def main():
    """
    Process multiple UU datasets by running the conversion script for each one.
    """

    # List of datasets to process
    datasets_to_process: List[str] = [
        "Fe12Co4N2",
        "Fe16N2",
        "Fe3Y",
        "Fe2.78Nb0.22Y",
        "Fe2.78Ta0.22Y",
        "Fe2.78V0.22Y",
        "Fe2.78Zn0.22Y",
    ]

    # Configuration
    datasets_dir = "./datasets"
    conversion_script = "run_convert_UU_data_to_h5_and_creat_yaml.py"

    # Check if the conversion script exists
    if not os.path.exists(conversion_script):
        print(f"Error: Conversion script '{conversion_script}' not found.")
        sys.exit(1)

    # Check if datasets directory exists
    if not os.path.exists(datasets_dir):
        print(f"Error: Datasets directory '{datasets_dir}' not found.")
        sys.exit(1)

    print(f"Starting batch processing of {len(datasets_to_process)} datasets...")
    print(f"Datasets directory: {os.path.abspath(datasets_dir)}")
    print(f"Conversion script: {conversion_script}")
    print("-" * 60)

    successful_conversions = []
    failed_conversions = []

    for i, dataset in enumerate(datasets_to_process, 1):
        print(f"\n[{i}/{len(datasets_to_process)}] Processing dataset: {dataset}")
        print("-" * 40)

        # Check if dataset directory exists
        dataset_path = os.path.join(datasets_dir, dataset)
        if not os.path.isdir(dataset_path):
            print(f"Warning: Dataset directory '{dataset_path}' not found. Skipping...")
            failed_conversions.append(f"{dataset} (directory not found)")
            continue

        # Run the conversion script for this dataset
        cmd = [
            sys.executable,
            conversion_script,
            "--source_dir",
            dataset,
            "--datasets_dir",
            datasets_dir,
        ]

        try:
            print(f"Running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=False)
            print(f"✓ Successfully processed {dataset}")
            successful_conversions.append(dataset)

        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to process {dataset}: {e}")
            failed_conversions.append(f"{dataset} (processing error)")

        except Exception as e:
            print(f"✗ Unexpected error processing {dataset}: {e}")
            failed_conversions.append(f"{dataset} (unexpected error)")

    # Print summary
    print("\n" + "=" * 60)
    print("BATCH PROCESSING SUMMARY")
    print("=" * 60)
    print(f"Total datasets: {len(datasets_to_process)}")
    print(f"Successful: {len(successful_conversions)}")
    print(f"Failed: {len(failed_conversions)}")

    if successful_conversions:
        print("\n✓ Successfully processed datasets:")
        for dataset in successful_conversions:
            print(f"  - {dataset}")

    if failed_conversions:
        print("\n✗ Failed to process datasets:")
        for dataset in failed_conversions:
            print(f"  - {dataset}")

    # Exit with appropriate code
    if failed_conversions:
        print("\nSome datasets failed to process. Check the errors above.")
        sys.exit(1)
    else:
        print("\nAll datasets processed successfully!")
        sys.exit(0)


if __name__ == "__main__":
    main()
