#!/usr/bin/env python3
"""
Batch processing script for UU datasets.
Runs run_convert_UU_data_to_h5_and_creat_yaml.py for multiple datasets.
"""

import subprocess
import sys
import os
import argparse
from typing import List

from run_convert_UU_data_to_h5_and_creat_yaml import resolve_datasets_dir

# Small default set for quick test rebuilds (--quick): covers a plain folder
# name, a '_suffix' name and a multi-MC dataset. Run ALL datasets with no
# flag (or an explicit --only list) when preparing real uploads.
QUICK_DATASETS = ["Al2Fe6C", "Fe16N2_exp_str", "Co3Mn2Ge"]


def main():
    """
    Process multiple UU datasets by running the conversion script for each one.
    By default all subfolders of the datasets directory are processed
    (./uppsala-data/data if present, else ./datasets).
    """

    parser = argparse.ArgumentParser(
        description="Batch-process datasets: runs run_convert_UU_data_to_h5_and_creat_yaml.py "
        "once per subfolder of the datasets directory (one .h5 + .archive.yaml + upload "
        "zip each). Use the runner directly for a single dataset.",
        epilog="""examples:
  %(prog)s
      process ALL subfolders of the default datasets directory
  %(prog)s --datasets_dir ./uppsala-data/data
      process all subfolders of an explicit datasets directory
  %(prog)s --only Fe12Co4N2,Co2Fe2H4
      process only the listed datasets (comma-separated)
  %(prog)s --quick
      quick test rebuild: only the small QUICK_DATASETS set (currently
      Al2Fe6C, Fe16N2_exp_str, Co3Mn2Ge — plain name, '_suffix' name,
      multi-MC), skipping the slow per-dataset mammos-parser validation of
      all 23 datasets. Run without flags to rebuild ALL datasets.
Note: each dataset is processed with the runner's defaults (--cif off).
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick test rebuild: process only the small QUICK_DATASETS set "
        f"({', '.join(QUICK_DATASETS)}) instead of all subfolders; "
        "--only overrides this when both are given",
    )
    parser.add_argument(
        "--datasets_dir",
        nargs="?",
        default=None,
        help="Directory whose subfolders are the datasets (default: ./uppsala-data/data if present, else ./datasets)",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Comma-separated subset of dataset names to process (default: all subfolders)",
    )
    args = parser.parse_args()

    # Configuration
    datasets_dir = resolve_datasets_dir(args.datasets_dir)
    conversion_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "run_convert_UU_data_to_h5_and_creat_yaml.py",
    )

    # List of datasets to process: every subfolder of the datasets directory
    datasets_to_process: List[str] = sorted(
        d
        for d in os.listdir(datasets_dir)
        if os.path.isdir(os.path.join(datasets_dir, d))
    )
    if args.only:
        only = {name.strip() for name in args.only.split(",") if name.strip()}
        datasets_to_process = [d for d in datasets_to_process if d in only]
    elif args.quick:
        quick = [d for d in QUICK_DATASETS if d in datasets_to_process]
        missing = [d for d in QUICK_DATASETS if d not in datasets_to_process]
        if missing:
            print(f"Warning: QUICK_DATASETS not present in {datasets_dir}: {missing}")
        datasets_to_process = quick

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
