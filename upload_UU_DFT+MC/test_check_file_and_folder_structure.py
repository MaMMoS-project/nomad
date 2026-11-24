import os
import unittest
from packages.check_file_and_folder_structure import (
    check_readme_exists,
    check_subfolders_exist,
    check_structure_cif_file_exists,
    # check_momfile_file_exists,
    check_structure,
)

# Legacy old-layout (GS/Jij/MC) datasets; the test passes bare folder names, so
# it expects to run with those folders in the CWD (e.g. from ./datasets).
LEGACY_DATASETS = [
    "Co2Fe16Y6",
    "Fe3-xTaxY(x=0.22)",
    "Fe3-xVxY(x=0.22)",
    "Fe3-xZnxY(x=0.22)",
]


@unittest.skipUnless(
    all(os.path.isdir(d) for d in LEGACY_DATASETS),
    "legacy old-layout dataset folders not present in CWD (run from ./datasets); "
    "new uppsala-data layout is covered by TestUppsaladataStructure",
)
class TestCheckFileAndFolderStructure(unittest.TestCase):
    def setUp(self):
        self.datasets = LEGACY_DATASETS

    def test_check_readme_exists(self):
        all_results_true = []
        print("\n##############\nRunning test_check_readme_exists")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                result = check_readme_exists(dataset, verbose="debug")
                print(f"README check for {dataset}: {result}")
        if all(all_results_true):
            print("All README checks passed for all datasets.")
        else:
            print("Some README checks failed for some datasets.")
            raise AssertionError("Some README checks failed for some datasets.")

    def test_check_subfolders_exist(self):
        all_results_true = []
        print("\n##############\nRunning test_check_subfolders_exist")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                mandatory_exist, optional_exist = check_subfolders_exist(
                    dataset, check_optional_subf=True, verbose="debug"
                )
                print(f"Mandatory subfolders check for {dataset}: {mandatory_exist}")
                print(f"Optional subfolders check for {dataset}: {optional_exist}")
        if all(all_results_true):
            print("All subfolder checks passed for all datasets.")
        else:
            print("Some subfolder checks failed for some datasets.")
            raise AssertionError("Some subfolder checks failed for some datasets.")

    def test_check_structure_cif_file_exists(self):
        all_results_true = []
        print("\n##############\nRunning test_check_structure_cif_file_exists")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                result = check_structure_cif_file_exists(dataset, verbose="debug")
                print(f"Structure CIF file check for {dataset}: {result}")
                all_results_true.append(result)
        if all(all_results_true):
            print("All structure CIF file checks passed for all datasets.")
        else:
            print("Some structure CIF file checks failed for some datasets.")
            raise AssertionError(
                "Some structure CIF file checks failed for some datasets."
            )

    # def test_check_structure_momfile_file_exists(self):
    #     all_results_true = []
    #     print("\n##############\nRunning test_check_structure_momfile_file_exists")
    #     for dataset in self.datasets:
    #         with self.subTest(dataset=dataset):
    #             result = check_momfile_file_exists(
    #                 dataset, verbose='debug')
    #             print(f"'momfile' check for {dataset}: {result}")
    #             all_results_true.append(result)
    #     if all(all_results_true):
    #         print("All 'momfile' checks passed for all datasets.")
    #     else:
    #         print("Some 'momfile' checks failed for some datasets.")
    #         raise AssertionError(
    #             "Some 'momfile' checks failed for some datasets.")

    def test_check_structure(self):
        print("\n##############\nRunning test_check_structure")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                result = check_structure(
                    dataset,
                    check_optional_subf=True,
                    check_README=True,
                    check_subfolders=True,
                    check_structure_cif=True,
                    check_out_last_files=True,
                    verbose="debug",
                )
                self.assertTrue(result, f"Check failed for dataset: {dataset}")


class TestUppsaladataStructure(unittest.TestCase):
    """Tests for the new uppsala-data layout (RSPt/ + UppASD/)."""

    def _make_minimal_dataset(self, root):
        import os

        for f in ["metadata.yaml", "intrinsic_properties.yaml", "structure.cif"]:
            open(os.path.join(root, f), "w").close()
        for f in ["atomdens", "kmap", "spts", "symcof", "symt.inp"]:
            os.makedirs(os.path.join(root, "RSPt/common_input"), exist_ok=True)
            open(os.path.join(root, "RSPt/common_input", f), "w").close()
        for gs in ["gs_x", "gs_z"]:
            os.makedirs(os.path.join(root, "RSPt", gs), exist_ok=True)
            for f in ["data", "out_last", "out_MF"]:
                open(os.path.join(root, "RSPt", gs, f), "w").close()
        os.makedirs(os.path.join(root, "RSPt/Jij"), exist_ok=True)
        for f in ["data", "out_last"]:
            open(os.path.join(root, "RSPt/Jij", f), "w").close()
        os.makedirs(os.path.join(root, "UppASD/MC_1"), exist_ok=True)
        for f in ["jfile", "momfile", "posfile", "inpsd.dat", "thermal.csv", "thermal.dat"]:
            open(os.path.join(root, "UppASD/MC_1", f), "w").close()

    def test_is_uppsala_structure(self):
        import tempfile
        from packages.check_file_and_folder_structure import is_uppsala_structure

        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(is_uppsala_structure(tmp))
            self._make_minimal_dataset(tmp)
            self.assertTrue(is_uppsala_structure(tmp))

    def test_check_uppsala_structure_ok(self):
        import tempfile
        from packages.check_file_and_folder_structure import check_uppsala_structure

        with tempfile.TemporaryDirectory() as tmp:
            self._make_minimal_dataset(tmp)
            self.assertTrue(check_uppsala_structure(tmp, verbose="off"))

    def test_check_uppsala_structure_missing_file(self):
        import os
        import tempfile
        from packages.check_file_and_folder_structure import check_uppsala_structure

        with tempfile.TemporaryDirectory() as tmp:
            self._make_minimal_dataset(tmp)
            os.remove(os.path.join(tmp, "UppASD/MC_1/thermal.csv"))
            self.assertFalse(check_uppsala_structure(tmp, verbose="off"))


if __name__ == "__main__":
    unittest.main()
