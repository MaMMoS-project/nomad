import unittest
from check_file_and_folder_structure import (
    check_readme_exists,
    check_subfolders_exist,
    check_structure_cif_file_exists,
    # check_momfile_file_exists,
    check_structure
)


class TestCheckFileAndFolderStructure(unittest.TestCase):

    def setUp(self):
        self.datasets = ['Co2Fe16Y6', 'Fe3-xTaxY(x=0.22)', 'Fe3-xVxY(x=0.22)', 'Fe3-xZnxY(x=0.22)'] 

    def test_check_readme_exists(self):
        all_results_true = []
        print("\n##############\nRunning test_check_readme_exists")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                result = check_readme_exists(dataset, verbose='debug')
                print(f"README check for {dataset}: {result}")
        if all(all_results_true):
            print("All README checks passed for all datasets.")
        else:
            print("Some README checks failed for some datasets.")
            raise AssertionError(
                "Some README checks failed for some datasets.")

    def test_check_subfolders_exist(self):
        all_results_true = []
        print("\n##############\nRunning test_check_subfolders_exist")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                mandatory_exist, optional_exist = check_subfolders_exist(
                    dataset, check_optional_subf=True, verbose='debug')
                print(
                    f"Mandatory subfolders check for {dataset}: {mandatory_exist}")
                print(
                    f"Optional subfolders check for {dataset}: {optional_exist}")
        if all(all_results_true):
            print("All subfolder checks passed for all datasets.")
        else:
            print("Some subfolder checks failed for some datasets.")
            raise AssertionError(
                "Some subfolder checks failed for some datasets.")

    def test_check_structure_cif_file_exists(self):
        all_results_true = []
        print("\n##############\nRunning test_check_structure_cif_file_exists")
        for dataset in self.datasets:
            with self.subTest(dataset=dataset):
                result = check_structure_cif_file_exists(
                    dataset, verbose='debug')
                print(f"Structure CIF file check for {dataset}: {result}")
                all_results_true.append(result)
        if all(all_results_true):
            print("All structure CIF file checks passed for all datasets.")
        else:
            print("Some structure CIF file checks failed for some datasets.")
            raise AssertionError(
                "Some structure CIF file checks failed for some datasets.")

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
                result = check_structure(dataset, check_optional_subf=True, check_README=True,
                                         check_subfolders=True, check_structure_cif=True, check_out_last_files=True, verbose='debug')
                self.assertTrue(result, f"Check failed for dataset: {dataset}")


if __name__ == '__main__':
    unittest.main()
