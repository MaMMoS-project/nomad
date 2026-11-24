import os
import tempfile
import unittest

from packages.data_tests_functions import interpolate_thermal_property, read_thermal_csv

SAMPLE_CSV = """# mammos csv v3
#----------------------------------------
# Temperature-dependent quantities computed with UppASD
#----------------------------------------
ThermodynamicTemperature,SpontaneousMagneticPolarization
,Computed as derivative dE/dT.
https://w3id.org/emmo#EMMO_affe07e4_e9bc_4852_86c6_69e26182a17f,https://w3id.org/emmo/domain/magnetic-materials#EMMO_db6f7b13-1f1d-584f-9d73-47939b86a7cd
K,T
T,Js
0.0,1.0
300.0,0.9
600.0,0.6
"""


class TestReadThermalCsv(unittest.TestCase):
    def _write(self, tmp):
        path = os.path.join(tmp, "thermal.csv")
        with open(path, "w") as f:
            f.write(SAMPLE_CSV)
        return path

    def test_columns_and_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            df, meta = read_thermal_csv(self._write(tmp))
        self.assertEqual(meta["columns"], ["T", "Js"])
        self.assertEqual(list(df.columns), ["T", "Js"])
        self.assertEqual(len(df), 3)

    def test_ontology_meta(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, meta = read_thermal_csv(self._write(tmp))
        self.assertEqual(meta["ontology_labels"]["T"], "ThermodynamicTemperature")
        self.assertEqual(meta["ontology_labels"]["Js"], "SpontaneousMagneticPolarization")
        self.assertTrue(meta["ontology_iris"]["Js"].startswith("https://w3id.org/emmo"))
        self.assertEqual(meta["units"]["T"], "K")
        self.assertEqual(meta["units"]["Js"], "T")
        self.assertEqual(meta["descriptions"]["Js"], "Computed as derivative dE/dT.")

    def test_comment_strips_separator_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, meta = read_thermal_csv(self._write(tmp))
        self.assertEqual(
            meta["comment"],
            "mammos csv v3\nTemperature-dependent quantities computed with UppASD",
        )

    def test_interpolate(self):
        with tempfile.TemporaryDirectory() as tmp:
            df, _ = read_thermal_csv(self._write(tmp))
        self.assertAlmostEqual(interpolate_thermal_property(df, "Js", 300.0), 0.9)
        self.assertAlmostEqual(interpolate_thermal_property(df, "Js", 450.0), 0.75)


if __name__ == "__main__":
    unittest.main()
