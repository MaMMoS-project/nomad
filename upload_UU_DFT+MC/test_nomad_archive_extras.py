import os
import tempfile
import unittest

from packages.nomad_archive_extras import (
    chemical_formula_from_name,
    elements_from_formula,
    h5_results_map,
    mammos_entity_records,
    records_for_computed,
    structure_run_section,
)

MINIMAL_CIF = """data_test
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0.0 0.0 0.0
"""

SAMPLE_PROPS = {
    "Js": {
        "value": 2.003,
        "unit": "T",
        "ontology_label": "SpontaneousMagneticPolarization",
        "ontology_iri": "https://w3id.org/emmo/domain/magnetic-materials#EMMO_db6f7b13-1f1d-584f-9d73-47939b86a7cd",
        "description": "",
    },
    "Tc": {
        "value": 1410.0,
        "unit": "K",
        "ontology_label": "CurieTemperature",
        "ontology_iri": "https://w3id.org/emmo#EMMO_6b5af5a8_a2d8_4353_a1d6_54c9f778343d",
        "description": "Tc from peak in specific heat",
    },
}


class TestStructureRunSection(unittest.TestCase):
    def test_cif_parsed_to_run_section(self):
        with tempfile.TemporaryDirectory() as tmp:
            cif = os.path.join(tmp, "structure.cif")
            with open(cif, "w") as f:
                f.write(MINIMAL_CIF)
            run = structure_run_section(cif)
        self.assertEqual(run["program"]["name"], "Structure File Importer")
        atoms = run["system"][0]["atoms"]
        self.assertEqual(atoms["labels"], ["Fe"])
        self.assertEqual(atoms["species"], [26])
        # lattice converted to meters
        self.assertAlmostEqual(atoms["lattice_vectors"][0][0], 5.0e-10)
        self.assertAlmostEqual(atoms["lattice_vectors"][2][2], 6.0e-10)
        self.assertEqual(atoms["periodic"], [True, True, True])


class TestChemicalFormulaFromName(unittest.TestCase):
    def test_plain_and_comment_suffix(self):
        self.assertEqual(chemical_formula_from_name("Fe12Co4N2"), "Fe12Co4N2")
        self.assertEqual(chemical_formula_from_name("Nd1.5Ce0.5Fe14B"), "Nd1.5Ce0.5Fe14B")
        self.assertEqual(chemical_formula_from_name("Fe16N2_exp_str"), "Fe16N2")
        self.assertEqual(chemical_formula_from_name("Nd2Fe14B_method1"), "Nd2Fe14B")


class TestH5ResultsMap(unittest.TestCase):
    def test_mammos_collection_layout(self):
        import h5py

        with tempfile.TemporaryDirectory() as tmp:
            h5 = os.path.join(tmp, "x.h5")
            with h5py.File(h5, "w") as f:
                for k in ["Js", "Ku", "Tc", "Ms"]:
                    f.create_dataset(f"intrinsic_properties/{k}", data=1.0)
            m = h5_results_map(h5)
        self.assertEqual(
            m,
            {
                "Js": "/intrinsic_properties/Js",
                "Ku": "/intrinsic_properties/Ku",
                "Tc": "/intrinsic_properties/Tc",
                "Ms": "/intrinsic_properties/Ms",
            },
        )

    def test_old_layout_aliases(self):
        import h5py

        with tempfile.TemporaryDirectory() as tmp:
            h5 = os.path.join(tmp, "x.h5")
            with h5py.File(h5, "w") as f:
                for k in ["Js_0", "Js_300", "K1_300"]:
                    f.create_dataset(f"results/{k}", data=1.0)
            m = h5_results_map(h5)
        self.assertEqual(m["Js"], "/results/Js_0")
        self.assertEqual(m["Js_300"], "/results/Js_300")
        self.assertEqual(m["Ku"], "/results/K1_300")
        self.assertNotIn("Tc", m)

    def test_missing_file(self):
        self.assertEqual(h5_results_map("/nonexistent.h5"), {})


class TestElementsFromFormula(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(elements_from_formula("ZnFe2Ni9"), ["Fe", "Ni", "Zn"])
        self.assertEqual(elements_from_formula("Al2Fe6C"), ["Al", "C", "Fe"])

    def test_fractional(self):
        self.assertEqual(
            elements_from_formula("Nd1.5Ce0.5Fe14B"), ["B", "Ce", "Fe", "Nd"]
        )

    def test_unparseable_raises(self):
        with self.assertRaises(ValueError):
            elements_from_formula("Fe3-xZnxY(x=0.22)")


class TestMammosEntityRecords(unittest.TestCase):
    def test_records_from_intrinsic_properties(self):
        records = mammos_entity_records(SAMPLE_PROPS)
        by_name = {r["name"]: r for r in records}
        self.assertEqual(set(by_name), {"Js", "Tc"})
        tc = by_name["Tc"]
        # mammos-entity names are used as-is for name and value_key
        self.assertEqual(tc["value_key"], "Tc")
        self.assertEqual(tc["ontology_label"], "CurieTemperature")
        self.assertEqual(
            tc["ontology_iri"],
            "https://w3id.org/emmo#EMMO_6b5af5a8_a2d8_4353_a1d6_54c9f778343d",
        )
        self.assertEqual(tc["unit"], "K")
        self.assertEqual(by_name["Js"]["value_key"], "Js")

    def test_records_for_computed(self):
        records = records_for_computed(2.0, 1.8, 3.2)  # ku in MJ/m^3
        by_key = {r["value_key"]: r for r in records}
        self.assertEqual(set(by_key), {"Js", "Js_300", "Ku"})
        self.assertEqual(
            by_key["Js"]["ontology_label"], "SpontaneousMagneticPolarization"
        )
        self.assertEqual(by_key["Ku"]["ontology_label"], "UniaxialAnisotropyConstant")
        self.assertEqual(by_key["Ku"]["unit"], "MJ/m3")
        for r in records:
            self.assertTrue(r["ontology_iri"].startswith("https://w3id.org/emmo"))


if __name__ == "__main__":
    unittest.main()
