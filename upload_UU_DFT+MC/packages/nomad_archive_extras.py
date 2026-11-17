"""
Helpers to enrich the NOMAD .archive.yaml with ontology metadata and the
crystal structure, following the approach of MaMMoS-project/to-nomad
(branch feature/subcollection-subsections):

- structure_run_section(): parse a structure.cif locally with ASE and embed
  it as an inline top-level ``run`` section in the archive. Together with
  ``metadata: {domain: dft}`` the server-side OptimadeNormalizer then
  populates the 3-D structure preview in the NOMAD entry overview.
  (NOMAD never reads the .cif at upload time; the .cif is shipped in the
  zip as raw data/provenance only. The ElnWithStructureFile base section is
  deliberately NOT used - it fails on nomad-lab 2.0.0.)
- mammos_entity_records(): build the repeating ``mammos_entities``
  subsection records (name, value_key, source_type, ontology_label,
  ontology_iri, unit, description) for the archive, with the ontology
  metadata from the MaMMoS/EMMO magnetic-materials ontology as stored in
  uppsala-data intrinsic_properties.yaml (or from mammos_entity itself).
"""

import os

ANGSTROM_TO_METER = 1e-10

# Mapping from intrinsic_properties.yaml keys to the NOMAD template quantities
ME_VALUE_KEY = {"Js": "Js_0", "Ms": "Ms", "Ku": "K1_300", "Tc": "Tc"}


def structure_run_section(cif_path):
    """Parse a structure file (e.g. .cif) with ASE into a NOMAD run-section dict.

    Returns a dict with 'program' and 'system' keys for the archive-level
    'run' list. Lattice vectors and positions are converted to the runschema
    base unit (meter).
    """
    from ase.io import read

    structure = read(str(cif_path))
    return {
        "program": {"name": "Structure File Importer"},
        "system": [
            {
                "atoms": {
                    "lattice_vectors": (
                        structure.get_cell()[:] * ANGSTROM_TO_METER
                    ).tolist(),
                    "labels": list(structure.get_chemical_symbols()),
                    "positions": (
                        structure.get_positions() * ANGSTROM_TO_METER
                    ).tolist(),
                    "periodic": [bool(p) for p in structure.get_pbc()],
                    "species": [int(z) for z in structure.get_atomic_numbers()],
                }
            }
        ],
    }


def mammos_entity_records(props):
    """Build ``mammos_entities`` subsection records from an
    intrinsic_properties.yaml dict (see data_tests_functions.read_intrinsic_properties).

    Each record: name, value_key, source_type, ontology_label, ontology_iri,
    unit, description. value_key maps the entity to the NOMAD template
    quantity (ME_VALUE_KEY); entities without a mapping are skipped.
    """
    records = []
    for name, entry in props.items():
        value_key = ME_VALUE_KEY.get(name)
        if value_key is None:
            continue
        records.append(
            {
                "name": name,
                "value_key": value_key,
                "source_type": "Entity",
                "ontology_label": entry.get("ontology_label") or "",
                "ontology_iri": entry.get("ontology_iri") or "",
                "unit": entry.get("unit") or "",
                "description": entry.get("description") or "",
            }
        )
    return records


def _ontology_info(entity):
    """Extract ontology_label / ontology_iri / unit from a mammos_entity Entity."""
    return {
        "ontology_label": getattr(entity, "ontology_label", None) or "",
        "ontology_iri": getattr(entity, "ontology_iri", None) or "",
        "unit": str(getattr(entity, "unit", "") or ""),
    }


def records_for_computed(js_0, js_300, k1_300):
    """Build ``mammos_entities`` records for the old-layout pipeline, where
    Js_0/Js_300/K1_300 are computed from raw output files. Ontology metadata
    is taken from mammos_entity (MaMMoS/EMMO magnetic-materials ontology).
    """
    import mammos_entity as me

    records = []
    for name, value_key, entity in [
        ("Js_0", "Js_0", me.Js(js_0 if js_0 is not None else 0.0, "T")),
        ("Js_300", "Js_300", me.Js(js_300 if js_300 is not None else 0.0, "T")),
        ("K1_300", "K1_300", me.Ku(k1_300 if k1_300 is not None else 0.0, "J/m3")),
    ]:
        info = _ontology_info(entity)
        records.append(
            {
                "name": name,
                "value_key": value_key,
                "source_type": "Entity",
                "ontology_label": info["ontology_label"],
                "ontology_iri": info["ontology_iri"],
                "unit": info["unit"],
                "description": getattr(entity, "description", "") or "",
            }
        )
    return records
