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

ANGSTROM_TO_METER = 1e-10


def ontology_description(ontology_label="", ontology_iri="", extra=""):
    """Build a quantity description string from ontology metadata, e.g.
    'Ku from force theorem | ontology: UniaxialAnisotropyConstant | IRI: [...]'."""
    parts = []
    if extra:
        parts.append(extra)
    if ontology_label:
        parts.append(f"ontology: {ontology_label}")
    if ontology_iri:
        parts.append(f"IRI: [{ontology_iri}]({ontology_iri})")
    return " | ".join(parts)


def chemical_formula_from_name(name):
    """Best-effort chemical formula from an uppsala-data folder name.

    The folder-name convention (uppsala-data README) is a reduced chemical
    formula with an optional '_<description>' suffix, e.g. 'Nd2Fe14B',
    'Nd1.5Ce0.5Fe14B' or 'Fe16N2_exp_str' -> 'Fe16N2'.
    """
    return name.split("_", 1)[0]


def mammos_yaml_description(file_path):
    """Top-level 'description' of a mammos yaml v2 file (e.g.
    intrinsic_properties.yaml), or None if empty/absent."""
    import yaml

    with open(file_path, "rt") as f:
        content = yaml.safe_load(f)
    return (content.get("description") or "").strip() or None


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
    unit, description. The mammos-entity name is used as-is for both `name`
    and `value_key` (the NOMAD template quantity has the same name).
    """
    records = []
    for name, entry in props.items():
        records.append(
            {
                "name": name,
                "value_key": name,
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


def records_for_computed(js_0, js_300, ku_mj):
    """Build ``mammos_entities`` records for the old-layout pipeline, where
    Js/Js_300/Ku are computed from raw output files. Ontology metadata
    is taken from mammos_entity (MaMMoS/EMMO magnetic-materials ontology).
    ku_mj is the anisotropy constant in MJ/m^3 (template unit of Ku).
    """
    import mammos_entity as me

    records = []
    for name, entity, unit in [
        ("Js", me.Js(js_0 if js_0 is not None else 0.0, "T"), "T"),
        ("Js_300", me.Js(js_300 if js_300 is not None else 0.0, "T"), "T"),
        ("Ku", me.Ku(ku_mj if ku_mj is not None else 0.0, "MJ/m3"), "MJ/m3"),
    ]:
        info = _ontology_info(entity)
        records.append(
            {
                "name": name,
                "value_key": name,
                "source_type": "Entity",
                "ontology_label": info["ontology_label"],
                "ontology_iri": info["ontology_iri"],
                "unit": unit,
                "description": getattr(entity, "description", "") or "",
            }
        )
    return records
