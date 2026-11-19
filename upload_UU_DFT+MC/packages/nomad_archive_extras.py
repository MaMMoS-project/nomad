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
- dataset_collection() / raw_data_collection(): build one mammos-entity
  EntityCollection per uppsala-data dataset (sub-collections
  'intrinsic_properties', 'thermal_MC_*' and 'raw_data') for writing with
  EntityCollection.to_hdf5().
"""

import os

import numpy as np

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


# Template quantity -> candidate h5 dataset paths for HDF5Normalizer linking.
# New uppsala-data layout (mammos collection written with to_hdf5) stores the
# intrinsic properties under /intrinsic_properties/; the old hand-written
# layout stored them in /results/ (Js/Ku/Tc/Ms or Js_0/K1_300/Js_300).
H5_RESULT_ALIASES = {
    "Js": ["/intrinsic_properties/Js", "/results/Js", "/results/Js_0"],
    "Js_300": ["/results/Js_300"],
    "Ku": ["/intrinsic_properties/Ku", "/results/Ku", "/results/K1_300"],
    "Tc": ["/intrinsic_properties/Tc", "/results/Tc"],
    "Ms": ["/intrinsic_properties/Ms", "/results/Ms"],
}


def h5_results_map(h5_file):
    """Map template quantities to h5 dataset paths for HDF5Normalizer linking.

    Returns {quantity_name: h5_path} for every template quantity whose alias
    path exists in the h5 file (checked with h5py).
    """
    import h5py

    mapping = {}
    try:
        with h5py.File(h5_file, "r") as f:
            for qname, paths in H5_RESULT_ALIASES.items():
                for p in paths:
                    if p in f:
                        mapping[qname] = p
                        break
    except OSError:
        pass
    return mapping


def raw_data_collection(source_dir, exclude_DOSCAR=False):
    """Mirror a dataset folder tree into a nested mammos EntityCollection.

    Plain values are stored as plain h5 datasets by EntityCollection.to_hdf5
    (arrays as numeric datasets, texts as string datasets). .csv files are
    parsed as numeric arrays (mammos csv v3 headers handled via
    data_tests_functions.read_thermal_csv); all other files are read as text.
    """
    import mammos_entity as me

    from . import data_tests_functions as dtf

    def build(d):
        col = me.EntityCollection(description=f"raw data: {os.path.relpath(d, source_dir)}")
        for name in sorted(os.listdir(d)):
            p = os.path.join(d, name)
            if os.path.isdir(p):
                col[name] = build(p)
                continue
            if exclude_DOSCAR and name == "DOSCAR":
                continue
            if name.endswith(".csv"):
                try:
                    col[name] = np.loadtxt(p, delimiter=",")
                except ValueError:
                    df, _ = dtf.read_thermal_csv(p)
                    col[name] = df.to_numpy()
            else:
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    col[name] = fh.read()
        return col

    return build(source_dir)


def dataset_collection(source_dir, exclude_DOSCAR=False, include_raw=True):
    """Build one mammos-entity EntityCollection for an uppsala-data dataset:
    'intrinsic_properties' (from intrinsic_properties.yaml, via me.from_yaml),
    one 'thermal_MC_<n>' per UppASD/MC_*/thermal.csv (via me.from_csv), and
    'raw_data' (the folder tree). Written with collection.to_hdf5() it is
    readable back with me.from_hdf5() and linkable from NOMAD via
    m_annotations.hdf5.path (e.g. /intrinsic_properties/Js).
    """
    import mammos_entity as me

    intrinsic_file = os.path.join(source_dir, "intrinsic_properties.yaml")
    col = me.EntityCollection(
        description=f"UU DFT+MC dataset {os.path.basename(os.path.abspath(source_dir))}"
    )
    col["intrinsic_properties"] = me.from_yaml(intrinsic_file)
    uppasd = os.path.join(source_dir, "UppASD")
    if os.path.isdir(uppasd):
        for mc in sorted(d for d in os.listdir(uppasd) if d.startswith("MC_")):
            tcsv = os.path.join(uppasd, mc, "thermal.csv")
            if os.path.isfile(tcsv):
                col[f"thermal_{mc}"] = me.from_csv(tcsv)
    if include_raw:
        col["raw_data"] = raw_data_collection(source_dir, exclude_DOSCAR=exclude_DOSCAR)
    return col


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
