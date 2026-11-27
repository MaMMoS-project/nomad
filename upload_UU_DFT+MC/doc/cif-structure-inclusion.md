# CIF crystal-structure inclusion in NOMAD entries (research notes)

Tracked distillation of the local working notes `AI_cif_inclusion_intent.md` /
`AI_cif_inclusion_possible_solutions.md` (those stay gitignored, `AI_*.md`).
Related upstream: MaMMoS-project/to-nomad branch `feature/subcollection-subsections`;
reference demo upload `Nomad_cif demo_xVP6STRzQ-aYY1sSr-2cLg` (untracked copy in
this folder).

## Goal

Each UU dataset in `uppsala-data` ships a `structure.cif`. The NOMAD entry
overview should show the 3-D crystal structure, with the `.cif` also available
as raw data (provenance + download).

## What we do now (static workaround, implemented)

In `run_convert_UU_data_to_h5_and_creat_yaml.py` + `packages/nomad_archive_extras.py`:

1. `structure.cif` is copied to the **zip root** as `<dataset>_structure.cif`
   (raw data), in addition to staying under `raw_datasets/<dataset>/structure.cif`.
   This happens **independently of `--cif`** — the file is always shipped when present.
2. With `--cif auto/on` (default `off`), the cif is parsed **locally with ASE**
   (`structure_run_section()`) and embedded as an **inline top-level `run`
   section** in the `.archive.yaml` (`program: Structure File Importer`;
   `system[0].atoms` with lattice vectors/positions in meters, `labels`,
   `species` = atomic numbers, `periodic`), plus `metadata: {domain: dft}` so the
   server-side `OptimadeNormalizer` populates `metadata.optimade` — that is what
   the NOMAD GUI renders as the 3-D structure preview. With the run section
   embedded, `chemical_formula` is omitted (server-side ResultsNormalizer
   collision: "Could not populate compositional data as 'elements' is already
   defined", observed 2026-11-17 on a Fe12Co4N2 upload). With `--cif off`, the
   ELN `Chemical` base section's `chemical_formula` populates `results.material`.

## Limitation of the static approach

NOMAD **never reads the `.cif` at upload time**. The preview is the structure
embedded in the `.archive.yaml` at generation time: if the cif changes and only
the cif is re-uploaded, the preview does not change (the archive must be
regenerated), and the zip carries the structure twice (raw cif + parsed run
section), which can drift apart.

## Desired (non-static) behavior

NOMAD parses the uploaded `.cif` server-side (or a normalizer reads the raw
file from the upload) and builds the viewer from it — single source of truth.
Concretely either:

- fix/land `ElnWithStructureFile` (or an equivalent base section/normalizer) so an
  entry can reference `<dataset>_structure.cif` from the upload — the
  nomad-FAIR issue #2250 direction; or
- extend the `OptimadeNormalizer`/parser chain to parse raw `.cif` files present
  in the upload's raw files and populate `metadata.optimade` from them.

## Why `ElnWithStructureFile` is not used

On nomad-lab 2.0.0 that section fails deterministically during archive
persistence (`AttributeError: qualified_name`): the normalizer chain stores a
direct `atoms_ref` into a temporary run section that `ElnWithStructureFile`
then deletes, breaking serialization (unfixed on develop as of 2026-10-23;
upstream report drafted). The workaround matches the nomad-FAIR#2250 direction
(build the viewer from a System section embedded in the archive).

## Solution space (from the local notes)

1. **Upstream bug fix in FAIRmat-NFDI/nomad** — ideal long-term (fix the
   `ElnWithStructureFile` crash / add server-side cif parsing).
2. **Custom NOMAD Python schema plugin** — for Oasis/custom deployments; a
   plugin normalizer could parse raw cif files at processing time.
3. **Native multi-entry upload + reference linking** — entry references the cif
   entry instead of embedding.
4. **Current production workaround (implemented)** — client-side ASE parse into
   an inline `run` section at archive-generation time; the coupling is explicit
   (the embedded run section is generated from the same `structure.cif` that goes
   into the zip, in one step) and documented: changing the cif requires
   regenerating the archive.

## Checklist when testing an upload

- zip contains: `<dataset>.h5`, `<dataset>.archive.yaml`,
  `<dataset>_structure.cif` at the root, plus `raw_datasets/<dataset>/...` raw
  tree and `packages/*.py`;
- with `--cif on`: `.archive.yaml` has top-level `run` (atoms with lattice/
  positions in m) and `metadata: {domain: dft}`, and no `chemical_formula`;
- with `--cif off` (default): `data.chemical_formula` set (ELN Chemical base),
  no `run` section — Material section + element search are populated server-side;
- after upload on NOMAD, the entry overview shows the 3-D structure preview
  (`--cif on` only); compare with the reference demo upload
  `Nomad_cif demo_xVP6STRzQ-aYY1sSr-2cLg`.
