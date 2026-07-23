---
myst:
  html_meta:
    "description lang=en": |
      How PANORAMA's biological system detection algorithm works, step by step.
    "keywords": |
      Pangenomics, Python, Developer, Detection, Systems, Algorithm
---

# Detection Workflow

The `systems` subcommand detects biological systems in one or more pangenomes.
The entry point is {py:func}`panorama.systems.detection.launch`; the core logic lives in
{py:mod}`panorama.systems.detection`.

---

## Prerequisites

Detection requires that gene families have already been **annotated**: each
{py:class}`GeneFamily <panorama.geneFamily.GeneFamily>` must carry HMM-hit metadata
(written by `panorama annotation`) so the detector knows which families are candidates
for each functional unit.

The pre-flight check {py:func}`panorama.systems.detection.check_pangenome_detection`
verifies this and raises early if the annotation source is missing or systems have
already been computed (use `--force` to overwrite).

---

## Step-by-step pipeline

### 1. Load models

{py:class}`Models <panorama.systems.models.Models>` are parsed from the JSON model files
specified by `--models`. Each file yields one or more
{py:class}`Model <panorama.systems.models.Model>` objects, each containing a list of
{py:class}`FuncUnit <panorama.systems.models.FuncUnit>` objects and their associated
{py:class}`Family <panorama.systems.models.Family>` profiles.

### 2. Build the gene-context graph

For each {py:class}`Pangenome <panorama.pangenomes.Pangenome>`, PANORAMA calls PPanGGOLiN's
`ppanggolin.context.searchGeneContext.compute_gene_context_graph` to construct a
**gene-context graph** — a `networkx` graph in which:

- **nodes** are {py:class}`GeneFamily <panorama.geneFamily.GeneFamily>` objects
- **edges** connect families that co-occur within a genomic window (controlled by the
  `--window` and `--transitivity` parameters of the model)

This graph captures local genomic neighbourhood without committing to a fixed operon
boundary.

### 3. Find functional units in the graph

{py:func}`panorama.systems.detection.search_system_units` iterates over every
{py:class}`Model <panorama.systems.models.Model>` and calls
{py:func}`panorama.systems.detection.search_unit_in_context` for each
{py:class}`FuncUnit <panorama.systems.models.FuncUnit>`.

`search_unit_in_context` extracts **connected components** of the context graph that
contain at least one family matching the unit's profile, then calls:

- {py:func}`panorama.systems.detection.search_unit_in_cc` — tests each connected
  component for the presence of qualifying families
- {py:func}`panorama.systems.detection.search_unit_in_combination` — when a component
  contains multiple candidates, enumerates sub-combinations via
  {py:func}`panorama.systems.detection.get_subcombinations` to find the minimal
  satisfying set

The result is a dict `{FuncUnit → set[SystemUnit]}` of candidate detections.

### 4. Apply quorum rules

Before a full {py:class}`System <panorama.systems.system.System>` is assembled,
two filters are applied:

| Function | Rule checked |
|---|---|
| {py:func}`panorama.systems.detection.check_for_needed_units` | At least `min_mandatory` mandatory units and `min_total` total units must be present |
| {py:func}`panorama.systems.detection.check_for_forbidden_unit` | No `forbidden` units may be present |

If either check fails, the candidate is discarded.

### 5. Assemble System objects

{py:func}`panorama.systems.detection.search_for_system` takes the passing
`{FuncUnit → SystemUnit}` mapping and constructs a
{py:class}`System <panorama.systems.system.System>` object. Where multiple valid unit
combinations exist, {py:func}`panorama.systems.detection.get_system_unit_combinations`
enumerates them and each yields a separate `System`.

The full scan across all models and organisms is orchestrated by
{py:func}`panorama.systems.detection.search_systems` (single pangenome) and
{py:func}`panorama.systems.detection.search_systems_in_pangenomes` (parallel,
multi-pangenome via `ThreadPoolExecutor` + shared `Lock`).

### 6. Attach results to the pangenome

{py:func}`panorama.systems.detection.write_systems_to_pangenome` iterates the detected
{py:class}`System <panorama.systems.system.System>` objects and registers them on the
{py:class}`Pangenome <panorama.pangenomes.Pangenome>` via its `_system_getter` index.
{py:func}`panorama.systems.detection.write_systems_to_pangenomes` does the same across
all pangenomes in a {py:class}`Pangenomes <panorama.pangenomes.Pangenomes>` container.

### 7. Persist to HDF5

The `launch` function calls {py:func}`panorama.format.write_binaries.write_pangenome`,
which in turn calls {py:func}`panorama.format.write_binaries.write_systems_to_hdf5` to
serialise every `System` and `SystemUnit` into the pangenome's `.h5` file.
See [I/O Layer](io_layer.md) for the table schema.

---

## Data-flow summary

```
CLI args
  │
  ├─ check_detection_args()          validate & build need_info dict
  ├─ load_pangenomes()               read .h5 files (parallel)
  ├─ Models (from JSON files)        parse model definitions
  │
  └─ search_systems_in_pangenomes()  ── ThreadPoolExecutor ──►
        │  per pangenome:
        ├─ compute_gene_context_graph()   PPanGGOLiN context graph
        ├─ search_system_units()          find FuncUnit candidates
        ├─ check_for_needed_units()       quorum check
        ├─ check_for_forbidden_unit()     forbidden-unit check
        └─ search_for_system()            assemble System objects
              │
              └─ write_systems_to_pangenome()   attach to Pangenome
                    │
                    └─ write_pangenome()         persist to HDF5
```

---

## Key parameters and their effect

| Parameter | Model field | Effect |
|---|---|---|
| `--window` | `window` | Maximum genomic distance (in genes) between families in the same context |
| `--transitivity` | `transitivity` | Extends neighbourhood transitively — larger values widen the search |
| `--jaccard` | — | Minimum Jaccard similarity for context graph edges; lower values are more permissive |
| `--sensitivity` | — | Preset combining window + transitivity for convenience |

```{note}
Detection is run per-source: running `panorama systems` twice with different `--source`
values produces two independent sets of `System` objects stored under separate HDF5 groups.
Use `--force` to overwrite an existing source.
```
