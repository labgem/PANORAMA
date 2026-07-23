---
myst:
  html_meta:
    "description lang=en": |
      Description of PANORAMA's class architecture and its relationship with PPanGGOLiN
    "keywords": |
      Pangenomics, Python, Developer, Architecture, Classes, PPanGGOLiN
---

# Class Architecture

PANORAMA is built **on top of PPanGGOLiN**: every core domain object in PANORAMA extends a corresponding PPanGGOLiN
class rather than reimplementing it from scratch. This page maps the full class hierarchy, explains what each class
adds, and shows how the pieces fit together.

---

## Inheritance overview

```
PPanGGOLiN                          PANORAMA extension
──────────────────────────────────────────────────────
ppanggolin.pangenome.Pangenome  →   panorama.pangenomes.Pangenome
ppanggolin.geneFamily.GeneFamily →  panorama.geneFamily.GeneFamily
ppanggolin.region.Region        →   panorama.region.Region
ppanggolin.region.Spot          →   panorama.region.Spot
ppanggolin.region.Module        →   panorama.region.Module
ppanggolin.region.GeneContext   →   panorama.region.GeneContext
ppanggolin.metadata.MetaFeatures →  panorama.systems.system.SystemUnit
ppanggolin.metadata.MetaFeatures →  panorama.systems.system.System
panorama.pangenomes.Pangenome   →   panorama.pangenomes.Pangenomes  (container)
panorama.geneFamily.GeneFamily  →   panorama.geneFamily.Akin        (cross-pangenome link)
```

---

## Pangenome layer

### `Pangenome` ({py:mod}`panorama.pangenomes`)

Extends `ppanggolin.pangenome.Pangenome`.

PPanGGOLiN's `Pangenome` already manages organisms, contigs, genes, gene families, RGPs, spots, and modules. PANORAMA's
subclass {py:class}`panorama.pangenomes.Pangenome` adds the systems layer on top:

| Addition                          | Purpose                                                                                                                           |
|-----------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `_system_getter` / `_name2system` | Index detected {py:class}`System <panorama.systems.system.System>` objects by ID and name                                         |
| `systems_sources`                 | Track which annotation sources have been run                                                                                      |
| `add_file()` override             | Validates that the file is HDF5 and reads PANORAMA-specific status fields via {py:func}`panorama.format.read_binaries.get_status` |
| `status["systems"]`               | Extends PPanGGOLiN's status dict with a `"systems"` key                                                                           |

The constructor takes a `name` (used throughout outputs) and an optional `taxid`.

### `Pangenomes` ({py:mod}`panorama.pangenomes`)

A plain PANORAMA class (no PPanGGOLiN parent) that acts as the **multi-pangenome container**
used by every comparative workflow: {py:class}`panorama.pangenomes.Pangenomes`.

| Feature       | Details                                                                                                               |
|---------------|-----------------------------------------------------------------------------------------------------------------------|
| Storage       | Dict of `name →` {py:class}`Pangenome <panorama.pangenomes.Pangenome>`                                                |
| Thread safety | Accepts a shared `multiprocessing.Manager().Lock()` via {py:func}`panorama.utils.init_lock` for concurrent HDF5 reads |
| Loading       | Delegates per-pangenome loading to `ThreadPoolExecutor` (see {py:mod}`panorama.format.read_binaries`)                 |
| Input         | Populated from a TSV file via {py:func}`panorama.utils.check_tsv_sanity`                                              |

---

## Gene family layer

### `GeneFamily` ({py:mod}`panorama.geneFamily`)

Extends `ppanggolin.geneFamily.GeneFamily`.

Adds HMM-annotation state and cross-pangenome linkage on top of PPanGGOLiN's graph-node representation: {py:class}
`panorama.geneFamily.GeneFamily`.

| Addition                               | Purpose                                                                                                                                                        |
|----------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `_hmm`, `profile`, `optimized_profile` | Store the pyHMMER `HMM` object and search profiles used during annotation                                                                                      |
| `_units_getter` / `_systems_getter`    | Index the {py:class}`SystemUnit <panorama.systems.system.SystemUnit>` and {py:class}`System <panorama.systems.system.System>` objects that contain this family |
| `_akin`                                | Link to an {py:class}`Akin <panorama.geneFamily.Akin>` object representing equivalent families in other pangenomes                                             |

### `Akin` ({py:mod}`panorama.geneFamily`)

A pure PANORAMA class (no PPanGGOLiN parent) that groups {py:class}`GeneFamily <panorama.geneFamily.GeneFamily>` objects
from **different pangenomes** that have been determined to be functionally equivalent (by alignment or clustering). It
is the data structure behind cross-pangenome family comparison:
{py:class}`panorama.geneFamily.Akin`.

---

## Region layer

All four classes in {py:mod}`panorama.region` extend the matching PPanGGOLiN region class. At present they are thin
subclasses — their main purpose is to allow PANORAMA to inject additional attributes without patching PPanGGOLiN objects
directly.

| PANORAMA class                                        | Extends                          | Key additions                                                                 |
|-------------------------------------------------------|----------------------------------|-------------------------------------------------------------------------------|
| {py:class}`Region <panorama.region.Region>`           | `ppanggolin.region.Region` (RGP) | None — placeholder for future PANORAMA-specific state                         |
| {py:class}`Spot <panorama.region.Spot>`               | `ppanggolin.region.Spot`         | `pangenome` back-reference; `conserved_id` for cross-pangenome matching       |
| {py:class}`Module <panorama.region.Module>`           | `ppanggolin.region.Module`       | Tracks associated {py:class}`System <panorama.systems.system.System>` objects |
| {py:class}`GeneContext <panorama.region.GeneContext>` | `ppanggolin.region.GeneContext`  | Minimal extension used during system detection                                |

### `ConservedSpots` ({py:mod}`panorama.region`)

A pure PANORAMA class ({py:class}`panorama.region.ConservedSpots`) that groups {py:class}`Spot <panorama.region.Spot>`
objects from different pangenomes that occupy the same conserved genomic locus. It is the output of
`panorama compare_spots`.

---

## Systems layer

The systems layer is entirely PANORAMA-specific — PPanGGOLiN has no equivalent.

### Model definitions ({py:mod}`panorama.systems.models`)

These classes parse and validate model JSON files before detection runs.

| Class                                                   | Role                                                                                                                                            |
|---------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| {py:class}`Models <panorama.systems.models.Models>`     | Top-level container; maps `source → {name → Model}`                                                                                             |
| {py:class}`Model <panorama.systems.models.Model>`       | One system definition: name, quorum rules (`min_mandatory`, `min_total`), transitivity, window, and a list of `FuncUnit` objects                |
| {py:class}`FuncUnit <panorama.systems.models.FuncUnit>` | One functional unit within a model: presence type (`mandatory`, `accessory`, `forbidden`, `neutral`) and the `Family` profiles that can fill it |
| {py:class}`Family <panorama.systems.models.Family>`     | A named HMM/profile entry inside a `FuncUnit`; the leaf of the model tree                                                                       |

{py:class}`Model <panorama.systems.models.Model>` and {py:class}`FuncUnit <panorama.systems.models.FuncUnit>`
share behaviour via two private mixin bases (`_BasicFeatures`, `_ModFuFeatures`); {py:class}
`FuncUnit <panorama.systems.models.FuncUnit>` and {py:class}`Family <panorama.systems.models.Family>`
share a third (`_FuFamFeatures`).

### Detected objects ({py:mod}`panorama.systems.system`)

These classes hold the **results** of running detection against a pangenome.

| Class                                                               | Extends                            | Role                                                                                                                                                                                                                                                                      |
|---------------------------------------------------------------------|------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| {py:class}`SystemUnit <panorama.systems.system.SystemUnit>`         | `ppanggolin.metadata.MetaFeatures` | One detected functional unit: links a {py:class}`FuncUnit <panorama.systems.models.FuncUnit>` model to the {py:class}`GeneFamily <panorama.geneFamily.GeneFamily>` objects that satisfy it in a specific pangenome                                                        |
| {py:class}`System <panorama.systems.system.System>`                 | `ppanggolin.metadata.MetaFeatures` | One detected system: groups {py:class}`SystemUnit <panorama.systems.system.SystemUnit>` objects and tracks their associated {py:class}`Spot <panorama.region.Spot>`, {py:class}`Module <panorama.region.Module>`, and {py:class}`Region <panorama.region.Region>` objects |
| {py:class}`ClusterSystems <panorama.systems.system.ClusterSystems>` | —                                  | Groups {py:class}`System <panorama.systems.system.System>` objects from different pangenomes that represent the same biological system (used during comparison)                                                                                                           |

Both {py:class}`SystemUnit <panorama.systems.system.SystemUnit>` and {py:class}`System <panorama.systems.system.System>`
extend PPanGGOLiN's `MetaFeatures`, which provides the generic metadata key/value store used for TSV and HDF5
serialisation.

---

## How the layers connect

```
Pangenomes
 └─ Pangenome  (one per input .h5 file)
     ├─ GeneFamily  (extends ppanggolin GeneFamily)
     │    └─ Akin  (links equivalent families across pangenomes)
     ├─ Region / Spot / Module  (extend ppanggolin region classes)
     │    └─ ConservedSpots  (groups Spots across pangenomes)
     └─ System  (PANORAMA-only)
          └─ SystemUnit  (one per satisfied FuncUnit)
               └─ GeneFamily  (the families that fill the unit)

Models  (parsed before detection, not stored in Pangenome)
 └─ Model
      └─ FuncUnit
           └─ Family
```

A {py:class}`System <panorama.systems.system.System>` is created when the detection algorithm (in {py:mod}
`panorama.systems.detection`) finds that a set of gene families in a genomic window satisfies the quorum rules of a {py:
class}`Model <panorama.systems.models.Model>`. Each satisfied {py:class}`FuncUnit <panorama.systems.models.FuncUnit>`
becomes a {py:class}`SystemUnit <panorama.systems.system.SystemUnit>`; the {py:class}
`SystemUnit <panorama.systems.system.SystemUnit>` objects are collected into a {py:class}
`System <panorama.systems.system.System>` which is then attached to the {py:class}
`Pangenome <panorama.pangenomes.Pangenome>`.

---

## Key PPanGGOLiN touchpoints

| PPanGGOLiN symbol                                    | Where PANORAMA uses it                                                                                                      |
|------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `ppanggolin.pangenome.Pangenome`                     | Base of {py:class}`panorama.pangenomes.Pangenome`                                                                           |
| `ppanggolin.geneFamily.GeneFamily`                   | Base of {py:class}`panorama.geneFamily.GeneFamily`                                                                          |
| `ppanggolin.region.{Region,Spot,Module,GeneContext}` | Bases of the four classes in {py:mod}`panorama.region`                                                                      |
| `ppanggolin.metadata.MetaFeatures`                   | Base of {py:class}`SystemUnit <panorama.systems.system.SystemUnit>` and {py:class}`System <panorama.systems.system.System>` |
| `ppanggolin.context.searchGeneContext`               | Called by {py:mod}`panorama.systems.detection` to build the genomic-context graph used during detection                     |
| `ppanggolin.formats.*`                               | HDF5 read/write helpers extended by {py:mod}`panorama.format.read_binaries` and {py:mod}`panorama.format.write_binaries`    |

```{note}
There is a known circular-import risk between `pangenomes.py` and `systems/system.py`.
If you need to import one from the other, use a **function-local import** rather than a
top-level one. `tests/conftest.py` documents the current workaround.
```
