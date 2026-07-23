---
myst:
  html_meta:
    "description lang=en": |
      How PANORAMA reads and writes HDF5 pangenome files and produces flat text outputs.
    "keywords": |
      Pangenomics, Python, Developer, HDF5, I/O, Format, Read, Write
---

# I/O Layer

PANORAMA stores all computed data in the same **HDF5 `.h5` file** that PPanGGOLiN uses,
extending it with PANORAMA-specific groups and tables. Flat text outputs are written
separately by {py:mod}`panorama.format.write_flat` and {py:mod}`panorama.systems.write_systems`.

---

## Module overview

| Module | Role |
|---|---|
| {py:mod}`panorama.format.read_binaries` | Read HDF5 files; parallel multi-pangenome loading |
| {py:mod}`panorama.format.write_binaries` | Write / update / erase system data in HDF5 files |
| {py:mod}`panorama.format.write_flat` | Write flat TSV/CSV outputs for gene families and spots |
| {py:mod}`panorama.systems.write_systems` | Write flat TSV outputs for detected systems |

---

## HDF5 structure

PANORAMA reuses the PPanGGOLiN HDF5 layout and adds one new top-level group:

```
/                             (PPanGGOLiN root)
├── status/                   extended with PANORAMA attributes
│    ├── systems              bool — True if any systems have been written
│    └── systems_sources      set of source name strings
│
├── annotations/              PPanGGOLiN (unchanged)
├── gene_families/            PPanGGOLiN (unchanged)
├── RGP/                      PPanGGOLiN (unchanged)
├── spots/                    PPanGGOLiN (unchanged)
├── modules/                  PPanGGOLiN (unchanged)
│
└── systems/                  PANORAMA addition
     └── <source_name>/       one group per annotation source
          ├── systems         table — one row per System
          ├── units           table — one row per (System, SystemUnit, GeneFamily) triple
          ├── canonical       table — one row per canonical System
          └── canonical_units table — one row per canonical (System, SystemUnit, GeneFamily) triple
```

### Table schemas

#### `systems` table

| Column | Type | Description |
|---|---|---|
| `ID` | `Int64` | System identifier |
| `name` | `String` | System name (from the model) |
| `canonical_count` | `Int64` | Number of canonical representations |

#### `units` table

| Column | Type | Description |
|---|---|---|
| `ID` | `Int64` | SystemUnit identifier |
| `name` | `String` | FuncUnit name |
| `geneFam` | `String` | GeneFamily name |
| `metadata_id` | `Int64` | Metadata record identifier |
| `metadata_source` | `String` | Annotation source name |

The `canonical` and `canonical_units` tables follow the same schemas as `systems` and
`units` respectively.

String column widths are computed dynamically before the table is created by
{py:func}`panorama.format.write_binaries.calculate_system_table_sizes` to avoid wasting
space while still fitting all values.

---

## Reading HDF5 files

### Status check

{py:func}`panorama.format.read_binaries.get_status` is called by
{py:meth}`Pangenome.add_file() <panorama.pangenomes.Pangenome.add_file>` for every
pangenome. It calls PPanGGOLiN's own `get_status` first, then reads the PANORAMA-specific
`status.systems` and `status.systems_sources` attributes to populate
`pangenome.status["systems"]` and `pangenome.status["systems_sources"]`.

### Loading a single pangenome

{py:func}`panorama.format.read_binaries.read_pangenome` is the all-in-one loader for a
single `.h5` file. It accepts a `need_info` dict that controls which data groups are
read:

| `need_info` key | What is loaded |
|---|---|
| `need_annotations` | Organism / gene annotations |
| `need_families` | Gene family graph |
| `need_families_info` | Per-family metadata (HMM hits) |
| `need_rgp` | Regions of genomic plasticity |
| `need_spots` | Spot structures |
| `need_modules` | Module structures |
| `need_systems` | PANORAMA system objects |

Internally it delegates to the matching PPanGGOLiN reader functions for all standard
groups and to {py:func}`panorama.format.read_binaries.read_systems` for the PANORAMA
systems group.

### Parallel loading of multiple pangenomes

{py:func}`panorama.format.read_binaries.load_pangenomes` loads a
{py:class}`Pangenomes <panorama.pangenomes.Pangenomes>` container using a
`ThreadPoolExecutor`. A shared `multiprocessing.Manager().Lock()` (initialised via
{py:func}`panorama.utils.init_lock`) serialises all concurrent HDF5 table reads because
PyTables is not thread-safe for simultaneous reads.

```python
# Typical call pattern in a launch() function
need_info = check_my_command_args(args)
pangenomes = load_pangenomes(
    pangenomes_file=args.pangenomes,
    need_info=need_info,
    disable_bar=args.disable_prog_bar,
    cpu=args.cpus,
)
```

---

## Writing HDF5 files

### Writing systems

{py:func}`panorama.format.write_binaries.write_systems_to_hdf5` is the low-level writer.
It is called by {py:func}`panorama.format.write_binaries.write_pangenome`, which also
handles the PPanGGOLiN base tables and updates the status group.

The write sequence is:

1. {py:func}`panorama.format.write_binaries.calculate_system_table_sizes` — scan all
   systems to determine maximum string lengths and row counts
2. {py:func}`panorama.format.write_binaries.create_system_tables` — create the four HDF5
   tables with the computed column sizes
3. {py:func}`panorama.format.write_binaries.write_system_data_to_tables` — iterate
   systems and write one row per `(System, SystemUnit, GeneFamily)` triple
4. {py:func}`panorama.format.write_binaries.write_pangenome_status` — update
   `status.systems` and `status.systems_sources`

### Erasing systems

{py:func}`panorama.format.write_binaries.erase_pangenome` removes the HDF5 group for a
given source (called when `--force` is passed to re-run detection).

---

## Flat file outputs

### Gene families and spots ({py:mod}`panorama.format.write_flat`)

Invoked by the `write` subcommand. Writes TSV files describing:
- HMM annotation results per gene family
- Spot border compositions
- Conserved spot memberships

### Systems ({py:mod}`panorama.systems.write_systems`)

Invoked by the `write_systems` subcommand (and by `pansystems`). Writes TSV files
describing:
- Per-pangenome system summaries
- System-to-organism projections
- System-to-spot / system-to-module associations
- Partition heatmap data

```{note}
All flat files are regenerated from in-memory objects — they do not read back from HDF5.
Re-running `write_systems` with different filters produces different outputs without
modifying the `.h5` file.
```
