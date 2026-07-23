---
myst:
  html_meta:
    "description lang=en": |
      How PANORAMA's pangenome comparison workflows work for spots and systems.
    "keywords": |
      Pangenomics, Python, Developer, Comparison, Spots, Systems, GFRR
---

# Compare Workflow

PANORAMA provides two comparison subcommands that share a common graph-based strategy:

| Subcommand | Entry point | What it compares |
|---|---|---|
| `compare_spots` | {py:func}`panorama.compare.spots.launch` | Genomic islands (RGP spots) across pangenomes |
| `compare_systems` | {py:func}`panorama.compare.systems.launch` | Detected biological systems across pangenomes |

Both rely on a **Gene Family Relatedness Ratio (GFRR)** metric computed in
{py:mod}`panorama.compare.utils` and follow the same three-phase structure:
build a graph → score edges → cluster into conserved groups.

---

## Shared infrastructure ({py:mod}`panorama.compare.utils`)

### GFRR — the similarity metric

The GFRR between two objects (spots or systems) from different pangenomes is the fraction
of their gene families that are considered equivalent (via {py:class}`Akin <panorama.geneFamily.Akin>`
links established by `panorama align` or `panorama cluster`):

```
GFRR(A, B) = |shared Akin families| / aggregation(|families in A|, |families in B|)
```

The aggregation function is controlled by `--gfrr_metrics` (`min_gfrr`, `max_gfrr`,
`mean_gfrr`). `--gfrr_cutoff` sets the threshold below which an edge is dropped.

{py:func}`panorama.compare.utils.compute_gfrr` computes the score for a single pair;
it is called from both comparison modules.

### Clustering

{py:func}`panorama.compare.utils.cluster_on_gfrr` applies connected-component clustering
on the filtered graph. Each connected component becomes one
{py:class}`ConservedSpots <panorama.region.ConservedSpots>` or
{py:class}`ClusterSystems <panorama.systems.system.ClusterSystems>` group.

---

## Spots comparison ({py:mod}`panorama.compare.spots`)

(validate_load_spot)=
### 1. Validate and load

{py:func}`panorama.compare.spots.check_compare_spots_args` validates arguments and
returns the `need_info` dict. Pangenomes are loaded in parallel by
{py:func}`panorama.format.read_binaries.load_pangenomes` with spots and RGPs enabled.

### 2. Build per-pangenome spot graphs

{py:func}`panorama.compare.spots.create_pangenome_spots_graph` builds a `networkx`
graph for each pangenome where nodes are {py:class}`Spot <panorama.region.Spot>` objects
and edges connect spots that share gene families within that pangenome.

(build_cross_spot)=
### 3. Build the cross-pangenome graph

{py:func}`panorama.compare.spots.create_spots_graph` merges the per-pangenome graphs
and adds cross-pangenome edges by calling
{py:func}`panorama.compare.spots.compute_gfrr_edges`, which iterates over every pair of
pangenomes and scores each spot–spot combination using GFRR.

#### Optional: overlay systems

If `--systems` is passed, {py:func}`panorama.compare.spots.add_systems_info` annotates
each spot node with the systems detected within it, and
{py:func}`panorama.compare.spots.create_pangenome_system_graph` /
{py:func}`panorama.compare.spots.create_systems_graph` build a parallel systems layer
linked to the spot graph via
{py:func}`panorama.compare.spots.graph_systems_link_with_conserved_spots`.

(write_outputs_spots)=
### 4. Cluster and write outputs

{py:func}`panorama.compare.utils.cluster_on_gfrr` groups connected spots into
{py:class}`ConservedSpots <panorama.region.ConservedSpots>` objects.
{py:func}`panorama.compare.spots.write_conserved_spots` writes the TSV and optional
graph outputs.

(data_flow_spots)=
### Data-flow summary

```
load_pangenomes()
  │
  ├─ create_pangenome_spots_graph()   per-pangenome graph
  ├─ create_spots_graph()             merge + GFRR cross-pangenome edges
  │    └─ compute_gfrr_edges()        score every spot pair
  │
  ├─ [optional] add_systems_info()    overlay systems on spots
  │
  ├─ cluster_on_gfrr()                connected components → ConservedSpots
  └─ write_conserved_spots()          TSV + graph outputs
```

---

## Systems comparison ({py:mod}`panorama.compare.systems`)

(validate_load_systems)=
### 1. Validate and load

{py:func}`panorama.compare.systems.check_compare_systems_args` validates that the number
of `--sources` matches `--models` and returns the `need_info` dict. Pangenomes are loaded
with systems and family metadata enabled.

### 2. Build per-pangenome system graphs

{py:func}`panorama.compare.systems.create_pangenome_system_graph` builds a graph for
each pangenome where nodes are {py:class}`System <panorama.systems.system.System>` objects
and edges connect systems that share gene families within that pangenome.
{py:func}`panorama.compare.systems.add_system_metadata_to_graph` annotates nodes with
system metadata for downstream output.

(build_cross_systems)=
### 3. Build the cross-pangenome graph

{py:func}`panorama.compare.systems.create_systems_graph` merges per-pangenome graphs
and {py:func}`panorama.compare.systems.compute_gfrr_edges` scores every system–system
pair across pangenomes using GFRR.

(write_outputs_systems)=
### 4. Cluster and write outputs

{py:func}`panorama.compare.utils.cluster_on_gfrr` groups connected systems into
{py:class}`ClusterSystems <panorama.systems.system.ClusterSystems>` objects.
{py:func}`panorama.compare.systems.write_conserved_systems` writes TSV outputs.

Optionally, {py:func}`panorama.compare.systems.generate_heatmap` and
{py:func}`panorama.compare.systems.create_pangenome_systems_heatmaps` produce Bokeh
HTML heatmaps of system distributions across pangenomes.

(data_flow_systems)=
### Data-flow summary

```
load_pangenomes()
  │
  ├─ create_pangenome_system_graph()  per-pangenome graph
  ├─ create_systems_graph()           merge + GFRR cross-pangenome edges
  │    └─ compute_gfrr_edges()        score every system pair
  │
  ├─ cluster_on_gfrr()                connected components → ClusterSystems
  ├─ write_conserved_systems()        TSV outputs
  └─ [optional] generate_heatmap()    Bokeh HTML visualisation
```

---

## Parallel execution

Both workflows wrap their per-pangenome graph construction in a `ThreadPoolExecutor`
and use a shared `multiprocessing.Manager().Lock()` (via {py:func}`panorama.utils.init_lock`)
to serialise concurrent HDF5 reads. The number of workers is controlled by `--cpus`.

```{note}
Alignment or clustering (`panorama align` / `panorama cluster`) must be run before
comparison so that {py:class}`Akin <panorama.geneFamily.Akin>` links exist between
gene families from different pangenomes. Without them every GFRR score is 0.
```
