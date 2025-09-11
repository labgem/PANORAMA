# 🧬 Conserved Spots Comparison Across Pangenomes

The compare_spots command identifies and analyzes conserved genomic spots across multiple pangenomes by comparing their
gene family composition and genomic organization patterns.
This analysis builds upon previously
[computed spots](https://ppanggolin.readthedocs.io/en/latest/user/RGP/rgpAnalyses.html#spot-prediction) from individual
pangenomes and uses **Gene Family Relatedness Relationship (GFRR) metrics** to identify regions that are conserved
across different bacterial populations. Optionally, it can integrate systems detection results to analyze biological
systems within conserved regions.

## ⚙️ Conserved Spots Detection Workflow

The conserved spots comparison process runs as follows:

1. 📂 Load and Validate Pangenomes
    - Multiple pangenomes are loaded from .h5 files based on a .tsv file.
    - Each pangenome is validated to ensure that spots and RGPs have been computed.
2. 📊 Create Spots Graph
    - All spots from all pangenomes are represented as nodes in a unified NetworkX graph.
    - Each spot is characterized by its bordering gene families.
3. 🧮 Compute GFRR-based Edges
    - For each pair of spots from different pangenomes:
        - Gene families at spot borders are extracted and compared.
        - GFRR metrics (min_gfrr and max_gfrr) are computed based on shared families.
    - Edges are added between spots that exceed both GFRR cutoff thresholds.
4. 🔗 Cluster Conserved Spots
    - Graph clustering algorithms identify groups of similar spots that represent conserved genomic regions across
      pangenomes based on the selected GFRR metric.
5. 🧬 Systems Integration (Optional)
    - When [enabled](../detection.md), systems analysis creates linkage graphs showing relationships between biological
      systems through their association with conserved spots.

6. 💾 Write Results to Files
    - Conserved spots are saved as detailed TSV files and optional graph formats (GEXF, GraphML) for visualization.

## 🚀 Command Line Usage*

Basic conserved spots comparison:

```shell
panorama compare_spots \
--pangenomes pangenomes.tsv \
--output conserved_spots_results \
--gfrr_metrics min_gfrr \
--gfrr_cutoff 0.8 0.8 \
--threads 8
```

With system analysis enabled:

```shell
panorama compare_spots \
--pangenomes pangenomes.tsv \
--output conserved_spots_results \
--systems \
--models defense_systems.tsv \
--sources defense_finder \
--gfrr_cutoff 0.8 0.8 \
--graph_formats gexf graphml \
--threads 8
```

## 📋 Key Options

| Shortcut | Argument        | Type           | Required/Optional       | Description                                                                   |
|----------|-----------------|----------------|-------------------------|-------------------------------------------------------------------------------|
| -p       | --pangenomes    | File path      | Required                | TSV file listing .h5 pangenomes with computed spots                           |
| -o       | --output        | Directory path | Required                | Output directory for conserved spots results                                  |
| —        | --gfrr_metrics  | String         | Optional                | GFRR metric for clustering: 'min_gfrr' (conservative) or 'max_gfrr' (liberal) |
| —        | --gfrr_cutoff   | Float Float    | Optional                | Two thresholds for min_gfrr and max_gfrr values (default: 0.8 0.8)            |
| —        | --dup_margin    | Float          | Optional                | Minimum ratio for multigenic family detection (default: 0.05)                 |
| —        | --systems       | Flag           | Optional                | Enable systems analysis within conserved spots                                |
| -m       | --models        | File path(s)   | Required with --systems | Path(s) to system model files (required with --systems)                       |
| -s       | --sources       | String(s)      | Required with --systems | System source names corresponding to models (required with --systems)         |
| —        | --canonical     | Flag           | Optional with --systems | Include canonical systems in analysis                                         |
| —        | --graph_formats | String(s)      | Optional                | Export graph formats: gexf, graphml                                           |

## Advanced Configuration Arguments

| Shortcut | Argument           | Type                 | Optional | Description                                                                              |
|----------|--------------------|----------------------|----------|------------------------------------------------------------------------------------------|
| —        | --cluster          | str (file path)      | True     | Tab-separated file with pre-computed clustering results (cluster_name\tfamily_id format) |
| —        | --tmpdir           | str (directory path) | True     | Directory for temporary files (default: /tmp)                                            |
| —        | --keep_tmp         | bool (flag)          | True     | Keep temporary files after completion                                                    |
| -c       | --cpus             | int                  | True     | Number of CPU threads for parallel processing (default: 1)                               |
| —        | --verbose          | int (choice)         | True     | Verbose level: 0 (warnings/errors), 1 (info), 2 (debug) (default: 1)                     |
| —        | --log              | str (file path)      | True     | Log output file (default: stdout)                                                        |
| -d       | --disable_prog_bar | bool (flag)          | True     | Disable the progress bars                                                                |
| —        | --force            | bool (flag)          | True     | Force writing in output directory and pangenome file                                     |

:::{note}
PANORAMA can perform the clustering step first thing, but it's also possible to use pre-computed clustering results with
the `--cluster` argument.
If you use let PANORAMA perform the clustering, you can look at the [Clustering](../clustering.md) section for more
details about options.
:::

## 📊 GFRR Metrics

| Metric   | Formula                                               | Description                                                        |
|----------|-------------------------------------------------------|--------------------------------------------------------------------|
| min_gfrr | shared_families / min(families_spot1, families_spot2) | Conservative metric requiring high overlap relative to smaller set |
| max_gfrr | shared_families / max(families_spot1, families_spot2) | Liberal metric allowing partial overlap relative to larger set     |

## 🎯 Sensitivity Control

The dual cutoff system provides fine-grained control over conservation stringency:

| Cutoff Level | min_gfrr | max_gfrr | Behavior                              |
|--------------|----------|----------|---------------------------------------|
| Strict       | 0.8      | 0.8      | High-confidence conserved spots only  |
| Moderate     | 0.6      | 0.7      | Balanced sensitivity and specificity  |
| Permissive   | 0.4      | 0.5      | Detects distant conservation patterns |

## 🗂 Output

PANORAMA generates multiple outputs: detailed spot information files, summary tables, and optional graph visualizations.

### File Organization

```
output_directory/
├── conserved_spots/
│ ├── conserved_spot_1.tsv
│ ├── conserved_spot_2.tsv
| |── ....................
│ └── conserved_spot_N.tsv
├── all_conserved_spots.tsv
├── conserved_spots.gexf (optional)
├── conserved_spots.graphml (optional)
├── systems_link_with_conserved_spots_louvain.gexf (optional)
└── systems_link_with_conserved_spots_mst.gexf (optional)
```

### Individual Conserved Spot Files

Each `conserved_spot_X.tsv` contains detailed RGP-level information:

| Column        | Description                                       |
|---------------|---------------------------------------------------|
| Spot_ID       | Original spot identifier from source pangenome    |
| Pangenome     | Source pangenome name                             |
| RGP_Name      | Region of Genomic Plasticity name within the spot |
| Gene_Families | Comma-separated list of gene families in the RGP  |

### Summary File

all_conserved_spots.tsv provides an overview of all conserved spots:

| Column            | Description                                      |
|-------------------|--------------------------------------------------|
| Conserved_Spot_ID | Unique identifier for the conserved spot group   |
| Spot_ID           | Individual spot identifier from source pangenome |
| Pangenome         | Source pangenome name                            |
| Num_RGPs          | Number of RGPs in this spot                      |
| Num_Gene_Families | Total number of gene families in this spot       |

### Conserved spots Graph Files (Optional)

When `--graph_formats` is enabled, additional graph files are generated:

- Conserved spots graph in GEXF format
- Conserved spots graph in GraphML format

Node attributes include conserved spot ID, pangenome name, spot ID, the number of gene families, and the number of
RGPs.
Edge attributes include GFRR metric and the number of shared gene families.

[PLACEHOLDER: Example conserved spots visualization across pangenomes]

### Systems Analysis Files (Optional)

When --systems is specified, genereate `systems_link_with_conserved_spots_louvain.gexf/graphml` Network graphs of
conserved system clusters. These graphs are generated using the Louvain algorithm.

[PLACEHOLDER: Systems linkage graph showing relationships through conserved spots]

