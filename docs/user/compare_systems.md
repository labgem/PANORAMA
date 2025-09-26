# Systems Comparison Across Pangenomes 🧬

The `compare_systems` command identifies and analyzes conserved biological systems across multiple pangenomes by
comparing their gene family composition and computing similarity metrics.
This analysis builds upon previously [detected systems from individual pangenomes](detection.md) and uses **Gene
Family Relatedness Relationship (GFRR) metrics** to identify systems that are conserved across different bacterial
populations. The analysis generates visualizations showing system distribution patterns and creates graphs of conserved
system clusters.

## Systems Comparison Workflow ⚙️

The systems comparison process runs as follows:

1. **Load and Validate Pangenomes**
    - Multiple pangenomes are loaded from .h5 files based on a .tsv file.
    - Each pangenome is validated to ensure that systems have been detected for the specified sources.

2. **Create Systems**
    - All systems from all pangenomes are represented as nodes in a
      unified [NetworkX](https://networkx.org/documentation/stable/) graph.
    - Each system is characterized by its gene families and model families for similarity assessment.

3. Compute GFRR-based Edges
    - For each pair of systems from different pangenomes:
        - Model gene families are compared using GFRR metrics.
        - If model families exceed thresholds, all gene families are compared.
    - Edges are added between systems that exceed both GFRR cutoff thresholds.

4. Cluster Conserved Systems

   Graph clustering algorithms
   ([Louvain](https://networkx.org/documentation/stable/reference/algorithms/community.html#module-networkx.algorithms.community.louvain))
   identify groups of similar systems that represent conserved biological systems across pangenomes based on the
   selected GFRR metric.

5. Generate Visualizations

   Heatmaps showing system distribution patterns across pangenomes are generated in HTML format for interactive
   exploration.

6. Write Results to Files

   Conserved systems are saved as graph files (GEXF, GraphML) and summary tables for further analysis and visualization.

## System comparison command Line Usage 🚀

Basic systems comparison with heatmap generation:

```shell
panorama compare_systems \
--pangenomes pangenomes.tsv \
--models defense_systems.tsv \
--sources defense_finder \
--output systems_comparison_results \
--heatmap \
--threads 8
```

Full analysis with conserved systems clustering:

```shell
panorama compare_systems \
--pangenomes pangenomes.tsv \
--models defense_systems.tsv cas_systems.tsv \
--sources defense_finder CasFinder \
--output systems_comparison_results \
--heatmap \
--gfrr_metrics min_gfrr_models \
--gfrr_cutoff 0.8 0.8 \
--gfrr_models_cutoff 0.2 0.2 \
--graph_formats gexf graphml \
--threads 8
```

### Key Options 📋

| Shortcut | Argument             | Type                   | Optional | Description                                                                                         |
|----------|----------------------|------------------------|----------|-----------------------------------------------------------------------------------------------------|
| -p       | --pangenomes         | str (file path)        | False    | TSV file listing .h5 pangenomes with detected systems                                               | 
| -m       | --models             | List[str] (file paths) | False    | Path(s) to system model files (must match --sources order)                                          |
| -s       | --sources            | List[str]              | False    | Name(s) of systems sources (must match --models order)                                              |
| -o       | --output             | str (directory path)   | False    | Output directory for comparison results                                                             |
| —        | --gfrr_cutoff        | List[float] (2 values) | True     | Two thresholds for min_gfrr and max_gfrr values (default: 0.5 0.8)                                  |
| —        | --heatmap            | bool (flag)            | True     | Generate heatmaps showing system distribution across pangenomes                                     |
| —        | --gfrr_metrics       | str (choice)           | True     | GFRR metric for clustering conserved systems (min_gfrr_models, max_gfrr_models, min_gfrr, max_gfrr) |
| —        | --gfrr_models_cutoff | List[float] (2 values) | True     | GFRR thresholds for model gene families (default: 0.4 0.6)                                          |
| —        | --graph_formats      | List[str]              | True     | Export graph formats: gexf, graphml                                                                 |
| —        | --canonical          | bool (flag)            | True     | Include canonical system versions in analysis                                                       |

### Advanced Configuration Arguments

| Shortcut | Argument           | Type                 | Optional | Description                                                                              |
|----------|--------------------|----------------------|----------|------------------------------------------------------------------------------------------|
| —        | --cluster          | str (file path)      | True     | Tab-separated file with pre-computed clustering results (cluster_name\tfamily_id format) |
| —        | --method           | str (choice)         | True     | MMSeqs2 clustering method: linclust or cluster (default: linclust)                       |
| —        | --tmpdir           | str (directory path) | True     | Directory for temporary files (default: /tmp)                                            |
| —        | --keep_tmp         | bool (flag)          | True     | Keep temporary files after completion                                                    |
| -c       | --cpus             | int                  | True     | Number of CPU threads for parallel processing (default: 1)                               |
| —        | --verbose          | int (choice)         | True     | Verbose level: 0 (warnings/errors), 1 (info), 2 (debug) (default: 1)                     |
| —        | --log              | str (file path)      | True     | Log output file (default: stdout)                                                        |
| -d       | --disable_prog_bar | bool (flag)          | True     | Disable the progress bars                                                                |
| —        | --force            | bool (flag)          | True     | Force writing in output directory and pangenome file                                     |

```{note}
PANORAMA can perform the clustering step first thing, but it's also possible to use pre-computed clustering results with
the `--cluster` argument.
If you use let PANORAMA perform the clustering, you can look at the [Clustering](../clustering.md) section for more
details about options.
```

### GFRR Metrics for Systems

| Metric          | Target Families     | Description                                        |
|-----------------|---------------------|----------------------------------------------------|
| min_gfrr_models | Model families only | Conservative metric using core functional families |
| max_gfrr_models | Model families only | Liberal metric using core functional families      |
| min_gfrr        | All families        | Conservative metric using complete gene repertoire |
| max_gfrr        | All families        | Liberal metric using complete gene repertoire      |

### Cutoff Configuration

The dual-cutoff system provides hierarchical filtering:

| Filtering Stage | Cutoffs            | Purpose                                         |
|-----------------|--------------------|-------------------------------------------------|
| Model families  | gfrr_models_cutoff | Primary filter using core functional genes      |
| All families    | gfrr_cutoff        | Secondary filter using complete gene repertoire |

### Recommended settings

- Strict: gfrr_models_cutoff=[0.5, 0.5], gfrr_cutoff=[0.8, 0.8]
- Moderate: gfrr_models_cutoff=[0.3, 0.3], gfrr_cutoff=[0.6, 0.7]
- Permissive: gfrr_models_cutoff=[0.2, 0.2], gfrr_cutoff=[0.4, 0.5]

## Output 📂

PANORAMA generates multiple outputs: interactive heatmaps, network graphs, and summary tables for comprehensive systems
analysis.

### File Organization

```
output_directory/
├── heatmap_number_systems.html
├── heatmap_normalized_systems.html
├── conserved_systems.gexf (optional)
├── conserved_systems.graphml (optional)
└── conserved_systems.tsv (optional)
```

### Files description

#### Heatmap Visualizations

Interactive HTML heatmaps showing system distribution patterns:

| File                            | Description                                       |
|---------------------------------|---------------------------------------------------|
| heatmap_number_systems.html     | Raw counts of each system type per pangenome      |
| heatmap_normalized_systems.html | Normalized percentages showing relative abundance |

[//]: # (Test)
<iframe src="../pictures/heatmap_number_sys.html"></iframe>

[PLACEHOLDER: Heatmap showing system distribution across multiple pangenomes]

[PLACEHOLDER: Normalized heatmap showing relative system abundance patterns]

#### Conserved System Clustering

##### Network Graphs
When `--gfrr_metrics` and `--graph_formats` are specified, genereate `conserved_systems.gexf/graphml` Network graphs of
conserved system clusters.
Node attributes include system metadata, pangenome information, and cluster assignments
Edge attributes contain GFRR similarity scores and the number of shared gene families.

[PLACEHOLDER: Network graph of conserved systems clusters with different colors]

##### Summary Tables
When conserved systems clustering is performed:
conserved_systems.tsv: Tabular summary of identified conserved system clusters
