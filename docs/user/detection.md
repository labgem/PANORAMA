# System Detection Based on Models
The `systems` command enables the detection of biological systems in pangenomes using predefined functional models.

This detection relies on gene family [annotations](annotation.md#gene-family-annotation), 
and a [model](../modeler/modeling.md#models) file defining the presence/absence of a specific function and the genomic organization.

## Model Detection Workflow

The detection process runs as follows:

1. Load Pangenomes
    
   Based on a `.tsv` file, `.h5` pangenomes are loaded with the required annotation and metadata sources.

2. Load System Models

   The models are parsed from a list provided by `--models`.

3. Search for System Units

   For each functional unit of each model:
   1. Gene families are matched based on annotation metadata. 
   2. A context graph is built based on the gene neighborhood (window, transitivity). 
   3. Jaccard similarity filters edges in the graph.
   4. Connected components are checked for necessary and forbidden families.

4. Assemble Systems

    Functional units are grouped into systems if they satisfy:
    - presence/abscence rules
    - restrain distance

5. Write Systems to File

   Detected systems are saved back into the pangenome `.h5` file, under the given source name.

## Command Line Usage

System detection command is used as such:
```shell
panorama systems \
--pangenomes pangenomes.tsv \
--models models.tsv \
--source defense_finder \
--annotation_sources  defense_finder CasFinder\
--jaccard 0.8 \
--sensitivity 3 \
--threads 8
```
## Key Options

| Shortcut | Argument             | Description                                                                |
|----------|----------------------|----------------------------------------------------------------------------|
| -p       | --pangenomes         | TSV file listing .h5 pangenomes                                            |
| -m       | --models             | Path to the models list file (see panorama utils --models)                 |
| -s       | --source             | Name of the source used to annotate and detect systems                     |
| —        | --annotation_sources | List of annotation sources to use (defaults to the one from --source)      |
| —        | --jaccard            | Minimum Jaccard similarity to keep an edge in context graph (default: 0.8) |
| —        | --sensitivity        | Sensitivity mode: 1, 2, or 3. Higher = more precise, slower (default: 3)   |
| —        | --threads            | Number of threads to use for parallel model evaluation                     |

<!--
## 🔍 Sensitivity Modes


| Level | Description                                                                            |

|-------|----------------------------------------------------------------------------------------|

| 1     | Global filtering of genomic context, faster, less sensitive                            |

| 2     | Global filtering context within each functional unit combination, moderate sensitivity |

| 3     | Local filtering for each combination (highest sensitivity, slowest)                    |
-->

## Output

PANORAMA integrates multiple outputs: textual, graph-based representations and figures to summarize results.

See the [write_systems](write_systems.md#systems-analysis-output)
