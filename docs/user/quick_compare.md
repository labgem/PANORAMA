# Pangenome Comparison Analysis

PANORAMA provides two specialized commands for comparing biological features across multiple pangenomes using [Gene
Family Relatedness Relationship (GFRR) metrics](#gfrr_overview).

---

## Compare spots across pangenomes

The `compare_spots` command identifies conserved genomic spots across multiple pangenomes by comparing gene family
composition at spot borders.

### Compare spots command-line usage 🚀

```bash
# Basic conserved spots comparison
panorama compare_spots \
    --pangenomes pangenomes.tsv \
    --output conserved_spots_results \
    --gfrr_cutoff 0.8 0.8 \
    --threads 8

# With systems analysis
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

### Compare spots key Arguments

| Argument         | Description                                                 |
|------------------|-------------------------------------------------------------|
| `--pangenomes`   | TSV file listing .h5 pangenomes with computed spots         |
| `--output`       | Output directory for conserved spots results                |
| `--gfrr_cutoff`  | Two thresholds for min_gfrr and max_gfrr (default: 0.8 0.8) |
| `--gfrr_metrics` | GFRR metric for clustering: `min_gfrr` or `max_gfrr`        |
| `--systems`      | Enable systems analysis within conserved spots              |
| `--models`       | Path(s) to system model files (required with --systems)     |
| `--sources`      | System source names (required with --systems)               |

For complete parameter details and GFRR metrics explanation, see the [compare_spots documentation](compare_spots.md).

---

## Compare systems across pangenomes

The `compare_systems` command identifies conserved biological systems across multiple pangenomes and generates
distribution heatmaps.

### Quick Usage

```bash
# Basic systems comparison with heatmaps
panorama compare_systems \
    --pangenomes pangenomes.tsv \
    --models defense_systems.tsv \
    --sources defense_finder \
    --output systems_comparison_results \
    --heatmap \
    --threads 8

# Full analysis with conserved systems clustering
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

### Compare systems key Arguments

| Argument               | Description                                                 |
|------------------------|-------------------------------------------------------------|
| `--pangenomes`         | TSV file listing .h5 pangenomes with detected systems       |
| `--models`             | Path(s) to system model files                               |
| `--sources`            | Name(s) of systems sources                                  |
| `--output`             | Output directory for comparison results                     |
| `--heatmap`            | Generate distribution heatmaps                              |
| `--gfrr_cutoff`        | Two thresholds for min_gfrr and max_gfrr (default: 0.5 0.8) |
| `--gfrr_models_cutoff` | GFRR thresholds for model gene families (default: 0.4 0.6)  |
| `--gfrr_metrics`       | GFRR metric for clustering conserved systems                |

For complete parameter details and GFRR metrics explanation, see
the [compare_systems documentation](compare_systems.md).

---
(gfrr_overview)=
## GFRR Metrics Overview

Both commands use Gene Family Relatedness Relationship metrics to assess conservation:

| Metric   | Formula                                       | Usage                               |
|----------|-----------------------------------------------|-------------------------------------|
| min_gfrr | shared_families / min(families_A, families_B) | Conservative: high overlap required |
| max_gfrr | shared_families / max(families_A, families_B) | Liberal: partial overlap accepted   |

### Sensitivity Control

| Cutoff Level | min_gfrr | max_gfrr | Behavior                              |
|--------------|----------|----------|---------------------------------------|
| Strict       | 0.8      | 0.8      | High-confidence conservation only     |
| Moderate     | 0.6      | 0.7      | Balanced sensitivity and specificity  |
| Permissive   | 0.4      | 0.5      | Detects distant conservation patterns |

---

## Output Summary

### Compare Spots Outputs

- Individual conserved spot files (`conserved_spot_X.tsv`)
- Summary file (`all_conserved_spots.tsv`)
- Optional graph files (GEXF, GraphML)
- Optional systems linkage graphs

### Compare Systems Outputs

- Interactive heatmaps showing system distribution
- Optional network graphs of conserved system clusters
- Optional summary tables of conserved systems

For detailed output descriptions, see the respective documentation:

- [Compare Spots Output](compare_spots.md#output)
- [Compare Systems Output](compare_systems.md#output)

---

## Prerequisites

### For Compare Spots

- Pangenomes must
  have [computed spots and RGPs](https://ppanggolin.readthedocs.io/en/latest/user/RGP/rgpAnalyses.html#spot-prediction)
- Optional: [detected systems](detection.md) if using `--systems` flag

### For Compare Systems

- Pangenomes must have [detected systems](detection.md) for the specified sources
- Models files corresponding to the systems sources
