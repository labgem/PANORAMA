# System Projection on Genomes

The `write_systems` command enables the projection of systems, previously detected at the pangenome level
(see [`systems` command](detection.md)), onto individual genomes. Projection relies on system detection results and
the genomic context of gene families within organisms.

## Projection Workflow

The projection process has been optimized and proceeds as follows:

### 1. Load Detected Systems and Metadata

- Detected systems from the .h5 pangenome file are loaded
- Required metadata and gene families are retrieved
- System-to-family mappings are established for efficient processing

### 2. Build Gene Context Components

For each organism and functional unit, the workflow uses a **component-based approach** instead of graph
construction:

1. **Identify Model Genes**: Extract genes belonging to system families in each organism
2. **Group by Contig**: Organize genes by their chromosomal/plasmid location
3. **Extract Windows**: Use `extract_contig_window()` to identify genomic regions containing system genes within the
   specified window size
4. **Create Components**: Each window becomes a component containing all genes (model + context) within that region

This approach directly identifies co-localized gene clusters.

### 3. Project System Units

Each system unit is evaluated in organisms through the following steps:

#### Unit Requirements Validation

- **Family Requirements**: Check if required families from the model are present
- **Completeness Calculation**: Determine what fractions of model families are found
- **Context Analysis**: Identify additional families within the same genomic context

#### System State Classification

Components are classified into three genomic organization states:

- **strict**: All model families are found within the same connected component/window
- **split**: Model families are present but spread across multiple disconnected components
- **extended**: All model families are in the same context with additional intervening families

#### Gene Categorization

Each projected gene is categorized as:

- **model**: Gene belongs to a family defined in the system model
- **context**: Gene is co-localized with model genes but not part of the system definition
- **filtered**: Gene was excluded during filtering steps

### 4. Aggregate and Filter Projections

The projection includes advanced filtering options:

#### Standard Projection

- Collects all valid projections for each organism
- Calculates completeness metrics
- Maintains full system context information

#### One-Unit-Per-Family Filtering

New optimization that handles overlapping system units:

- **Overlap Resolution**: When multiple units contain the same gene family, keeps only the unit with the highest
  completeness
- **Overlapping Units Tracking**: Records information about filtered units in `overlapping_units` column
- **System Elimination Options**:
    - `eliminate_filtered_systems`: Remove entire systems if any model families were filtered
    - `eliminate_empty_systems`: Remove systems with no remaining model families

### 5. Write Output

Projection results are written as TSV files with improved organization and metadata.
See [Output Files](#output-files) for details on the organization and contents.

## Projection command Line Usage

### Basic Projection

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder \
    --projection \
    --threads 8 \
    --output results/
```

## Advanced Options

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder immune_system \
    --projection \
    --association RGPs spots \                # Associate systems with RGPs and hotspots
    --partition \                             # Write partition heatmap files
    --canonical \                             # Project canonical versions of systems
    --organisms organism_A organism_B \       # Project only these organisms
    --threads 16 \
    --force \
    --output results/
```

## Projection command Line Arguments

### Projection-specific keys

| Argument        | Type | Default | Description                                            |
|-----------------|------|---------|--------------------------------------------------------|
| `--projection`  | flag | False   | Enable the projection of systems onto genomes          |
| `--organisms`   | list | None    | List of organisms to project (defaults to all)         |
| `--canonical`   | flag | False   | Also project canonical versions of systems             |

### Required Arguments

| Argument       | Type | Description                                     |
|----------------|------|-------------------------------------------------|
| `--pangenomes` | Path | TSV file listing pangenome .h5 files to process |
| `--output`     | Path | Output directory for projection results         |
| `--models`     | Path | Path(s) to model list files                     |
| `--sources`    | str  | Name(s) of the systems sources                  |

### Optional Arguments

| Argument        | Type | Default | Description                                            |
|-----------------|------|---------|--------------------------------------------------------|
| `--projection`  | flag | False   | Enable the projection of systems onto genomes          |
| `--organisms`   | list | None    | List of organisms to project (defaults to all)         |
| `--canonical`   | flag | False   | Also project canonical versions of systems             |
| `--threads`     | int  | 1       | Number of parallel threads to use                      |
| `--force`       | flag | False   | Overwrite existing projection files                    |

## Projection Output Files

Output is organized in the specified `--output` directory with subdirectories for each pangenome and source combination:

```
output/
├── pangenome_1/
│   └── source_1/
│       ├── systems.tsv                    # Pangenome summary
│       └── projection/
│           ├── organism_A.tsv            # Per-organism detailed results
│           ├── organism_B.tsv
│           └── ...
└── pangenome_2/
    └── source_1/
        ├── systems.tsv
        └── projection/
            └── ...
```

## 1. Pangenome Systems Summary (`systems.tsv`)

This file provides a high-level summary of all detected systems across the pangenome:

| Column               | Description                                                                              |
|----------------------|------------------------------------------------------------------------------------------|
| system number        | Unique numeric ID for the system                                                         |
| system name          | Name of the system (corresponds to model name)                                           |
| functional unit name | Name of the functional unit within the system                                            |
| organism             | Organism name where the system is detected                                               |
| model_GF             | Comma-separated list of gene families encoding system functions                          |
| context_GF           | Comma-separated list of gene families found in genomic context but not part of the model |
| partition            | Pangenome partition of the system (persistent, shell, cloud, or combinations)            |
| completeness         | Average proportion of model families found across organisms (0.0-1.0)                    |
| strict               | Number of organisms with strict genomic organization                                     |
| split                | Number of organisms with split genomic organization                                      |
| extended             | Number of organisms with extended genomic organization                                   |

**Additional columns** (when using `--association`):

- **RGPs**: Associated Regions of Genomic Plasticity
- **spots**: Associated hotspots of genome evolution
- **modules**: Associated functional modules

## 2. Organism Projection Files (`projection/<organism>.tsv`)

Each organism gets a detailed file with gene-level projections:

| Column               | Description                                                       |
|----------------------|-------------------------------------------------------------------|
| system number        | Unique system ID                                                  |
| system name          | System name from the model                                        |
| functional unit name | Functional unit name                                              |
| subsystem number     | ID for the genomic component/subgraph                             |
| organism             | Organism name                                                     |
| gene family          | Gene family identifier                                            |
| partition            | Pangenome partition (persistent/shell/cloud)                      |
| annotation           | Functional annotation from metadata                               |
| secondary_names      | Alternative names for the gene family                             |
| gene.ID              | Unique gene identifier                                            |
| gene.name            | Gene name/locus tag                                               |
| contig               | Contig/chromosome name                                            |
| start                | Gene start position                                               |
| stop                 | Gene stop position                                                |
| strand               | Gene orientation (+/-)                                            |
| is_fragment          | Whether gene is fragmented                                        |
| category             | Gene category: `model`, `context`, or `filtered`                  |
| genomic organization | System organization: `strict`, `split`, or `extended`             |
| completeness         | Proportion of model families present in this organism             |
| product              | Gene product description                                          |
| overlapping_units    | Information about overlapping units (format: `unit:completeness`) |

**Additional columns** (when using `--association`):

| Column     | Description                |
|------------|----------------------------|
| RGPs       | Associated RGP identifier  |
| spots name | Associated spot identifier |
