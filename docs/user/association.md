# System Association with Pangenome Elements

The `write_systems` command enables the creation of associations between systems and other pangenome elements such as
RGPs (Regions of Genomic Plasticity), spots, and modules. Association analysis provides correlation matrices and
visualizations to analyze relationships between systems and various pangenome components.

## Association Workflow

The association process creates comprehensive correlation matrices and coverage analyses:

### 1. Load Systems and Elements

- Detected systems from the .h5 pangenome file are loaded
- Target pangenome elements (RGPs, spots, modules) are retrieved
- System-to-element mappings are established for correlation analysis

### 2. Build Association Matrices

For each association type, the workflow:

1. **Extract Element Relationships**: Map systems to their associated RGPs, spots, or modules
2. **Create Correlation Matrix**: Generate system × element matrices showing co-occurrence patterns
3. **Calculate Coverage**: Compute how well each element is covered by associated systems
4. **Compute Frequencies**: Determine element occurrence frequencies in organisms in which systems also exist

### 3. Generate Visualizations

The association analysis produces interactive correlation matrix heatmaps with:

- **Main Heatmap**: Color-coded correlation matrix showing system-element associations
- **Coverage Plot**: Visual representation of how systems cover pangenome elements
- **Frequency Plot**: Display of element frequencies across genomes
- **Bar Charts**: Summary statistics for systems and elements
- **Color Bars**: Legend and scaling information

## Association Command Line Usage

### Basic Association Analysis

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder \
    --association RGPs \
    --threads 8 \
    --output results/
```

### Multiple Association Types

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder \
    --association all \          # Analyze RGPs, spots, and modules
    --threads 16 \
    --output results/
```

### Combined with Other Analyses

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder immune_system \
    --projection \               # Also project systems onto genomes
    --association spots modules \ # Associate with spots and modules
    --partition \               # Write partition heatmaps
    --canonical \               # Include canonical versions
    --threads 16 \
    --output results/
```

For more information on other analysis options, here is the documentation about [projection](projection.md).

## Association command line arguments

### Association-Specific key

| Argument        | Type | Choices                           | Description                                  |
|-----------------|------|-----------------------------------|----------------------------------------------|
| `--association` | list | `RGPs`, `spots`, `modules`, `all` | Pangenome elements to associate with systems |

**Association Types:**

- **RGPs**: Regions of Genomic Plasticity - variable genomic regions
- **spots**: Hotspots of genome evolution - frequently variable loci
- **modules**: Functional modules - co-evolving gene clusters
- **all**: Analyze all three association types

### Required Arguments

| Argument       | Type | Description                                     |
|----------------|------|-------------------------------------------------|
| `--pangenomes` | Path | TSV file listing pangenome .h5 files to process |
| `--output`     | Path | Output directory for projection results         |
| `--models`     | Path | Path(s) to model list files                     |
| `--sources`    | str  | Name(s) of the systems sources                  |

### Optional Arguments

| Argument           | Type | Default  | Description                               |
|--------------------|------|----------|-------------------------------------------|
| `--output_formats` | list | ["html"] | Visualization output format customization |
| `--threads`        | int  | 1        | Number of parallel threads to use         |
| `--force`          | flag | False    | Overwrite existing projection files       |

**Output Format Options**

The visualization outputs can be customized (currently the HTML format is default):

- **HTML**: Interactive Bokeh plots (default)
- **PNG**: Static high-resolution images (if requested)

## Association Output Files

Association analysis creates organized output in the specified `--output` directory:

```
output/
├── pangenome_1/
│   └── source_1/
│       ├── association.tsv                    # Main association matrix
│       ├── correlation_RGPs.html              # RGP correlation visualization
│       ├── correlation_RGPs.png               # RGP correlation visualization (PNG)
│       ├── correlation_spots.html             # Spot correlation visualization  
│       ├── correlation_spots.png              # Spot correlation visualization  (PNG)
│       ├── correlation_modules.html           # Module correlation visualization
│       ├── correlation_modules.png            # Module correlation visualization (PNG)
│       ├── rgp_to_systems.tsv                 # RGP-system mappings
│       ├── spot_to_systems.tsv                # Spot-system mappings
│       └── module_to_systems.tsv              # Module-system mappings
└── pangenome_2/
    └── source_1/
        └── ...
```

### 1. Main Association File (`association.tsv`)

Primary file containing system-element associations:

| Column        | Description                                                                |
|---------------|----------------------------------------------------------------------------|
| system number | Unique numeric ID for the system (index)                                   |
| system name   | Name of the system from the model                                          |
| families      | Comma-separated list of gene families in the system                        |
| RGPs          | Comma-separated list of associated RGP names (if `--association RGPs`)     |
| spots         | Comma-separated list of associated spot IDs (if `--association spots`)     |
| modules       | Comma-separated list of associated module IDs (if `--association modules`) |

### 2. Element-to-Systems Mapping Files

#### RGP Mapping (`rgp_to_systems.tsv`)

| Column       | Description                                     |
|--------------|-------------------------------------------------|
| name         | RGP identifier (index)                          |
| systems_ID   | Comma-separated list of associated system IDs   |
| systems_name | Comma-separated list of associated system names |
| coverage     | Proportion of RGP families covered by systems   |
| frequency    | Frequency of RGP occurrence across organisms    |

#### Spot Mapping (`spot_to_systems.tsv`)

| Column       | Description                                     |
|--------------|-------------------------------------------------|
| name         | Spot identifier (index, format: `spot_<ID>`)    |
| systems_ID   | Comma-separated list of associated system IDs   |
| systems_name | Comma-separated list of associated system names |
| coverage     | Proportion of spot families covered by systems  |
| frequency    | Frequency of spot occurrence across organisms   |

#### Module Mapping (`module_to_systems.tsv`)

| Column       | Description                                      |
|--------------|--------------------------------------------------|
| name         | Module identifier (index, format: `module_<ID>`) |
| systems_ID   | Comma-separated list of associated system IDs    |
| systems_name | Comma-separated list of associated system names  |
| coverage     | Proportion of module families covered by systems |
| frequency    | Frequency of module occurrence across organisms  |

### 3. Interactive Visualization Files

#### Correlation Matrix Plots (`correlation_<type>.html`)

Interactive Bokeh visualizations containing:

**Main Components:**

- **Central Heatmap**: System × Element correlation matrix with hover tooltips
- **Left Bar Chart**: System occurrence counts across elements
- **Top Bar Chart**: Element occurrence counts across systems
- **Right Color Bar**: Legend showing correlation intensity scale

**Lower Panels:**

- **Frequency Plot**: Element frequency across genomes (blue color scale)
- **Coverage Plot**: System coverage of elements (red color scale)
- **Color Bar Legends**: Scaling information for frequency and coverage plots

**Interactive Features:**

- **Hover Tooltips**: Show detailed information on mouseover
- **Zoom/Pan**: Navigate large correlation matrices
- **Save Tools**: Export plots or data
- **Responsive Layout**: Adapts to different screen sizes

## Data Interpretation

### Coverage Metrics

**Coverage** represents the number of element gene families covered by the associated systems:

- **1.0 (dark red)**: Associated systems cover all families in the element
- **0.5 (light red)**: Systems explain half of the element's families
- **0.0 (white)**: No overlap between element and system families

Depending on the element type, coverage is computed differently:

- **RGP**: use the intersection between the gene of the RGP and gene of an associated system present in the genome of the RGP
- **Spot**: use gene families between borders of the spot
- **Module**: use the intersection of the module's gene families and the gene families of an associated system

### Frequency Metrics

**Frequency** indicates how common elements and associated systems appear across genomes:

- **High Frequency (dark blue)**: Element found in many genomes
- **Low Frequency (light blue)**: Element found in few genomes
- **Interpretation**: Helps identify core vs. accessory genomic features

### Correlation Intensity

**Correlation values** in the heatmap show system-element co-occurrence:

- **High Values (dark colors)**: high number of system types associated with the element
- **Low Values (light colors)**: Weak or no association
- **Pattern Analysis**: Reveals functional relationships and genomic organization

## Technical Details

### Visualization Technology

- **Bokeh Framework**: Interactive web-based visualizations
- **Responsive Design**: Adapts to correlation matrix size
- **Color Palettes**: Carefully chosen for accessibility and clarity
- **Export Options**: Support for HTML and PNG formats

### Performance Optimization

The association analysis includes several optimizations:

- **Parallel Processing**: Multi-threaded computation of associations
- **Memory Efficiency**: Streaming processing of large pangenomes
- **Vectorized Operations**: Pandas-based matrix operations
- **Component-based Analysis**: Direct mapping without graph construction