# System Partition Analysis ➗

The `write_systems` command enables the creation of partition visualizations for pangenome systems, showing how genetic
systems are distributed across different partition categories (persistent, shell, cloud, etc.) in pangenomes. Partition
analysis provides interactive heatmap visualizations to understand the conservation patterns and evolutionary dynamics
of genetic systems across organisms.

## Partition Workflow 🧪

The partition analysis creates comprehensive heatmap visualizations showing system distribution patterns:

### 1. Load Systems and Partition Data

- Detected systems from the .h5 pangenome file are loaded
- System projection data with partition information is retrieved
- Partition categories are reconciled for systems spanning multiple partitions

### 2. Build Partition Matrix

For partition analysis, the workflow:

1. **Preprocess Partition Data**: Create pivot tables mapping organisms × systems to partition categories
2. **Reconcile Partitions**: Handle systems that may span multiple partition types using the `conciliate_partition`
   function
3. **Create Matrix**: Generate organism × system matrix with partition assignments
4. **Handle Missing Data**: Fill gaps with "Not_found" for systems absent in specific organisms

### 3. Generate Interactive Visualizations

The partition analysis produces comprehensive interactive heatmaps with:

- **Main Heatmap**: Color-coded matrix showing system partition assignments across organisms
- **Left Bar Chart**: Stacked bars showing partition distribution by system
- **Top Bar Chart**: Bar chart showing system counts per organism
- **Color Bar**: Legend showing partition type color coding
- **Hover Tooltips**: Interactive information display

## Partition Command Line Usage 🚀

### Basic Partition Analysis

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder \
    --partition \
    --threads 8 \
    --output results/
```

### Combined with Other Analyses

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models models.tsv \
    --sources defense_finder \
    --projection \               # Also project systems onto genomes
    --partition \               # Write partition heatmaps
    --association RGPs \        # Associate with RGPs
    --canonical \               # Include canonical versions
    --threads 16 \
    --output results/
```

### Multiple Sources Analysis

```bash
panorama write_systems \
    --pangenomes pangenomes.tsv \
    --models defense_finder.tsv padloc.tsv \
    --sources defense_finder padloc \
    --partition \
    --output_formats html png \ # Generate both HTML and PNG outputs
    --threads 16 \
    --output results/
```

For more information on other analysis options, see the documentation about [projection](projection.md)
and [association](association.md).

## Partition command Line Arguments ⚙️

### Partition-Specific Key 🔑

| Argument      | Type | Description                                                      |
|---------------|------|------------------------------------------------------------------|
| `--partition` | flag | Enable partition heatmap generation for systems across organisms |

### Required Arguments

| Argument       | Type | Description                                     |
|----------------|------|-------------------------------------------------|
| `--pangenomes` | Path | TSV file listing pangenome .h5 files to process |
| `--output`     | Path | Output directory for partition results          |
| `--models`     | Path | Path(s) to model list files                     |
| `--sources`    | str  | Name(s) of the systems sources                  |

### Optional Arguments

| Argument           | Type | Default  | Description                               |
|--------------------|------|----------|-------------------------------------------|
| `--output_formats` | list | ["html"] | Visualization output format customization |
| `--threads`        | int  | 1        | Number of parallel threads to use         |
| `--force`          | flag | False    | Overwrite existing partition files        |
| `--canonical`      | flag | False    | Include canonical versions of systems     |

**Output Format Options**

The visualization outputs can be customized:

- **HTML**: Interactive Bokeh plots with hover tooltips and zoom capabilities (default)
- **PNG**: Static high-resolution images for publications and reports

## Partition Categories 📊

The partition analysis classifies systems into five main categories based on the partition of gene families coding
a function in the system:

| Partition Category        | Color Code | Description                                               |
|---------------------------|------------|-----------------------------------------------------------|
| **Persistent**            | 🟠 Orange  | All systems gene families are persistent                  |
| **Persistent\|Accessory** | 🔴 Red     | Systems gene families spanning both all paritions         |
| **Accessory**             | 🟣 Purple  | Systems gene families spanning cloud and shell partitions |
| **Shell**                 | 🟢 Green   | All systems gene families are shell                       |
| **Cloud**                 | 🔵 Blue    | All systems gene families are cloud                       |

## Partition Output Files 📁

Partition analysis creates organized output in the specified `--output` directory:

```
output/
├── pangenome_1/
│   └── source_1/
│       ├── partition.html                     # Interactive partition heatmap
│       ├── partition.png                      # Static partition heatmap (if requested)
└── pangenome_2/
    └── source_1/
        └── ...
```

### Interactive Visualization Files

#### Partition Heatmap (`partition.html`)

Interactive Bokeh visualization containing:

**Main Components:**

- **Central Heatmap**: Organism × System matrix with partition-colored cells
- **Left Stacked Bar Chart**: Partition distribution by system (shows composition of systems across organisms)
- **Top Bar Chart**: System count per organism (green bars showing total systems per organism)
- **Right Color Bar**: Legend showing partition category colors and labels

**Interactive Features:**

- **Hover Tooltips**: Display system name, organism, and partition category
- **Zoom/Pan Tools**: Navigate large matrices with many organisms/systems
- **Selection Tools**: Highlight specific regions of interest
- **Export Tools**: Save plots or extract data
- **Responsive Layout**: Automatically adjusts to matrix dimensions

**Visual Design:**

- **Color Scheme**: Carefully chosen colors for accessibility and biological interpretation
- **Cell Borders**: White borders for clear cell separation
- **Axis Labels**: Rotated organism names for readability
- **Legend**: Comprehensive partition category explanations

### Stacked Bar Interpretation

**Left Stacked Bars (Systems):**

- **Composition Analysis**: Shows how each system is distributed across partition categories
- **System Flexibility**: Systems with multiple partition types are evolutionarily plastic
- **Conservation Assessment**: Predominantly red systems are highly conserved, blue systems are rare

**Top Bars (Organisms):**

- **System Richness**: Total height indicates the diversity of systems in each organism
- **Comparative Analysis**: Organisms can be compared for their system repertoire size
- **Ecological Adaptation**: System-rich organisms may inhabit complex or challenging environments

## Technical Details 🔧

### Visualization Technology

- **Bokeh Framework**: Web-based interactive plotting library
- **Responsive Design**: Automatically scales with data matrix size
- **Performance Optimization**: Efficient rendering of large matrices
- **Cross-platform**: Works in any modern web browser

### Color Mapping

- **Categorical Mapper**: Uses `CategoricalColorMapper` for discrete partition categories
- **Accessibility**: Colors chosen for colorblind-friendly visualization
- **Consistent Encoding**: Same colors used across all components for coherent interpretation
