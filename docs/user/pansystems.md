# Pansystems - Complete Biological Systems Prediction Workflow

The `pansystems` command provides a comprehensive, all-in-one workflow to detect biological systems in multiple
pangenomes. This command combines [gene family annotation](annotation.md), [system detection](detection.md),
and [results visualization](write_systems.md) into a single streamlined process.

---

## Pansystems Workflow 🧪

The pansystems process integrates three main steps:

1. **Gene Family Annotation**: Uses either pre-computed metadata tables or HMM database searches.
   See [annotation documentation](annotation.md) for detailed information on annotation modes and options.

2. **System Detection**: Matches annotated families against functional models using genomic context analysis.
   See [detection documentation](detection.md) for details on the detection algorithm and sensitivity modes.

3. **Results Writing**: Generates comprehensive outputs including projection, association, and partition analysis.
   See [write_systems documentation](write_systems.md) for complete information on output options and file formats.

---

## Command Line Usage 🚀

### HMM-based Annotation and Detection

```bash
panorama pansystems \
    --pangenomes pangenomes.tsv \
    --source defense_finder \
    --hmm hmms.tsv \
    --models models.tsv \
    --mode fast \
    --k_best_hit 3 \
    --jaccard 0.8 \
    --projection \
    --association RGPs spots \
    --partition \
    --output results/ \
    --threads 16
```

### Table-based Annotation and Detection

```bash
panorama pansystems \
    --pangenomes pangenomes.tsv \
    --source KEGG \
    --table annotations.tsv \
    --models models.tsv \
    --jaccard 0.8 \
    --projection \
    --output results/ \
    --threads 8
```

---

## Key Arguments 🔑

### Required Arguments

| Shortcut | Argument       | Description                                            |
|----------|----------------|--------------------------------------------------------|
| `-p`     | `--pangenomes` | TSV file listing .h5 pangenome files                   |
| `-s`     | `--source`     | Name of the annotation source (e.g., `defense_finder`) |
| `-o`     | `--output`     | Output directory for results                           |

### Annotation Mode (Mutually Exclusive)

Must provide either `--table` or `--hmm`:

| Argument  | Description                                                  |
|-----------|--------------------------------------------------------------|
| `--table` | TSV linking pangenome names to pre-computed annotation files |
| `--hmm`   | HMM metadata TSV file (from `panorama utils --hmm`)          |

For detailed annotation options, see the [annotation command documentation](annotation.md#key-options).

### System Detection

| Argument        | Type  | Default | Description                                        |
|-----------------|-------|---------|----------------------------------------------------|
| `--models`      | Path  | -       | **Required.** Path to models list file             |
| `--jaccard`     | float | 0.8     | Minimum Jaccard similarity for context graph edges |
| `--sensitivity` | int   | 3       | Detection sensitivity (1-3, higher = more precise) |

For more detection parameters, see the [systems command documentation](detection.md#key-options).

### Output Options

| Argument        | Description                                                      |
|-----------------|------------------------------------------------------------------|
| `--projection`  | Project systems onto individual organisms                        |
| `--partition`   | Write partition heatmap showing system distribution              |
| `--association` | Associate systems with pangenome elements (RGPs, spots, modules) |
| `--proksee`     | Write proksee-compatible visualization files                     |

```{note}
For complete output options and file formats, see:

- [Projection documentation](projection.md)
- [Association documentation](association.md)
- [Partition documentation](partition.md)
```

### Performance Options

| Argument     | Type | Default | Description                |
|--------------|------|---------|----------------------------|
| `--threads`  | int  | 1       | Number of parallel threads |
| `--tmp`      | Path | auto    | Temporary directory        |
| `--keep_tmp` | flag | False   | Keep temporary files       |

---

## Output Structure

The pansystems command creates the same organized output structure as the individual commands. See
the [write_systems output documentation](write_systems.md#systems-analysis-output) for complete details on:

- System summary files
- Per-organism projections
- Association matrices and visualizations
- Partition heatmaps
- Interactive plots

---

## When to Use Pansystems vs Individual Commands

**Use `pansystems` when:**

- Running the complete workflow from annotation to visualization
- Processing multiple pangenomes with the same annotation source
- Need streamlined parameter validation and optimized data flow

**Use individual commands when:**

- Only need specific analysis steps
- Working with pre-annotated pangenomes
- Need fine-grained control over intermediate outputs
- Experimenting with different parameters for each step

---

## Integration Notes

The pansystems command ensures optimal integration between analysis steps:

- Validates parameter compatibility across all workflow components
- Eliminates needs for intermediate file management
- Provides unified progress reporting and error handling
- Maintains consistent source naming throughout the pipeline