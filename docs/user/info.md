# Extract and Visualize Pangenome Information ℹ️

The info subcommand extracts summary information from PPanGGOLiN **.h5 pangenome files** and generates interactive HTML
reports. These reports support quick content comparison of each pangenome.

## Info command line usage ️

```shell
panorama info -i <pangenome_list.tsv> -o <output_directory> [--status] [--content]
```

## Output

- [**status_info.html**](#status-info) — Status of pangenome processing steps (annotation, clustering, etc.)

- [**content_info.html**](#content-info) — Numerical summary: genomes, genes, gene families, modules, etc.

## Key options

| Option       | Description                                                  |
|--------------|--------------------------------------------------------------|
| --status     | Extract and export the status (booleans) of each pangenome.  |
| --content    | Extract and export structural and numerical content metrics. |

default If no flags are provided, all (status, content, parameters, metadata) are extracted.

```{warning}
`--parameters` and `--metadata` outputs are not settled yet. We are currently working on a useful output.
Please add `--status` and/or `--content` to don't get an error.
```

## Exploring the Reports

### Status info

Shows whether each processing step was completed:

Example columns: Genomes_Annotated, Genes_Clustered, RGP_Predicted, etc.

Features:

- Radio button filters for boolean values.
- TSV download of filtered results.

### Content info

Displays statistics such as:

- Number of genes, genomes, gene families, modules.
- Frequencies and standard deviations.
- Partition composition (persistent, shell, cloud).

Features:

- Column visibility toggles.
- Range sliders for numeric filtering.
- TSV export of filtered view.
