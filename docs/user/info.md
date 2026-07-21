# Extract and Visualize Pangenome Information

The info subcommand extracts summary information from PPanGGOLiN **.h5 pangenome files** and generates interactive HTML
reports. These reports support quick content comparison of each pangenome.

## Info command line usage

```shell
panorama info -i <pangenome_list.tsv> -o <output_directory> [--status] [--content]
```

## Output

- [**status_info.html**](#status-info) — Status of pangenome processing steps (annotation, clustering, etc.)

- [**content_info.html**](#content-info) — Numerical summary: genomes, genes, gene families, modules, etc.

## Key options

| Option    | Description                                                  |
|-----------|--------------------------------------------------------------|
| --status  | Extract and export the status (booleans) of each pangenome.  |
| --content | Extract and export structural and numerical content metrics. |

Default: if no flags are provided, both `--status` and `--content` are extracted.

```{note}
`--parameters` and `--metadata` are not yet available. Use `--status` and/or `--content` for now.
```

## Exploring the Reports

(status-info)=

### Status info

Shows whether each processing step was completed:

Example columns: Genomes_Annotated, Genes_Clustered, RGP_Predicted, etc.

Features:

- Radio button filters for boolean values.
- TSV download of filtered results.

(content-info)=

### Content info

Displays statistics such as:

- Number of genes, genomes, gene families, modules.
- Frequencies and standard deviations.
- Partition composition (persistent, shell, cloud).

Features:

- Column visibility toggles.
- Range sliders for numeric filtering.
- TSV export of filtered view.
