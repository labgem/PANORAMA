# 🧬 Gene Family Alignment Across Pangenomes

The align command performs sequence alignment of gene families between multiple pangenomes using MMseqs2 to identify
homologous relationships and sequence similarities across different bacterial populations. This analysis supports both
targeted inter-pangenome comparisons (excluding intra-pangenome alignments) and comprehensive all-against-all alignments
that capture both inter- and intra-pangenome relationships.

## ⚙️ Alignment Workflow

The gene family alignment process runs as follows:

1. 📂 Load and Validate Pangenomes
    - Multiple pangenomes are loaded from .h5 files based on a .tsv file.
    - Each pangenome is validated to ensure gene families have been clustered and sequences are available.
2. 📝 Extract Gene Family Sequences
    - Gene family sequences are extracted from each pangenome and written to individual FASTA files.
    - Sequences are compressed and organized in temporary directories for processing.
3. 🗃️ Create MMseqs2 Databases
    - Individual sequence databases are created for each pangenome using MMseqs2.
    - For all-against-all mode, sequences are combined into a single unified database.
4. 🔍 Perform Sequence Alignments
    - **Inter-pangenome mode**: Pairwise alignments between all pangenome combinations, excluding self-alignments.
    - **All-against-all mode**: Comprehensive alignment including both inter- and intra-pangenome comparisons.
    - MMseqs2 search algorithms apply identity and coverage thresholds to identify significant matches.
5. 📊 Process and Merge Results
    - Alignment results are converted from binary format to human-readable TSV files.
    - Multiple alignment files are merged into consolidated output files.
6. 💾 Write Results to Files
    - Final alignment results are saved as detailed TSV files containing sequence similarity metrics.

## 🚀 Command Line Usage

Basic inter-pangenome alignment:

```shell
panorama align \
--pangenomes pangenomes.tsv \
--output alignment_results \
--inter_pangenomes \
--align_identity 0.8 \
--align_coverage 0.8 \
--threads 8
```

Comprehensive all-against-all alignment:

```shell
panorama align \
--pangenomes pangenomes.tsv \
--output alignment_results \
--all_against_all \
--align_identity 0.5 \
--align_coverage 0.8 \
--align_cov_mode 0 \
--threads 8 \
--keep_tmp
```

## 📋 Key Options

| Shortcut | Argument           | Type           | Required/Optional | Description                                                            |
|----------|--------------------|----------------|-------------------|------------------------------------------------------------------------|
| -p       | --pangenomes       | File path      | Required          | TSV file listing .h5 pangenomes with gene families and sequences       |
| -o       | --output           | Directory path | Required          | Output directory for alignment results                                 |
| —        | --inter_pangenomes | Flag           | Required (either) | Align gene families between pangenomes only (excludes intra-pangenome) |
| —        | --all_against_all  | Flag           | Required (either) | Align all gene families including intra-pangenome comparisons          |

## MMseqs2 Alignment Parameters

| Shortcut | Argument         | Type  | Optional | Description                                                                               |
|----------|------------------|-------|----------|-------------------------------------------------------------------------------------------|
| —        | --align_identity | Float | True     | Minimum identity percentage threshold (0.0-1.0, default: 0.5)                             |
| —        | --align_coverage | Float | True     | Minimum coverage percentage threshold (0.0-1.0, default: 0.8)                             |
| —        | --align_cov_mode | Int   | True     | Coverage mode: 0=query, 1=target, 2=shorter seq, 3=longer seq, 4=both, 5=all (default: 0) |

## Advanced Configuration Arguments

| Shortcut | Argument           | Type                 | Optional | Description                                                    |
|----------|--------------------|----------------------|----------|----------------------------------------------------------------|
| —        | --tmpdir           | str (directory path) | True     | Directory for temporary files (default: system temp directory) |
| —        | --keep_tmp         | bool (flag)          | True     | Keep temporary files after completion (useful for debugging)   |
| —        | --threads          | int                  | True     | Number of CPU threads for parallel processing (default: 1)     |

## 🎯 Alignment Modes

### Inter-Pangenome Alignment

This mode performs alignments **only between** different pangenomes, excluding intra-pangenome comparisons:

- **Use case**: Identifying shared gene families between populations
- **Results**: Focus on inter-population relationships

### All-Against-All Alignment

This mode performs comprehensive alignments including both inter- and intra-pangenome comparisons:

- **Use case**: Complete similarity analysis including within-population diversity
- **Results**: Complete gene family relationship matrix

## 📊 Parameter Guidelines

### Identity Thresholds

| Threshold | Use Case                          |
|-----------|-----------------------------------|
| 0.9-1.0   | Nearly identical sequences        |
| 0.7-0.9   | Highly similar homologs           |
| 0.5-0.7   | Moderate similarity               |
| 0.3-0.5   | Low similarity (use with caution) |

### Coverage Thresholds

| Threshold | Description               |
|-----------|---------------------------|
| 0.8-1.0   | High coverage requirement |
| 0.6-0.8   | Moderate coverage         |
| 0.4-0.6   | Permissive coverage       |

### Coverage Modes

| Mode | Target Coverage           |
|------|---------------------------|
| 0    | Query coverage            |
| 1    | Target coverage           |
| 2    | Shorter sequence coverage |

## 🗂 Output Files

PANORAMA generates alignment results in standardized TSV format with detailed similarity metrics.

### File Organization

```
output_directory/
├── inter_pangenomes.tsv        (inter-pangenome mode)
└── all_against_all.tsv         (all-against-all mode)
```

### Alignment Results Format

Each alignment file contains the following columns:

| Column    | Description                            | Example     |
|-----------|----------------------------------------|-------------|
| query     | Query gene family identifier           | PG1_FAM_001 |
| target    | Target gene family identifier          | PG2_FAM_045 |
| identity  | Percentage sequence identity (0.0-1.0) | 0.85        |
| qlength   | Query sequence length in amino acids   | 150         |
| tlength   | Target sequence length in amino acids  | 145         |
| alnlength | Alignment length in amino acids        | 142         |
| e_value   | E-value of the alignment               | 1.2e-45     |
| bits      | Bit score of the alignment             | 185.2       |

[//]: # (### Results Interpretation)

[//]: # (## 🔬 Biological Applications)
