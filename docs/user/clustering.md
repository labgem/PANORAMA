(clustering)=

# Gene Family Clustering Across Pangenomes

The cluster command groups related gene families from multiple pangenomes into homologous clusters based on sequence
similarity using MMseqs2. This analysis identifies gene families that share common evolutionary origins across different
bacterial populations, enabling comparative genomics studies and the construction of meta-pangenomes. The clustering
supports both fast (linclust) and sensitive (cluster) methods to accommodate different accuracy and performance
requirements.

## Clustering Workflow

The gene family clustering process runs as follows:

1. Load and Validate Pangenomes
    - Multiple pangenomes are loaded from .h5 files based on a .tsv file.
    - Each pangenome is validated to ensure gene families have been clustered and sequences are available.
2. Extract Gene Family Sequences
    - Gene family sequences are extracted from each pangenome and written to compressed FASTA files.
    - All sequences are combined while maintaining pangenome-specific identifiers.
3. Create Unified MMseqs2 Database
    - All gene family sequences are combined into a single MMseqs2 database.
    - Database indexing optimizes subsequent clustering operations.
4. Perform Sequence Clustering
    - **Linclust method**: Fast linear-time clustering suitable for large datasets with moderate sensitivity
      requirements.
    - **Cluster method**: Sensitive clustering with comprehensive sequence comparisons for high-accuracy results.
    - MMseqs2 applies identity and coverage thresholds to group similar sequences.
5. Process Clustering Results
    - Binary clustering results are converted to human-readable TSV format.
    - Cluster IDs are assigned to group-related gene families.
6. Write Results to Files
    - Final clustering results are saved with cluster assignments and membership details.

## Clustering command Line Usage

- Fast clustering with linclust:

   ```shell
   panorama cluster \
   --pangenomes pangenomes.tsv \
   --output clustering_results \
   --method linclust \
   --cluster_identity 0.8 \
   --cluster_coverage 0.8 \
   --threads 8
   ```

- Sensitive clustering with comprehensive parameters:

   ```shell
   panorama cluster \
   --pangenomes pangenomes.tsv \
   --output clustering_results \
   --method cluster \
   --cluster_identity 0.5 \
   --cluster_coverage 0.8 \
   --cluster_sensitivity 7.5 \
   --cluster_max_seqs 500 \
   --threads 8 \
   --keep_tmp
   ```

### Key Options

| Shortcut | Argument     | Type           | Required/Optional | Description                                                      |
|----------|--------------|----------------|-------------------|------------------------------------------------------------------|
| -p       | --pangenomes | File path      | Required          | TSV file listing .h5 pangenomes with gene families and sequences |
| -o       | --output     | Directory path | Required          | Output directory for clustering results                          |
| -m       | --method     | String         | Required          | Clustering method: 'linclust' (fast) or 'cluster' (sensitive)    |

### MMseqs2 Core Clustering Parameters

| Shortcut | Argument           | Type  | Optional | Description                                                                       |
|----------|--------------------|-------|----------|-----------------------------------------------------------------------------------|
| —        | --cluster_identity | Float | True     | Minimum sequence identity threshold (0.0-1.0, default: 0.5)                       |
| —        | --cluster_coverage | Float | True     | Minimum coverage threshold (0.0-1.0, default: 0.8)                                |
| —        | --cluster_cov_mode | Int   | True     | Coverage mode: 0=query, 1=target, 2=shorter, 3=longer, 4=both, 5=all (default: 0) |
| —        | --cluster_eval     | Float | True     | E-value threshold for sequence similarity (default: 1e-3)                         |

### Cluster method-specific parameters

| Argument               | Type  | Description                                             |
|------------------------|-------|---------------------------------------------------------|
| --cluster_sensitivity  | Float | Search sensitivity (higher = more sensitive but slower) |
| --cluster_max_seqs     | Int   | Maximum sequences per cluster representative            |
| --cluster_min_ungapped | Int   | Minimum ungapped alignment score                        |

### Advanced Parameters (Both Methods)

| Argument                 | Type | Description                                                                             |
|--------------------------|------|-----------------------------------------------------------------------------------------|
| --cluster_comp_bias_corr | Int  | Compositional bias correction: 0=disabled, 1=enabled                                    |
| --cluster_kmer_per_seq   | Int  | Number of k-mers per sequence                                                           |
| --cluster_align_mode     | Int  | Alignment mode: 0=auto, 1=score only, 2=extended, 3=both, 4=fast                        |
| --cluster_max_seq_len    | Int  | Maximum sequence length                                                                 |
| --cluster_max_reject     | Int  | Maximum number of rejected sequences                                                    |
| --cluster_mode           | Int  | Clustering algorithm: 0=Set Cover, 1=Connected Component, 2=Greedy, 3=Greedy Low Memory |

### Advanced Configuration Arguments

| Shortcut | Argument   | Type                 | Optional | Description                                                    |
|----------|------------|----------------------|----------|----------------------------------------------------------------|
| —        | --tmpdir   | str (directory path) | True     | Directory for temporary files (default: system temp directory) |
| —        | --keep_tmp | bool (flag)          | True     | Keep temporary files after completion (useful for debugging)   |
| —        | --threads  | int                  | True     | Number of CPU threads for parallel processing (default: 1)     |

## Clustering Methods Comparison

### Linclust (Fast Method)

**Algorithm**: Linear time complexity clustering
**Performance**: Fast execution suitable for large datasets
**Sensitivity**: Moderate - may miss distant relationships

| Parameter Range   | Recommendation           |
|-------------------|--------------------------|
| Identity: 0.7-0.9 | High confidence clusters |
| Identity: 0.5-0.7 | Balanced approach        |
| Identity: 0.3-0.5 | Permissive clustering    |

### Cluster (Sensitive Method)

**Algorithm**: Comprehensive pairwise comparisons
**Performance**: Slower but more thorough
**Sensitivity**: High - detects distant homologs

| Parameter   | Recommended Range | Impact                    |
|-------------|-------------------|---------------------------|
| Sensitivity | 4.0-6.0           | Standard sensitivity      |
| Sensitivity | 6.0-8.0           | High sensitivity (slower) |
| Max Seqs    | 100-300           | Balanced performance      |
| Max Seqs    | 500-1000          | Comprehensive but slow    |

## Parameter Guidelines

### Identity and Coverage Combinations

| Identity | Coverage | Clustering Behavior | Biological Interpretation  |
|----------|----------|---------------------|----------------------------|
| 0.9      | 0.9      | Very strict         | Nearly identical sequences |
| 0.8      | 0.8      | Strict              | Highly conserved families  |
| 0.6      | 0.8      | Moderate            | Homologous families        |
| 0.5      | 0.7      | Permissive          | Distant homologs           |
| 0.3      | 0.5      | Very permissive     | Potential remote homology  |

### Coverage Mode Selection

| Mode | Coverage Target  | Best Use Case                                |
|------|------------------|----------------------------------------------|
| 0    | Query coverage   | Find complete matches to smaller sequences   |
| 1    | Target coverage  | Find sequences that completely match targets |
| 2    | Shorter sequence | Balanced approach (recommended)              |
| 3    | Longer sequence  | Conservative clustering                      |

## Output Files

PANORAMA generates clustering results with detailed cluster membership information.

### File Organization

```
output_directory/
└── clustering.tsv
```

### Clustering Results Format

The clustering results file contains the following columns:

| Column     | Description                         | Example     |
|------------|-------------------------------------|-------------|
| cluster_id | Unique cluster identifier (integer) | 0           |
| referent   | Representative gene family ID       | PG1_FAM_001 |
| in_clust   | Gene family member of this cluster  | PG2_FAM_045 |

### Results Interpretation

- **cluster_id**: Sequential integer identifying each cluster group
- **referent**: The representative (typically first or most abundant) family in the cluster
- **in_clust**: All gene families assigned to this cluster (including the referent)

Example clustering results:

```
cluster_id  referent      in_clust
0           PG1_FAM_001   PG1_FAM_001
0           PG1_FAM_001   PG2_FAM_045
0           PG1_FAM_001   PG3_FAM_078
1           PG1_FAM_002   PG1_FAM_002
1           PG1_FAM_002   PG2_FAM_123
```
