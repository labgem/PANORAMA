## 🔭 System Projection on Genomes

The `write_systems` command enables the projection of systems, previously detected at the pangenome level 
(see [`systems` command](./detection.md), onto individual genomes. 
Projection relies on system detection results and the genomic context of gene families within organisms.

### 🧪 Projection Workflow
The projection process proceeds as follows:

1. 📂 Load Detected Systems and Metadata

    Detected systems from the .h5 pangenome file are loaded, along with required metadata and gene families.

2. 🧱 Build systems organisms contexts graphs

    For each organism and functional unit:
    1. Build a gene context graph using adjacency in the genome (based on window size and strand orientation). 
    2. Nodes represent genes; edges represent local gene neighborhood.

3. 🧠 Project System Units

    Each system unit is evaluated in organisms:
    1. Genes belonging to system families are checked for co-localization and graph connectivity. 
    2. Subsystems are classified as:
       - strict: All families in the model are connected by a maximum distance equal to the *transitivity* parameter. 
       - extended: All model families are found in the same context but additional families between them. 
       - split: All model families are found in the organisms but in different connected components. 
          Such fragmentation can occur in a subset of genomes due to rearrangements, insertion events (e.g., insertion sequences), or assembly breaks.

4. 🔄 Aggregate Projections

    Results are collected at:
   - Organism-level: Individual gene projections with system metadata.
   - Pangenome-level: Summarized presence/absence and completeness across all organisms.

5. 📝 Write Output (more detail in [output section](#-output-files))

    Projection results are written as TSV files:
    - **systems.tsv**: One summary file for the pangenome
    - **[organism_name].tsv**: One detailed file per organism.

### 🚀 Command Line Usage

Projection is triggered via:

```shell
panorama write_systems \
--pangenomes pangenomes.tsv \
--models models.tsv \
--sources defense_finder \
--projection \
--threads 8
```

To project systems only to a specific set of organisms:

```shell
--organisms organism_A organism_B
```

#### ⚙️ Key Options for Projection
| Argument        | Description                                               |
|-----------------|-----------------------------------------------------------|
| `--projection`  | Enable the projection of systems onto genomes             |
| `--organisms`   | List of organisms to project systems to (defaults to all) |
| `--canonical`   | Also project canonical versions of systems                |
| `--association` | Add associations with RGPs, spots, or modules             |
| `--threads`     | Number of parallel threads to use                         |
| `--force`       | Overwrite existing projection files                       |

(output-projection)=
#### 📄 Output Files

Output is written inside the specified `--output` directory. 
For each pangenome, a directory based on the pangenome name is created to save the results.

##### 1. Pangenome systems summary

First a TSV file summarizing all systems detected at the pangenome level is written.


| Column               | Description                                                                                                                |
|----------------------|----------------------------------------------------------------------------------------------------------------------------|
| system number        | Unique ID for the system                                                                                                   |
| system name          | name of the system corresponding to the model name                                                                         |
| functional unit name | List of functional units that composed the system (from the model)                                                         |
| organisms            | List of organisms where the system is detected                                                                             |
| model_GF             | List of gene families that encoded for the system                                                                          |
| context_GF           | List of gene families that are found in the context of the system but did not encoded any function described in the model  |
| partition            | system partition conciliate with model gene families                                                                       |
| completeness         | Average of the number of model gene families found in organisms compared with the number of families making up the system. |
| strict               | Number of organism where the system is found in a strict genomic organisation                                              |
| extended             | Number of organism where the system is found in an extended genomic organisation                                           |
| split                | Number of organism where the system is found in a split genomic organisation                                               |

:::{note}
If you use the `--association` option of the command line. Additional columns will be added, 
corresponding to the system association with RGPs, spots, or modules.   
:::

##### 2. Organism projection file

projection/<organism>.tsv: One file per organism with all matching genes and metadata.

Each organism file includes:

| Column               | Description                                  |
|----------------------|----------------------------------------------|
| system number        | Unique ID for the system                     |
| subsystem number     | ID for the subgraph (connected component)    |
| functional unit name | Name of the functional unit (from the model) |
| gene family          | Projected gene family                        |
| annotation           | Matched annotation (from metadata)           |
| genomic organization | strict, extended, or split                   |
| ...                  | Genomic location and feature annotations     |



🧠 Notes
Systems are projected independently per unit; partial matches are retained and classified.

Context-aware graph-building ensures precision in localizing system units.

Projection supports RGPs, spots, and modules if enabled via --association.