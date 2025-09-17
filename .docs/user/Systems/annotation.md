##  🔢 Gene Family annotation

The `annotation` command adds **functional annotations** to gene families in pangenomes.
You can choose between:

- 📄 A **TSV file** with metadata
- 🧬 A **HMM database**, searched with `pyhmmer`

---

### ⚙️ Annotation Modes

#### 📄 1. TSV-based annotation

This mode injects gene family metadata from a `.tsv` file.

**Expected format**:
A TSV file where each row lists a pangenome name and the path to its gene family annotation file.

These annotation files contain functional details (e.g., protein name, accession, score, etc.).
The only mandatory column is `families`, which correspond to the gene families identifier.
See [metadata format](https://ppanggolin.readthedocs.io/en/latest/user/metadata.html#metadata-format) PPanGGOLiN documentation,
for more information.

#### 🧬 2. HMM-based annotation

To annotate with a HMM database, you must provide a HMM metadata file (TSV format), containing:

| Column               | Description                                                            | Type   | Mandatory |
|----------------------|------------------------------------------------------------------------|--------|-----------|
| name                 | The name of the HMM                                                    | string | True      |     
| accession            | Identifier of the HMM                                                  | string | True      |
| path                 | Path to the HMM file                                                   | string | True      |    
| length               | Length of the profile. Automatically recover by pyhmmer if necessary   | int    | False     |
| protein_name         | Name of the protein/function corresponding to the HMM                  | string | True      |
| secondary_name       | Secondary name of the protein                                          | string | True      |
| score_threshold      | Threshold used on the score to filter the profile                      | float  | False     |
| eval_threshold       | Threshold used on the E-value to filter the profile                    | float  | False     |
| ieval_threshold      | Threshold used on the iE-value to filter the profile                   | float  | False     |
| hmm_cov_threshold    | Threshold used on the HMM covering to filter the profile               | float  | False     |
| target_cov_threshold | Threshold used on the target covering to filter the profile            | float  | False     |
| description          | Description of the HMM, its protein function, or any other information | float  | False     |

:::{warning}
Not all the columns need to be filled with value as indicated by the mandatory column, but they should exist in the metadata file.
:::
:::{tips}
To keep all assignations possible between a profile and a gen family, you can let the threshold columns empty.
:::
:::{note}
You can generate the input files expected by PANORAMA using `panorama utils --hmm`.

[//]: # (TODO Ajouter le lien vers la documentation quand écrit)
:::

To align gene families against a HMM database, you can use different modes:

| Mode        | Description                                                |
|-------------|------------------------------------------------------------|
| `fast`      | Aligns representative sequences of each family to the HMMs |
| `profile`   | Builds HMMs for each family from MSAs                      |
| `sensitive` | Aligns **all genes** from each family to the HMMs          |

---

### 🚀 Command Line Usage

To annotate gene families with precomputed metadata, do as such:

```bash
panorama annotation \
  --pangenomes pangenomes.tsv \
  --source KEGG \
  --table annotations.tsv
  --threads 8
```

To annotate with a HMM database, do as such:

```bash
panorama annotation \
  --pangenomes pangenomes.tsv \
  --source defensefinder \
  --hmm hmms.tsv \
  --mode sensitive \
  --k_best_hit 3 \    # <-- or use the alias -b to keep only the best hit
  --output results/ \
  --threads 8
```
:::{tips}
More options are available to annotate with a HMM database. See below.
::::
:::{note}
source name should not contain a special character. They could interfere with the h5 writing.
:::
:::{warning}
You **must provide either** `--table` **or** `--hmm`, **but not both**. These options are mutually exclusive.
:::

#### 🔑 Key options

| Shortcut | Argument             | Description                                                                           |
|----------|----------------------|---------------------------------------------------------------------------------------|
| `-p`     | `--pangenomes`       | TSV file listing `.h5` pangenomes                                                     |
| `-s`     | `--source`           | Name of the annotation source (e.g. `KO2024`, `Pfam`)                                 |
| —        | `--table`            | **Mutually exclusive with `--hmm`**. TSV linking pangenome names to annotation files  |
| —        | `--hmm`              | **Mutually exclusive with `--table`**. HMM metadata TSV (from `panorama utils --hmm`) |
| —        | `--mode`             | Required with `--hmm`. Alignment strategy: `fast`, `profile`, or `sensitive`          |
| —        | `--msa`              | (Used only in `profile` mode) TSV listing MSAs per gene family                        |
| `-b`     | `--only_best_hit`    | Equivalent to `--k_best_hit 1`                                                        |
| —        | `--k_best_hit`       | Keep up to `k` best hits per gene family                                              |
| —        | `--output`           | Output directory for HMM result files (optional, used with `--save_hits`)             |
| —        | `--save_hits`        | Save HMM alignment results in formats: `tblout`, `domtblout`, `pfamtblout`            |
| —        | `--tmp`              | Temporary directory (used with HMM mode)                                              |
| —        | `--keep_tmp`         | Keep temporary files after HMM alignment                                              |
| —        | `--Z`                | Custom Z value for e-value scaling (advanced HMMER option)                            |
| —        | `--msa-format`       | Format of MSA files (default: `afa`) — rarely changed                                 |
| —        | `--disable_prog_bar` | Disable progress bars in multi-threaded execution                                     |
| —        | `--force`            | Overwrite existing annotation in the pangenome                                        |
| —        | `--threads`          | Number of threads to use                                                              |


### 🧪 Annotation Workflow
1. 🧬 Load pangenomes

    Pangenomes are loaded from .h5 files. Only necessary information is retrieved based on the mode.

2. 📂 Retrieve annotations

   - With `--table`: loads metadata from TSV

   - With `--hmm`: aligns families via annot_with_hmm() from hmm_search.py

3. (only for hmm option) 🎯 Filter HMM hits

    Each hit is filtered using the thresholds defined in the HMM metadata:
    - e-value
   - i-evalue
   - score
   - target coverage
   - HMM coverage

:::{tips}
Prefer to use the score instead of the e-value or the i-evalue to ensure reproducibility of the results even if the size of your targets changes.
:::

4. 💾 Write annotations

   Filtered annotations are stored in the .h5 files, under the given --source name.

:::{note}
Annotations can be viewed or reused with panorama, ppanggolin, or custom tools (*e.g.*, [vitables](https://vitables.org/index.html)).
:::

### 🧠 HMM Search Details

Annotation relies on the [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/index.html) Python API.

Depending on sequence size, PANORAMA chooses the best method:

| Method    | Use case                               |
|-----------|----------------------------------------|
| hmmsearch | In-memory, fast                        |
| hmmscan   | Streaming, used when memory is limited |

:::{note}
If sequences exceed 10% of available RAM, PANORAMA uses hmmscan, as recommended by pyhmmer documentation 
[here](https://pyhmmer.readthedocs.io/en/stable/examples/performance_tips.html#Performance-tips-and-tricks)
:::

### 📦 Minimal example

#### Annotate gene families based on the reference sequence with COG HMM

```bash
panorama annotation \
  -p pangenomes.tsv \
  -s COG \
  --hmm hmms.tsv \
  --mode fast \
  --only_best_hit   # <-- or use the alias: -b
```
