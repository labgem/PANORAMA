# PANORAMA: A robust pangenome-based method for predicting and comparing biological systems across species

[![Actions](https://img.shields.io/github/actions/workflow/status/labgem/PANORAMA/main.yml?branch=dev&event=pull_request&label=build&logo=github)](https://github.com/labgem/PANORAMA/actions/workflows/main.yml)
[![License](https://anaconda.org/bioconda/ppanggolin/badges/license.svg)](http://www.cecill.info/licences.fr.html)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/panorama?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/panorama)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/labgem/PANORAMA/)
[![GitHub issues](https://img.shields.io/github/issues/labgem/panorama.svg?style=flat-square&maxAge=600)](https://github.com/labgem/panorama/issues)
[![Docs](https://img.shields.io/readthedocs/panorama/latest?style=flat-square&maxAge=600)](https://panorama.readthedocs.io)
[![Downloads](https://anaconda.org/bioconda/panorama/badges/downloads.svg)](https://bioconda.github.io/recipes/panorama/README.html#download-stats)

PANORAMA is a software suite used to analyse and compare partitioned pangenomes graph provided. It benefits from
methods for the reconstruction and analysis of pangenome graphs, thanks to
the [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN)
software suite. It is designed to perform pangenome comparison at high-throughtup level.

# Quick Installation

PANORAMA is easily installed with [conda](https://docs.conda.io/projects/conda/en/latest/index.html) and
[pip](https://pip.pypa.io/en/stable/). Follow the next step to install panorama.

```shell
# 1. Clone the Repository
git clone https://github.com/labgem/PANORAMA.git
cd PANORAMA

# 2. Create and Configure the Conda Environment
conda create -n panorama
conda config --add channels bioconda
conda config --add channels conda-forge
conda activate panorama
conda env update --file panorama.yml
```

[//]: # (You can find more information on the installation [here]&#40;link_read_the_doc&#41;)

# PANORAMA overview

## Input files

PANORAMA can process multiple pangenomes at a time.
The common input file of most of the commands is a *TSV* file with 2 columns.
In the following we will name this file *pangenomes.tsv*

| Name       | Path               |
|------------|--------------------|
| Pangenome1 | path/to/pangenome1 |
| ...        | ...                |
| PangenomeX | path/to/pangenomeX |

*NB: We recommand to use an absolute path in this file to avoid errors.
You can use the path from your current directory or the path from the input file as relative path to find pangenomes*

## Biological systems detection

PANORAMA allow to detect systems in pangenomes by using models.
A model is an exhaustive and specific representation of a system.
PANORAMA models are flexible to describe any models provide by user.
PANORAMA provide a command to perform the complete detection workflow as follow:

```shell
panorama pansystems \
-p pangenomes.tsv \
--hmm /PATH/TO/HMM/LIST/FILE/hmm_list.tsv \
 -m /PATH/TO/MODELS/LIST/FILE/models_list.tsv \
 -s system_model_source_name \
-o PATH/TO/OUPUT/DIRECTORY \
--projection
```

## Pangenome comparison

PANORAMA allow to compare pangenomes based on pangenome structure previously detected.

### Pangenome comparison based on spots

Basic conserved spots comparison:

```shell
panorama compare_spots \
--pangenomes pangenomes.tsv \
--output conserved_spots_results \
--gfrr_metrics min_gfrr \
--gfrr_cutoff 0.8 0.8 \
--threads 8
```

### Pangenome comparison based on systems

Basic systems comparison with heatmap generation:

```shell
panorama compare_systems \
--pangenomes pangenomes.tsv \
--models defense_systems.tsv \
--sources defense_finder \
--output systems_comparison_results \
--heatmap \
--threads 8
```

## Utilities

### Alignment & clustering of gene families

To perform the comparison of pangenomes, gene families are aligned and cluster.
PANORAMA provide commands to perform the alignment and clustering of gene families before the comparison of pangenomes.

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

Fast clustering with linclust:

```shell
panorama cluster \
--pangenomes pangenomes.tsv \
--output clustering_results \
--method linclust \
--cluster_identity 0.8 \
--cluster_coverage 0.8 \
--threads 8
```

## 💬 Feedback & Contribution

**Give us feedback**  
PANORAMA is still in early development.

Have suggestions, ideas, or bug reports?  
👉 [Open an issue on GitHub](https://github.com/labgem/PANORAMA/issues)

We cannot correct bugs if we do not know about them, and will try to help you the best we can.
