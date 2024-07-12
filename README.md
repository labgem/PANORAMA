# PANORAMA: Pangenomics toolbox for comparison and analyses of prokaryotic pangenomes

PANORAMA is a software suite used to analyse and compare partitioned pangenomes graph provided. It benefits from 
methods for the reconstruction and analysis of pangenome graphs, thanks to the [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN)
software suite. It is designed to perform pangenome comparison at high-throughtup level.

# Quick Installation
PANORAMA is easily installed with [conda](https://docs.conda.io/projects/conda/en/latest/index.html) and
[pip](https://pip.pypa.io/en/stable/). Follow the next step to install panorama.

```shell
# Configuration of conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Requirements installation
conda install --file requirements.txt

# PANORAMA installation
pip install .
```

[//]: # (You can find more information on the installation in mettre le lien vers ReadTheDocs)

# PANORAMA overview
## Input files
PANORAMA can process multiple pangenome in one time. The common input file of most of the commands is a *TSV* file with 
2 columns. In the following we will name this file *organisms_pangenome.list*

| Name       |               Path |
|:-----------|-------------------:|
| Pangenome1 | path/to/pangenome1 |
| ...        |                ... |
| PangenomeX | path/to/pangenomeX |

*NB: We recommand to use an absolute path in this file to avoid errors. 
You can use the path from your current directory or the path from the input file as relative path to find pangenomes*

## Biological systems detection
PANORAMA allow to detect systems in pangenomes by using models. A model is an exhaustive and specific representation of
a system. PANORAMA models are flexible to describe any models provide by user. 

First of all you need to annotate the gene families. For this you can use a list of TSV (like for pangenomes) 
or a HMM database. Here we are going to use HMM contain in hmm_directory.

```shell
panorama annotation --pangenomes organisms_pangenome.list --hmm hmm_directory --source source
```

To annotate the pangenome you must provide an annotation source. This will be used as key to save annotations in the
pangenome. In the next step to detect systems you must provide the same annotation source to use those annotations.

It's now possible to detect systems with PANORAMA.

```shell
panorama systems --pangenomes organisms_pangenome.list --models models_directory --source source
```
System are now saved in the pangenome. 


# Issues, Questions, Remarks

If you have any question or issue, please do not hesitate to post an issue ! 
We cannot correct bugs if we do not know about them, and will try to help you the best we can.




