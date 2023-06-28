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

```{note}
You can find more information on the installation in [install section](#install-section)
```


