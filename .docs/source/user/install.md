(install-section)=
# Installation

## Conda install
```{note}
For the moment, there is no recipe available on bioconda. We hope we can add one soon
```

First of all you need to create or activate your conda environment. For more information you can look at [conda documentation](https://docs.conda.io/projects/conda/en/latest/index.html)

Your environment should have a python version greater than 3.8 (included)

The first step is to add the required channels in your environment. 
```shell
# Configuration of conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

You need to clone [PANORAMA GitHub repository](https://github.com/labgem/PANORAMA) localy by running the next line
```shell
git clone https://github.com/labgem/PANORAMA.git 
```
You can now install the required packages with the as follows:
```console
# Requirements installation
cd PANORAMA
conda install --file requirements.txt
```

You will need now to use pip to install PANORAMA. There is a *pyproject.toml* file to setup the installation.

```shell
# in the git repository 
pip install .
```

Congratulation, PANORAMA is installed !!!

```{admonition} Give us feed back
We are at the beginning of our project development. For the time we will add PANORAMA to bioconda.
If you have any suggestion for enhancement or idear for another package repository where PANORAMA could be added
post an issue on our GitHub.  
```