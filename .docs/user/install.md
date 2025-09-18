# PANORAMA – Installation Guide 📦⚒

## Latest version
### Installation via Conda 🐍

```{note}
PANORAMA is not yet available on Bioconda. We hope to provide a recipe soon.
```

```{important}
Before starting, create and activate a Conda environment. 
For more information, see the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/index.html).
```

* * *

### Installing from source code (GitHub) 🐙
(with_conda_env)=
#### Within a Conda environmnent 🐍

##### 1. Clone the Repository

```shell
git clone https://github.com/labgem/PANORAMA.git
cd PANORAMA
```

##### 2. Create and Configure the Conda Environment

```shell
conda create -n panorama
conda config --add channels bioconda
conda config --add channels conda-forge
conda activate panorama
conda env update --file panorama.yml
```

Alternatively, in one line:  
`conda create -n panorama -c bioconda -c conda-forge --file panorama.yml`

##### 3. Install PANORAMA

```shell
pip install .
```


##### 4. Test the Installation

```shell
panorama --help
panorama --version
```

#### 🛠️ Manual Installation (without Conda)

##### 1. Clone the Repository

```shell
git clone https://github.com/labgem/PANORAMA.git
cd PANORAMA
```

##### 2. Install Dependencies

Install **MMSeqs2 (≥ 13.45111)** manually:

* With [Homebrew](https://github.com/Homebrew/brew): `brew install mmseqs2`
* On Debian-based systems: `sudo apt install mmseqs2`

```{note}
Here is the complete installation guide of MMSeqs2: 
[https://github.com/soedinglab/MMseqs2/wiki#installation](https://github.com/soedinglab/MMseqs2/wiki#installation)
```

##### 3. Install PANORAMA and Python Dependencies

Dependencies are declared in the `pyproject.toml` file:

```shell
pip install .
```

* * *

## 🧪 Development Version

### 1. Get the Development branch 

To clone the `dev` Branch

```shell
git clone --branch dev https://github.com/labgem/PANORAMA.git
```


### 2. Install Dependencies

Follow the same steps as in the standard installation. 
With a [conda environment](#with_conda_env) or with a [manual installation](#manual-installation-without-conda).

### 3. Install the Development Version
```shell
pip install .
```
    
***

## Troubleshooting

If you encounter any problems, please check the known issues below.
If you still can't solve your problem, please [open an issue on GitHub](https://github.com/labgem/PANORAMA/issues).

**No issues yet**
