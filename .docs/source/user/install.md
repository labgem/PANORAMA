# 🧬 PANORAMA – Installation Guide
## 🐍 Installation via Conda

:::{note}
PANORAMA is not yet available on Bioconda. We hope to provide a recipe soon.
:::

> 💡 _Before starting, create and activate a Conda environment. For more information, see the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/index.html)._

* * *

## 🚀 Installing from source code (GitHub)
### 🐍 Within a Conda environmnent

1. Clone the Repository

```shell
git clone https://github.com/labgem/PANORAMA.git
cd PANORAMA
```

2. Create and Configure the Conda Environment
```shell
conda create -n panorama
conda config --add channels bioconda
conda config --add channels conda-forge
conda activate panorama
conda env update --file panorama.yml
```

Alternatively, in one line:  
`conda create -n panorama -c bioconda -c conda-forge --file panorama.yml`

3. Install PANORAMA
```shell
pip install .
```
    

4. Test the Installation
```shell
panorama --help
panorama --version
```

### 🛠️ Manual Installation (without Conda)

1. Clone the Repository

```shell
git clone https://github.com/labgem/PANORAMA.git
cd PANORAMA
```

2. Install Dependencies

Install **MMSeqs2 (≥ 13.45111)** manually:

* On Debian-based systems: `sudo apt install mmseqs2`
    

3. Install PANORAMA and Python Dependencies

Dependencies are declared in the `pyproject.toml` file:

```shell
pip install .
```
    
* * *

## 🧪 Development Version
1. Clone the `dev` Branch

```shell
git clone --branch dev https://github.com/labgem/PANORAMA.git
```
    

2. Install Dependencies

Follow the same steps as in the standard installation.

3. Install the Development Version
```shell
pip install .
```
    
* * *

## 💬 Feedback & Contribution

**Give us feedback**  
PANORAMA is still in early development. We plan to submit it to Bioconda soon.  
  
Have suggestions, ideas, or bug reports?  
👉 [Open an issue on GitHub](https://github.com/labgem/PANORAMA/issues)

