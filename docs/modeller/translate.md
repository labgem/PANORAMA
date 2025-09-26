# 🧬 PANORAMA – Model Translation Guide

The model translation module in PANORAMA allows you to convert system detection models from other tools into the
PANORAMA format.
Also it will prepare all required input files, so they can be used directly for annotation and analysis in PANORAMA.

**Supported Sources**

- PADLOC (e.g., antiviral defense systems),
- MacSyFinder-based tools like:
    - DefenseFinder (e.g., restriction-modification, retrons, etc.),
    - CASFINDER ()
    - CONJScan (conjugative systems),
    - TXSScan (type secretion systems),
    - TFFScan (type IV pili),

## 🎯 When should you use it?

```{important}
Translated and updated versions of models are already available [here](https://github.com/PANORAMA-models).
```

1. When a new version exists, and it's not yet available on [PANORAMA models](https://github.com/PANORAMA-models),
2. When you want to use a specif anterior version not available on [PANORAMA models](https://github.com/PANORAMA-models)
3. When you have your own version of models that use the same grammar than supported sources.

## 📂 What do you need?

You should have:

| Tool/Source          | HMM profiles                        | models                        | metadata                               |
|----------------------|-------------------------------------|-------------------------------|----------------------------------------|
| PADLOC               | Hidden Markov Model profiles (.hmm) | Models in YAML format (.yaml) | File containing HMM information (.txt) |
| MacSyModel (or like) | Hidden Markov Model profiles (.hmm) | Models in XML format (.xml)   |                                        | 

## 🚀 How to run it?

Use the command line interface (CLI):

```shell
panorama utils \
--translate /path/to/models_directory \
--source name_of_the_source \
--output /path/to/output_directory \
--binary \
--hmm_coverage 0.6 \
--target_coverage 0.8
```

### 🔑 Key options

| Parameter           | Type     | Description                                                               |
|---------------------|----------|---------------------------------------------------------------------------|
| `--translate`       | Required | The folder with PADLOC, DefenseFinder or MacSyFinder models.              |
| `--source`          | Required | One of the following: padloc, defense-finder, CONJScan, TXSScan, TFFscan. |
| `--output`          | Required | Where you want PANORAMA to save its output.                               |
| `--binary`          | Optional | Use this to speed up annotations (recommended).                           |
| `--hmm_coverage`    | Optional | Thresholds for alignment coverage on HMM                                  |
| `--target_coverage` | Optional | Thresholds for alignment coverage on target.                              |

```{tip}
Defense Finder and PADLOC use different threshold for the HMM and the target coverage.
Defense Finder use 0.4 for HMM and 0 for target for all profiles.
PADLOC define in the metadata a specific threshold for each profile.
This default behavior is tunable with the `--hmm_coverage` and `--target_coverage` options 
that will affect a threshold for all profiles.
```

## 📁 What files are created?

After translation, you’ll get:

- **models/**: All models in PANORAMA .json format,
- **models_list.tsv**: A list of all translated models,
- **hmm/**: All HMM profile compatible with the hmm_list TSV file, binary or not,
- **hmm_list.tsv**: A list of all hmm with the metadata needed for the annotation.

These files can be used directly in PANORAMA commands.

## 🛠️ Example

### Translate PADLOC

```shell
git clone https://github.com/padlocbio/padloc-db.git
mkdir padloc_translate
panorama utils \
--translate padloc-db/ \
--source padloc \
--output padloc_translate/ \
--binary
```
```{tip}
PANORAMA automatically find all required files on its own, as the structure is already known.
```

### Translate DefenseFinder

```shell
git clone https://github.com/mdmparis/defense-finder-models.git
mkdir dfinder_translate
panorama utils \
--translate defense-finder-models/ \
--source defense-finder \
--output dfinder_translate/ \
--binary
```
```{note}
Defense Finder team as the tendency to reorganize the repository.
If you're trying to translate a new version this could not work.
In this case report an issue [here](https://github.com/PANORAMA-models/Defense-Finder)
```
