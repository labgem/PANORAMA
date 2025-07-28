# 🧬 PANORAMA System Models Overview

PANORAMA detects macromolecular systems in pangenomes using user-defined models. These models are an exhaustive and
specific representation of a system. They describe the presence/abscence rules of protein families constituting a
system, how they are organized in genomes, and what genomic constraints govern their presence.

System models are written in JSON and provide a flexible, hierarchical structure that captures both essential and
optional elements of complex systems.

---

## 🧭 Definition and Purpose

A **system model** in PANORAMA serves two main purposes:

1. **Represent the system in a unique way**  
   Models specify which protein families are required, optional, or excluded for a given system.

2. **Guide the detection process**  
   The model governs how PANORAMA searches for gene co-occurrence and proximity in the pangenome context graph.

Models allow users to capture:

- Core and variable components of a system,
- Modular architecture (via *Functional Units*),
- Flexible quorum rules,
- Genome organization constraints (e.g., clustering, strand consistency),
- Inhibitory or exchangeable elements.

---

## 🧱 Structure (Brief)

Each model file is a single JSON object composed of:

- **`name`**: The system’s identifier.
- **`parameters`**: Detection parameters at the model level (e.g. quorum, transitivity).
- **`func_units`**: A list of *Functional Units* (modules), each of which groups multiple protein families.

The structure is hierarchical:

```
Model
├── FunctionalUnit
│ ├── Family
│ └── ...
├── FunctionalUnit
│ ├── Family
│ └── ...
```

Each **Functional Unit** and **Family** has a `presence` type (`mandatory`, `accessory`, `neutral`, `forbidden`) that
governs how it contributes to system detection.

Detection rules are defined at **both the model level** and the **unit level**, using parameters such as:

- `min_mandatory`
- `min_total`
- `transitivity`
- `same_strand`

For full details, see:

- 📄 [Model structure](./modeling.md#-model-structure)
- 🧩 [Presence Types](./modeling.md#-presence-types-explained)
- 🧰 [Detection rules](./modeling.md#-detection-parameters)

---

## PANORAMA models gallery

All existing PANORAMA models are available on the [PANORAMA models repository](https://github.com/PANORAMA-models).

Most of them are directly translated from different sources.

We welcome contributions and would like to provide PANORAMA users with a variety of models for a multitude of analyses.
You can contribute to the creation of new models by following the
guide [here](./contribute.md#how-to-contribute-to-panorama-models).

## Translate models

PANORAMA can translate several sources (see [PANORAMA models repository](https://github.com/PANORAMA-models)).

We ary trying to keep our models updated. But in case we miss it,
you can translate the model yourself with `panorama utils --translate` command.
Look here for more information: [PANORAMA – Model Translation Guide](./translate.md) 