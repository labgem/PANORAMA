# 🧬 PANORAMA System Modeling

## 🧰 Definition

PANORAMA detects macromolecular systems in pangenomes using user-defined models. These models specify what protein
families constitute a system and how they are expected to be organized genomically.

While inspired by the [MacSyFinder](https://macsyfinder.readthedocs.io/en/latest/modeler_guide/modeling.html)
or [PADLOC](https://github.com/padlocbio/padloc) modeling format, PANORAMA introduces several innovations:

- Components are defined as protein families, not individual genes.

- Families are grouped into Functional Units, which capture biologically meaningful modules.

- A canonical model mechanism allows detection of incomplete or hypothetical systems.

:::{important}
PANORAMA system models are written in JSON format.
:::

## 🧱 Model Structure

A model is composed of:

- One or more **Functional Units**, each containing one or more **Families**
- A set of parameters to specify the *quorum* and the co-localization rules

### 🔎 Example

```json
{
  "name": "NanoDefense_V1",
  "parameters": {
    "transitivity": 4,
    "window": 5,
    "min_mandatory": 2,
    "min_total": 3
  },
  "func_units": [
    {
      "name": "DetectionUnit",
      "presence": "mandatory",
      "parameters": {
        "min_mandatory": 1,
        "min_total": 1
      },
      "families": [
        {
          "name": "ND-SensorA",
          "presence": "mandatory"
        },
        {
          "name": "ND-SensorB",
          "presence": "accessory"
        },
        {
          "name": "ND-Disruptor",
          "presence": "forbidden"
        }
      ]
    },
    {
      "name": "ResponseUnit",
      "presence": "mandatory",
      "parameters": {
        "min_mandatory": 1,
        "min_total": 2
      },
      "families": [
        {
          "name": "ND-Toxin1",
          "presence": "mandatory",
          "exchangeable": [
            "ND-Toxin2"
          ]
        },
        {
          "name": "ND-Delivery",
          "presence": "accessory"
        }
      ]
    },
    {
      "name": "ControlUnit",
      "presence": "accessory",
      "parameters": {
        "min_mandatory": 1,
        "min_total": 1
      },
      "families": [
        {
          "name": "ND-Regulator",
          "presence": "mandatory"
        }
      ]
    },
    {
      "name": "InsertUnit",
      "presence": "neutral",
      "families": [
        {
          "name": "NS-md",
          "presence": "mandatory"
        },
        {
          "name": "NV-ac",
          "presence": "accessory"
        }
      ]
    }
  ]
}
```

This fictional system represents a defense mechanism composed of sensor, effector, and regulatory units.
It models a modular system architecture using three functional units:

- Detection unit (mandatory): Detects environmental signals or threats.
    - ND-SensorA: mandatory — at least one detection family must be present.
    - ND-SensorB: accessory — may appear in some variants, enhancing specificity.
    - ND-Disruptor: forbidden — if present, disqualifies the system (may represent a mobile element or anti-system
      gene).
    - To validate this unit, one mandatory and a total of one family is required.
- Response unit (mandatory): Delivers a toxic response to eliminate the threat.
    - ND-Toxin1: mandatory, but exchangeable with ND-Toxin2 — either can fulfill the same role.
    - ND-Delivery: accessory — a delivery mechanism that might enhance toxin efficiency.
    - To validate this unit, one mandatory and two total families are required.
- Control unit (accessory): An optional regulatory unit that may fine-tune the response.
    - ND-Regulator: mandatory within the unit, but the whole unit is accessory.
    - To validate this unit, one mandatory and a total of one family is required.
- Insert unit (neutral): A neutral unit often found in the system context but without known function.
    - NS-md: mandatory within the unit, if the unit is found this family is too
    - NV-ac: accessory within the unit, not always present in the unit.

:::{seealso}
This example shows a fairly complete and specific model.
In the next section, we'll look at how to create more simplified models.
:::

### 🧩 Components

#### Model

#### Functional Units

Functional Units represent a set of genes/families that together perform a system function. Each has:

| Field      | Description                                    | Required/Optional | Possible Values                                                                   |
|------------|------------------------------------------------|-------------------|-----------------------------------------------------------------------------------|
| name       | Unique name identifying the functional unit    | 🔴 Required       | String (annotation identifier)                                                    |
| presence   | Role of the unit                               | 🔴 Required       | `mandatory`, `accessory`, `neutral`, `forbidden`                                  |
| families   | List of protein families included in this unit | 🔴 Required       | List of family objects                                                            |
| parameters | List of specific rules to detect the unit      | 🟡 Optional       | Dictionary with fields describe in [detection parameters](#-detection-parameters) |

A functional unit could biologically represent a functional module, such as isoenzyme, or subunit of protein dimers.

:::{attention}
A model must include at least one functional unit.
:::

#### Families

Families correspond to isofunctional protein families used to search the pangenome annotations. Each family has:

| Field        | Description                                               | Required/Optional | Possible Values                                  |
|--------------|-----------------------------------------------------------|-------------------|--------------------------------------------------|
| name         | Identifier matching the annotation (e.g. from HMM source) | 🔴 Required       | String (annotation identifier)                   |
| presence     | represent the contribution of the family                  | 🔴 Required       | `mandatory`, `accessory`, `neutral`, `forbidden` |
| duplicate    | Max number of times this family can occur                 | 🟡 Optional       | Integer (positive number)                        |
| exchangeable | List of other families that can substitute this one       | 🟡 Optional       | Array of strings (family names)                  |
| multi_system | Can be used in multiple predicted systems                 | 🟡 Optional       | Boolean (`true`/`false`)                         |
| multi_model  | Can be shared across models                               | 🟡 Optional       | Boolean (`true`/`false`)                         |

:::{attention}
A family must be included in a functional unit.
:::
:::{warning}
A family can theoretically be in multiple unit, but this feature has never been tested.
:::

### 🎯 Presence Types Explained

Each family and each functional unit in a PANORAMA model must be assigned a presence type.
This type determines how the element contributes to system detection and scoring.

The presence type influences:

- Whether the element must be present or absent,
- Whether it affects the detection quorum,
- Whether it can be used for connectivity or context analysis only.

Below is a complete reference:

| Presence Type | Applies To   | Required for Detection? | Affects Quorum? | May Link Components in Graph? | Notes                                                                                                       |
|---------------|--------------|-------------------------|-----------------|-------------------------------|-------------------------------------------------------------------------------------------------------------|
| mandatory     | Family, Unit | 🔴 Required             | ✔ Yes           | ✔ Yes                         | Must be present for the system or unit to be considered valid.                                              |
| accessory     | Family, Unit | 🟡 Optional             | ✔ Yes           | ✔ Yes                         | Helpful but not required; counted toward min_total. Often lost or divergent in evolution.                   |
| forbidden     | Family, Unit | 🚫 Must be absent       | ❌ No            | 🚫 No                         | If present, it disqualifies the system or unit. Useful to distinguish similar systems or detect inhibitors. |
| neutral       | Family, Unit | 〰 Ignored               | ❌ No            | ✔ Yes                         | Ignored for scoring, but included in the graph. Helps connect elements that are close in genomic context.   |

### ⚙️ Detection rules

Parameters are defined at the model or functional unit level, such as:

```json
"parameters": {
"min_mandatory": 2,
"min_total": 3,
"transitivity": 5,
"window": 6
"same_strand": false,
}
```

#### Quorum rules

The same model can represent systems with the same function, but a different composition.
Quorum parameters are used to define the quantity of elements required to guarantee a functional system.

| Parameter     | Description                                      | Required/Optional | Possible Values           |
|---------------|--------------------------------------------------|-------------------|---------------------------|
| min_mandatory | Minimum number of mandatory elements             | 🔴 Required       | Integer (positive number) |
| min_total     | Minimum number of mandatory + accessory elements | 🔴 Required       | Integer (positive number) |

#### Genomic organisation

Genomic organization parameters define the space in which families should be located and how far apart they should be.

| Parameter    | Description                                                          | Required/Optional | Possible Values           |
|--------------|----------------------------------------------------------------------|-------------------|---------------------------|
| transitivity | Max distance in the pangenome graph between gene families            | 🔴 Required       | Integer (positive number) |
| window       | Size of genomic window examined (default: transitivity + 1)          | 🟡 Optional       | Integer (positive number) |
| same_strand  | true if all families must be on the same DNA strand (default: False) | 🟡 Optional       | Boolean (`true`/`false`)  |

#### Parameter inheritance

All the parameters described before should be set at the model level.
All functional units will, by default, inherit the parameters from the model.
However, it's possible to redefine one or more parameters for a functional unit.

Example:

```json
{
  "parameters": {
    "min_mandatory": 2,
    "min_total": 3,
    "transitivity": 5
  },
  "func_units": [
    {
      "name": "FA",
      "parameters": {},
      "families": []
    },
    {
      "name": "FB",
      "parameters": {
        "min_mandatory": 1,
        "transitivity": 2
      },
      "families": []
    }
  ]
}
```

Here, the functional unit _FA_ inherits all the parameters from the model, whereas _FB_ redefines the _min_mandatory_
and the _transitivity_.

:::{attention}
Either it's possible to don't precise all parameters in functional unit, the `parameters` field must exist.
To let the functional unit inherits all parameter you can let the dictionary empty.
:::

### 🧪 Canonical Models

PANORAMA supports the concept of canonical models, which represent partial, hypothetical,
or computationally predicted systems.
They Are only used if no non-canonical model explains the matching families.
Canonical model can help detect system variants or candidates for novel systems.

To define a canonical model, include:

```json
"canonical": [
        non-canonical_model_A,
        non-canonical_model_B,
        non-canonical_model_C
]
```

## ✅ Validation Rules

Model files are automatically validated and loaded using the PANORAMA engine (models.py). Checks include:

- Required keys (name, func_units, etc.)
- Valid presence types (mandatory, etc.)
- Consistent quorum thresholds (min_mandatory ≤ total mandatory)
- At least one mandatory family in each unit and one mandatory unit in each model

Models failing these checks will raise clear exceptions.

## 📝 Notes

1. Each model must be saved in its own .json file.
2. Names are case-sensitive. 
3. Families must match the name given during the annotation step (see [annotation command](../user/annotation.md#-gene-family-annotation)). 
4. Exchangeable families inherit the parameters (presence, etc.) of their reference unless specified otherwise.