---
myst:
  html_meta:
    "description lang=en": |
      Documentation for users who wish to use PANORAMA
---

# User Documentation

PANORAMA is organised around three main workflows: **biological system prediction** (annotate gene families, detect
systems, write results), **pangenome comparison** (compare conserved spots and detected systems across pangenomes),
and a set of **utilities** (pangenome info, alignment, clustering) that support both. If you are new to PANORAMA,
start with the installation guide and the quick usage page below, then jump into the workflow that matches your goal.

```{toctree}
:caption: 'Get Started:'
:maxdepth: 2

install
quick_usage
issues
citation
```

```{toctree}
:caption: 'Biological system Prediction:'
:maxdepth: 2

annotation
detection
write_systems
```

```{toctree}
:caption: 'Pangenomes comparison:'
:maxdepth: 2

compare_spots
compare_systems
```

```{toctree}
:caption: 'PANORAMA utilities:'
:maxdepth: 2

info
align
clustering
```