---
myst:
  html_meta:
    "description lang=en": |
      Step-by-step guide for adding a new subcommand to the PANORAMA CLI.
    "keywords": |
      Pangenomics, Python, Developer, CLI, Subcommand, argparse
---

# Adding a Subcommand

Every PANORAMA feature is exposed as an argparse subcommand. Adding a new one requires
three steps: implement the module, wire it into {py:mod}`panorama.main`, and register
the common arguments.

---

## 1. Implement the module

Create a new file (e.g. `panorama/myfeature/myfeature.py`). Every subcommand module
must expose exactly two callables at the module level:

### `subparser(subparsers)`

Registers the subcommand and its arguments. Return the created parser so `main.py` can
attach common arguments to it.

```python
import argparse

def subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "my_feature",
        description="One-line description shown in --help.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_my_feature(parser)
    return parser


def parser_my_feature(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-p", "--pangenomes",
        required=True,
        type=Path,
        help="TSV file listing pangenome .h5 files.",
    )
    required.add_argument(
        "-o", "--output",
        required=True,
        type=Path,
        help="Output directory.",
    )
```

Separating the parser construction into a `parser_my_feature` helper keeps the file
consistent with the rest of the codebase and makes unit-testing the argument logic easier.

### `launch(args)`

Called by `main()` with the parsed `argparse.Namespace`. This is where the feature runs.

```python
from panorama.format.read_binaries import load_pangenomes
from panorama.utils import mkdir, set_verbosity_level

def launch(args: argparse.Namespace) -> None:
    # 1. validate args and build need_info
    need_info = check_my_feature_args(args)

    # 2. load pangenomes
    pangenomes = load_pangenomes(
        pangenomes_file=args.pangenomes,
        need_info=need_info,
        disable_bar=args.disable_prog_bar,
        cpu=args.cpus,
    )

    # 3. create output directory
    mkdir(args.output, force=args.force)

    # 4. run the feature
    my_feature(pangenomes, args.output)
```

#### Validating arguments

Write a `check_my_feature_args` function that validates user input and returns the
`need_info` dict controlling which data groups are loaded from each `.h5` file:

```python
from typing import Dict

def check_my_feature_args(args: argparse.Namespace) -> Dict:
    return {
        "need_annotations": True,
        "need_families": True,
        "need_families_info": False,
        "need_rgp": False,
        "need_spots": False,
        "need_systems": False,
    }
```

Only set keys to `True` for data your feature actually uses — unnecessary loads slow
down startup on large pangenomes.

---

## 2. Register common arguments

{py:func}`panorama.utils.add_common_arguments` injects the following flags into every
subparser automatically (called by `main.py` for every registered subparser):

| Flag | Type | Description |
|---|---|---|
| `--verbose` / `-v` | int | Logging verbosity level |
| `--log` | Path | Write log to file |
| `-d` / `--disable_prog_bar` | flag | Disable tqdm progress bars |
| `--force` | flag | Overwrite existing outputs |

Do **not** declare these flags yourself — they are added by `main.py` after you return
your parser from `subparser()`.

---

## 3. Wire into `main.py`

Open {py:mod}`panorama.main` and make three additions:

### Import the module

```python
from panorama.myfeature.myfeature import launch as my_feature_launcher
from panorama.myfeature.myfeature import subparser as my_feature_subparser
```

### Register the subparser

Add your subparser to the `subs` list inside `cmd_line()`:

```python
subs = [
    info_subparser(subparsers),
    annotate_subparser(subparsers),
    ...
    my_feature_subparser(subparsers),   # ← add here
]
```

### Dispatch in `main()`

Add an `elif` branch to the dispatch chain:

```python
elif args.subcommand == "my_feature":
    my_feature_launcher(args)
```

Also add a one-line description to the `desc` string at the top of `cmd_line()` so your
subcommand appears in `panorama --help`.

---

## 4. Logging convention

Use the single named logger — do **not** create a per-module logger:

```python
import logging
logger = logging.getLogger("PANORAMA")
```

Call {py:func}`panorama.utils.set_verbosity_level` only once in `main()` before
dispatch; never call it from `launch()`.

---

## Checklist

- [ ] `subparser()` returns the parser
- [ ] `parser_my_feature()` separates argument definitions from registration
- [ ] `launch()` calls `check_my_feature_args()` to validate and build `need_info`
- [ ] `need_info` requests only the data groups actually needed
- [ ] Common arguments (`--force`, `--verbose`, etc.) are **not** declared manually
- [ ] Logger is `logging.getLogger("PANORAMA")`
- [ ] Import and register in `panorama/main.py` (imports + `subs` list + `elif` dispatch)
- [ ] One-line description added to `cmd_line()` help text
- [ ] Unit test for `check_my_feature_args()` and at least one functional test via `run_command()`
