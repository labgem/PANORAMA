#!/usr/bin/env python3
"""
Checks that duplicated CLI argument tables in the user documentation stay in sync.

`panorama/utils.py:add_common_arguments` injects --verbose, --log, --disable_prog_bar,
and --force into every subcommand, and several docs/user/*.md pages repeat the same
"--pangenomes/--output/--models/--sources" required-argument rows. Rather than a single
included partial (which breaks Markdown table rendering when merged with page-specific
rows), each page keeps its own literal copy. This script diffs those copies so they
cannot silently drift apart.

Exits non-zero and prints a diff-like report if any duplicated row differs between files.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

DOCS_USER = Path(__file__).resolve().parents[2] / "docs" / "user"

# Files whose "Required Arguments" table is expected to have identical
# --pangenomes / --output / --models / --sources rows.
REQUIRED_ARGS_FILES = ["projection.md", "association.md", "partition.md"]
REQUIRED_ARGS = ["--pangenomes", "--output", "--models", "--sources"]

# Files whose "common optional arguments" rows (injected by add_common_arguments) use the
# same 4-column (Argument/Type/Default/Description) table shape and are expected to be
# identical. compare_spots.md and compare_systems.md document the same flags but with a
# different 5-column table shape (they add a "Shortcut" column), so they are intentionally
# left out of this row-for-row comparison.
COMMON_OPTIONAL_ARGS_FILES = [
    "projection.md",
    "association.md",
    "partition.md",
    "pansystems.md",
]
COMMON_OPTIONAL_ARGS = ["--verbose", "--log", "--disable_prog_bar", "--force"]


def parse_table_row_cells(line: str) -> list[str]:
    """Split a Markdown table row into its stripped cell contents."""
    cells = line.strip().split("|")
    # A well-formed "| a | b | c |" row splits into ['', ' a ', ' b ', ' c ', ''].
    return [c.strip() for c in cells[1:-1]]


def find_row_for_argument(lines: list[str], argument: str) -> tuple[int, list[str]] | None:
    """Find the first table row mentioning `argument` (e.g. as `--force`), return (line_no, cells)."""
    pattern = re.compile(rf"`?{re.escape(argument)}`?")
    for i, line in enumerate(lines, start=1):
        if line.strip().startswith("|") and pattern.search(line):
            return i, parse_table_row_cells(line)
    return None


def check_argument_group(files: list[str], arguments: list[str]) -> list[str]:
    errors = []
    for argument in arguments:
        rows: dict[str, tuple[int, list[str]]] = {}
        for filename in files:
            path = DOCS_USER / filename
            if not path.exists():
                errors.append(f"{filename}: file not found")
                continue
            lines = path.read_text().splitlines()
            found = find_row_for_argument(lines, argument)
            if found is None:
                errors.append(f"{filename}: no row found for {argument}")
                continue
            rows[filename] = found

        if not rows:
            continue

        reference_filename, (_, reference_cells) = next(iter(rows.items()))
        for filename, (line_no, cells) in rows.items():
            if cells != reference_cells:
                errors.append(
                    f"{argument} row differs between {reference_filename} and {filename} "
                    f"(docs/user/{filename}:{line_no}):\n"
                    f"    {reference_filename}: {reference_cells}\n"
                    f"    {filename}: {cells}"
                )
    return errors


def main() -> int:
    errors = []
    errors += check_argument_group(REQUIRED_ARGS_FILES, REQUIRED_ARGS)
    errors += check_argument_group(COMMON_OPTIONAL_ARGS_FILES, COMMON_OPTIONAL_ARGS)

    if errors:
        print("Common argument tables have drifted out of sync:\n")
        for error in errors:
            print(f"- {error}")
        print(
            "\nUpdate the affected docs/user/*.md files so the duplicated rows match "
            "again, or update this script's expectations in "
            ".github/scripts/check_common_arguments_docs.py if the CLI itself changed."
        )
        return 1

    print("Common argument tables are consistent across all checked pages.")
    return 0


if __name__ == "__main__":
    sys.exit(main())