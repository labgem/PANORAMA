from __future__ import annotations

import argparse
import json

from pathlib import Path
from rules import System, Systems


def read_files(systems_path, systems=Systems()):
    """
    Read all json files in the directory

    :param systems_path: path of systems directory
    :param systems: class Systems with all systems
    """
    for file in systems_path.glob("*.json"):
        with open(file.resolve().as_posix()) as json_file:
            data = json.load(json_file)
            system = System()
            system.read_system(data)
            systems.add_sys(system)
    systems.print_systems()


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    systems_path = Path(args.systems)
    read_files(systems_path)


def parser_annot(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-s', '--systems',
                          required=True,
                          type=str,
                          help="Path to systems directory")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser_annot(parser)
    args = parser.parse_args()
    launch(args)
