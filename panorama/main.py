#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys

if sys.version_info < (3, 9):  # minimum is python3.9
    raise AssertionError(
        "Minimum python version to run Panorama is 3.9. Your current python version is "
        + ".".join(map(str, sys.version_info))
    )

import argparse
from importlib.metadata import distribution

# local modules
from panorama.utils import (
    set_verbosity_level,
    add_common_arguments,
    RawTextArgumentDefaultsHelpFormatter,
)
import panorama.utility
import panorama.info
import panorama.annotate
import panorama.systems
import panorama.alignment
import panorama.compare
import panorama.workflow
import panorama.format.write_flat

version = distribution("panorama").version
opening = """
    ____     ___     _   __   ____     ____     ___     __  ___    ___ 
   / __ \   /   |   / | / /  / __ \   / __ \   /   |   /  |/  /   /   |
  / /_/ /  / /| |  /  |/ /  / / / /  / /_/ /  / /| |  / /|_/ /   / /| |
 / ____/  / ___ | / /|  /  / /_/ /  / _, _/  / ___ | / /  / /   / ___ |
/_/      /_/  |_|/_/ |_/   \____/  /_/ |_|  /_/  |_|/_/  /_/   /_/  |_|

"""
opening_full = (
    opening
    + """

   ______                                       __  _               ____                                               _     
  / ________  ____ ___  ____  ____ __________ _/ /_(__   _____     / __ \____ _____  ____ ____  ____  ____  ____ ___  (______
 / /   / __ \/ __ `__ \/ __ \/ __ `/ ___/ __ `/ __/ | | / / _ \   / /_/ / __ `/ __ \/ __ `/ _ \/ __ \/ __ \/ __ `__ \/ / ___/
/ /___/ /_/ / / / / / / /_/ / /_/ / /  / /_/ / /_/ /| |/ /  __/  / ____/ /_/ / / / / /_/ /  __/ / / / /_/ / / / / / / / /__  
\____/\____/_/ /_/ /_/ .___/\__,_/_/   \__,_/\__/_/ |___/\___/  /_/    \__,_/_/ /_/\__, /\___/_/ /_/\____/_/ /_/ /_/_/\___/  
                    /__                        __              __     __          /____/                                     
  ____ _____  ____ _/ __  __________  _____   / /_____  ____  / _____/ /_  ____  _  __                                       
 / __ `/ __ \/ __ `/ / / / / ___/ _ \/ ___/  / __/ __ \/ __ \/ / ___/ __ \/ __ \| |/_/                                       
/ /_/ / / / / /_/ / / /_/ (__  /  __(__  )  / /_/ /_/ / /_/ / (__  / /_/ / /_/ _>  <                                         
\__,_/_/ /_/\__,_/_/\__, /____/\___/____/   \__/\____/\____/_/____/_.___/\____/_/|_|                                         
                   /____/                                                                                                                                                                   

"""
)
epilog = f"""
By Jérôme Arnoux <jarnoux@genoscope.cns.fr> 
PANORAMA ({version}) is an opensource bioinformatic tools under CeCILL FREE SOFTWARE LICENSE AGREEMENT
LABGeM
"""


def cmd_line():
    """
    Configures and executes the command-line interface for the Panorama software package.

    Returns:
        argparse.Namespace: Parsed command-line arguments and options for execution.

    Raises:
        SystemExit: If invalid arguments are provided, or the help message needs to be displayed.
    """
    desc = "\n"
    desc += (
        "All of the following subcommands have their own set of options. To see them for a given subcommand,"
        " use it with -h or --help, as such:\n"
    )
    desc += "panorama <subcommand> -h\n"
    desc += "\n"
    desc += "     Analyzes:\n"
    desc += "         annotation          Annotate pangenome gene families with HMM or TSV file\n"
    desc += "         systems             Detect systems in pangenome based on one annotation source\n"
    desc += "         write               Writes 'flat' files representing pangenomes that can be used with other software\n"
    desc += "         utility             Some utility command to run analyses more easily\n"
    desc += "         write_systems       Writes 'flat' files about systems detected in pangenomes\n"
    desc += "         pansystems          A workflow to annotate gene families, detect systems and write flat files associated\n"
    desc += "\n"
    desc += "     Compare:\n"
    desc += (
        "         align               Align gene families from multiple pangenomes\n"
    )
    desc += (
        "         cluster             Cluster gene families from multiple pangenomes\n"
    )
    desc += "         compare_context     Compare contexts among pangenomes\n"
    desc += "         compare_systems     Compare contexts among pangenomes\n"
    desc += "         compare_spots       Compare spots among pangenomes\n"
    desc += "\n"
    desc += "     Utils:\n"
    desc += "         info                Provide and compare information through pangenomes\n"
    desc += "         utils               Some utility command to run analyses more easily\n"
    desc += "\n"

    parser = argparse.ArgumentParser(
        description=opening_full,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + version
    )
    subparsers = parser.add_subparsers(
        metavar="", dest="subcommand", title="subcommands", description=desc
    )
    subparsers.required = True  # because python3 sent subcommands to hell apparently

    # print help if no subcommand is specified
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    subs = [
        panorama.info.subparser(subparsers),
        panorama.annotate.subparser(subparsers),
        panorama.systems.detection.subparser(subparsers),
        panorama.alignment.align.subparser(subparsers),
        panorama.alignment.cluster.subparser(subparsers),
        panorama.compare.context.subparser(subparsers),
        panorama.compare.systems.subparser(subparsers),
        panorama.compare.spots.subparser(subparsers),
        panorama.format.write_flat.subparser(subparsers),
        panorama.systems.write_systems.subparser(subparsers),
        panorama.workflow.subparser(subparsers),
        panorama.utility.subparser(subparsers),
    ]

    for sub in subs:  # add options common to all subcommands
        add_common_arguments(sub)
        if sub.prog.split()[1:] == sys.argv[1:]:
            sub.print_help()
            exit(1)
        sub.formatter_class = RawTextArgumentDefaultsHelpFormatter
        sub.description = (
            opening + "\n\n" + sub.description if sub.description else opening
        )
        sub.epilog = epilog

    args = parser.parse_args()
    return args


def main():
    """
    The main function is the entry point for the panorama command line tool.
    It parses arguments and calls subcommands as appropriate.

    :return: The exit status
    """
    args = cmd_line()

    set_verbosity_level(args)

    if args.subcommand == "info":
        panorama.info.launch(args)
    elif args.subcommand == "annotation":
        panorama.annotate.launch(args)
    elif args.subcommand == "systems":
        panorama.systems.detection.launch(args)
    elif args.subcommand == "align":
        panorama.alignment.align.launch(args)
    elif args.subcommand == "cluster":
        panorama.alignment.cluster.launch(args)
    elif args.subcommand == "compare_context":
        panorama.compare.context.launch(args)
    elif args.subcommand == "compare_systems":
        panorama.compare.systems.launch(args)
    elif args.subcommand == "compare_spots":
        panorama.compare.spots.launch(args)
    elif args.subcommand == "write":
        panorama.format.write_flat.launch(args)
    elif args.subcommand == "write_systems":
        panorama.systems.write_systems.launch(args)
    elif args.subcommand == "utils":
        panorama.utility.launch(args)
    elif args.subcommand == "pansystems":
        panorama.workflow.launch(args)


if __name__ == "__main__":
    main()
