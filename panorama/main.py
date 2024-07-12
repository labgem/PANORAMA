#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys

if sys.version_info < (3, 8):  # minimum is python3.8
    raise AssertionError("Minimum python version to run Panorama is 3.8. Your current python version is " +
                         ".".join(map(str, sys.version_info)))

import argparse
from importlib.metadata import distribution

# local modules
from panorama.utils import set_verbosity_level, add_common_arguments
import panorama.utility
import panorama.info
import panorama.annotate
import panorama.systems
import panorama.alignment
import panorama.compare
import panorama.workflow
import panorama.format.write_flat


version = distribution("panorama").version
epilog = f"""
By Jérôme Arnoux <jarnoux@genoscope.cns.fr> 
PANORAMA ({version}) is an opensource bioinformatic tools under CeCILL FREE SOFTWARE LICENSE AGREEMENT
LABGeM
"""


def cmd_line():
    # need to manually write the description so that it's displayed into groups of subcommands ....
    desc = "\n"
    desc += "All of the following subcommands have their own set of options. To see them for a given subcommand," \
            " use it with -h or --help, as such:\n"
    desc += "  panorama <subcommand> -h\n"
    desc += "\n"
    desc += "  Global:\n"
    desc += "       info            Provide and compare information through pangenomes\n"
    desc += "       annotation      Annotate pangenome gene families with HMM or TSV file\n"
    desc += "       systems       Detect systems in pangenome based on one annotation source\n"
    desc += "       compare         Pangenome comparison methods\n"
    desc += "       align           Align gene families from multiple pangenomes\n"
    desc += "       cluster         Cluster gene families from multiple pangenomes\n"
    desc += "       compare         Compare contexts and modules across multiple pangenomes\n"    
    desc += "       write           Writes 'flat' files representing pangenomes that can be used with other software\n"
    desc += "       utility         Some utility command to run analyses more easily\n"
    desc += "       write_systems   Writes 'flat' files about systems detected in pangenomes\n"
    desc += "       utils           Some utility command to run analyses more easily\n"
    desc += "       pansystems      A workflow to annotate gene families, detect systems and write flat files associated\n"
    desc += "\n"

    parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + version)
    subparsers = parser.add_subparsers(metavar="", dest="subcommand", title="subcommands", description=desc)
    subparsers.required = True  # because python3 sent subcommands to hell apparently

    # print help if no subcommand is specified
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    subs = [panorama.info.subparser(subparsers),
            panorama.annotate.subparser(subparsers),
            panorama.systems.detection.subparser(subparsers),
            panorama.alignment.align.subparser(subparsers),
            panorama.alignment.cluster.subparser(subparsers),
            panorama.compare.subparser(subparsers),
            panorama.format.write_flat.subparser(subparsers),
            panorama.systems.write_systems.subparser(subparsers),
            panorama.workflow.subparser(subparsers),
            panorama.utility.subparser(subparsers)]

    for sub in subs:  # add options common to all subcommands
        # common = sub._action_groups.pop(1)  # get the 'optional arguments' action group.
        add_common_arguments(sub)
        # launch help when no argument is given except the command
        # sub.prog content examples that trigger print_help:
        # panorama compare
        # panorama info
        # panorama compare context
        if sub.prog.split()[1:] == sys.argv[1:]:
            sub.print_help()
            exit(1)

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
    elif args.subcommand == "compare":
        panorama.compare.launch(args)
    elif args.subcommand == "write":
        panorama.format.write_flat.launch(args)
    elif args.subcommand == "write_systems":
        panorama.systems.write_systems.launch(args)
    elif args.subcommand == "utils":
        panorama.utility.launch(args)
    elif args.subcommand == "pansystems":
        panorama.workflow.launch(args)


if __name__ == '__main__':
    main()
