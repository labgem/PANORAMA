#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys

if sys.version_info < (3, 8):  # minimum is python3.8
    raise AssertionError("Minimum python version to run Panorama is 3.8. Your current python version is " +
                         ".".join(map(str, sys.version_info)))

import argparse
from importlib.metadata import distribution
import logging
from typing import TextIO

# local modules
import panorama.utility
import panorama.info
import panorama.annotate
import panorama.systems
import panorama.alignment
import panorama.compare
import panorama.format.write_flat

version = distribution("panorama").version
epilog = f"""
By Jérôme Arnoux <jarnoux@genoscope.cns.fr> 
PANORAMA ({version}) is an opensource bioinformatic tools under CeCILL FREE SOFTWARE LICENSE AGREEMENT
LABGeM
"""


def check_log(name: str) -> TextIO:
    """Check if the output log is writable

    :param name: Path to the log output

    :return: file object to write log
    """
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        return open(name, "w")


def add_common_arguments(subparser: argparse.ArgumentParser) -> None:
    """
    Add common argument to the input subparser.

    :param subparser: A subparser object from any subcommand.
    """
    common = subparser._action_groups.pop(1)  # get the 'optional arguments' action group.
    common.title = "Common arguments"
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    subparser._action_groups.append(common)


def set_verbosity_level(args: argparse.Namespace) -> None:
    """Set the verbosity level

    :param args: argument pass by command line
    """
    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

        if args.log != sys.stdout and not args.disable_prog_bar:  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True
        str_format = "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s"
        datefmt = '%Y-%m-%d %H:%M:%S'
        if args.log in [sys.stdout, sys.stderr]:
            # use stream
            logging.basicConfig(stream=args.log, level=level,
                                format=str_format,
                                datefmt=datefmt)
        else:
            # log is written in a files. basic condif uses filename
            logging.basicConfig(filename=args.log, level=level,
                                format=str_format,
                                datefmt=datefmt)
        logging.getLogger("PANORAMA").debug("Command: " + " ".join([arg for arg in sys.argv]))
        logging.getLogger("PANORAMA").debug(f"PANORAMA version: {distribution('panorama').version}")


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
    desc += "       write           Writes 'flat' files representing pangenomes some can be used with other software\n"
    desc += "       write_systems   Writes 'flat' files about systems detected in pangenomes\n"
    desc += "       utils           Some utility command to run analyses more easily\n"
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
            panorama.utility.subparser(subparsers)]

    for sub in subs:  # add options common to all subcommands
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


if __name__ == '__main__':
    main()
