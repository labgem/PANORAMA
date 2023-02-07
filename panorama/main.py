#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys

if sys.version_info < (3, 8):  # minimum is python3.8
    raise AssertionError("Minimum python version to run Panorama is 3.8. Your current python version is " +
                         ".".join(map(str, sys.version_info)))

import argparse
import pkg_resources

# local modules
import panorama.utils
from panorama.utils.overall import check_log, set_verbosity_level
import panorama.info
import panorama.annotate
import panorama.detection
import panorama.dbGraph
import panorama.format.write_flat


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
    desc += "       detection       Detect systems in pangenome based on one annotation source\n"
    desc += "       graph-db        Load pangenomes in Neo4J graph database and allow to perform some queries\n"
    desc += "       write           Writes 'flat' files representing pangenomes that can be used with other software\n"
    desc += "       utility         Some utility command to run analyses more easily\n"
    desc += "\n"

    parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + pkg_resources.get_distribution("panorama").version)
    subparsers = parser.add_subparsers(metavar="", dest="subcommand", title="subcommands", description=desc)
    subparsers.required = True  # because python3 sent subcommands to hell apparently

    subs = [panorama.info.subparser(subparsers),
            panorama.annotate.subparser(subparsers),
            panorama.detection.subparser(subparsers),
            panorama.dbGraph.subparser(subparsers),
            panorama.format.write_flat.subparser(subparsers),
            panorama.utils.subparser(subparsers)]

    for sub in subs:  # add options common to all subcommands
        common = sub._action_groups.pop(1)  # get the 'optional arguments' action group.
        common.title = "Common arguments"
        # common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
        #                     help="directory for storing temporary files")
        common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                            help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
        common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
        common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                            help="disables the progress bars")
        # common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
        common.add_argument('--force', action="store_true",
                            help="Force writing in output directory and in pangenome output file.")
        sub._action_groups.append(common)
        if len(sys.argv) == 2 and sub.prog.split()[1] == sys.argv[1]:
            sub.print_help()
            exit(1)

    args = parser.parse_args()
    return args


def main():
    args = cmd_line()

    set_verbosity_level(args)

    if args.subcommand == "info":
        panorama.info.launch(args)
    elif args.subcommand == "annotation":
        panorama.annotate.launch(args)
    elif args.subcommand == "detection":
        panorama.detection.launch(args)
    elif args.subcommand == "graph-db":
        panorama.dbGraph.launch(args)
    elif args.subcommand == "write":
        panorama.format.write_flat.launch(args)
    elif args.subcommand == "utility":
        panorama.utils.launch(args)


if __name__ == '__main__':
    main()
