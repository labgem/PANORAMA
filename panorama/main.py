#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys

if sys.version_info < (3, 8):  # minimum is python3.6
    raise AssertionError("Minimum python version to run Panorama is 3.8. Your current python version is " +
                         ".".join(map(str, sys.version_info)))

import argparse
import logging
import pkg_resources
# import tempfile
# import os

import panorama.PanFluidity


def check_log(name):
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        return open(name, "w")


def cmd_line():
    # need to manually write the description so that it's displayed into groups of subcommands ....
    desc = "\n"
    desc += "All of the following subcommands have their own set of options. To see them for a given subcommand," \
            " use it with -h or --help, as such:\n"
    desc += "  panorama <subcommand> -h\n"
    desc += "\n"
    desc += "  Global:\n"
    desc += "    panfluidity           Compute genome fluidity of all pangenomes and compare them\n"
    desc += "\n"

    parser = argparse.ArgumentParser(
        description="Comparative pangenomic  toolsbox",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + pkg_resources.get_distribution("panorama").version)
    subparsers = parser.add_subparsers(metavar="", dest="subcommand", title="subcommands", description=desc)
    subparsers.required = True  # because python3 sent subcommands to hell apparently

    subs = [panorama.PanFluidity.subparser(subparsers)]

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
        common.add_argument('-f', '--force', action="store_true",
                            help="Force writing in output directory and in pangenome output file.")
        sub._action_groups.append(common)
        if len(sys.argv) == 2 and sub.prog.split()[1] == sys.argv[1]:
            sub.print_help()
            exit(1)

        args = parser.parse_args()

        return args


def main():
    args = cmd_line()

    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

        if args.log != sys.stdout and not args.disable_prog_bar:  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True

        logging.basicConfig(stream=args.log, level=level,
                            format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.getLogger().info("Command: " + " ".join([arg for arg in sys.argv]))
        logging.getLogger().info("Panorama version: " + pkg_resources.get_distribution("panorama").version)

    if args.subcommand == "panfluidity":
        panorama.PanFluidity.launch(args)


if __name__ == '__main__':
    main()
