import argparse
import sys

from rhea_mapper import workflow
import pkg_resources

VERSION = pkg_resources.get_distribution("rhea_mapper").version
LICENSE = """Copyright (C) Dyliss
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
rhea_mapper is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.\n
"""
MESSAGE = """
Create Rhea sbml files from genomes.
"""
REQUIRES = """
Requires: Orthofinder for genome if you provide a fasta
"""

def main():
    parser = argparse.ArgumentParser(
        "rhea_mapper",
        description=MESSAGE + " For specific help on each subcommand use: rhea_mapper {cmd} --help",
        epilog=REQUIRES
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n" + LICENSE)

    # parent parser
    parent_parser_f = argparse.ArgumentParser(add_help=False)
    parent_parser_f.add_argument(
        "-f",
        "--folder",
        dest="folder",
        help="Input folder containing either annotation, genomes or both",
        required=False,
    )

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Output folder",
        required=False,
    )

    parent_parser_d = argparse.ArgumentParser(add_help=False)
    parent_parser_d.add_argument(
        "-d",
        "--database",
        dest="database",
        help="Database folder, if not existing it will be created",
        required=True,
    )

    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument(
        '-c',
        '--cpu',
        dest="cpu",
        help='cpu number for multi-process',
        required=False,
        type=int,
        default=1)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")

    genome_parser = subparsers.add_parser(
        "genome",
        help="metabolic network reconstruction from genome",
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_d, parent_parser_c
        ],
        description=
        "Run metabolic network reconstruction for each annotated genome of the input directory, using Pathway Tools"
    )

    args = parser.parse_args()

    if not args.cmd:
        print("rhea_mapper " + VERSION + "\n" + LICENSE)
        parser.print_help()
        sys.exit()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if args.cmd == "genome":
        workflow.rhea_mapper_workflow(args.folder, args.output, args.database, args.cpu)
    else:
        print("Invalid command given, rhea_mapper commands are: genome ")
        sys.exit()