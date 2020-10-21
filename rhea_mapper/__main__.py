import argparse
import sys

from rhea_mapper import create_database, workflow
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

    # Genome/SPARQL common arguments
    parent_parser_d = argparse.ArgumentParser(add_help=False)
    parent_parser_d.add_argument(
        "-d",
        "--database",
        dest="database",
        help="Database folder, if not existing it will be created",
        required=True,
    )

    # Genome command arguments
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

    # Sparql command arguments
    parent_parser_s = argparse.ArgumentParser(add_help=False)
    parent_parser_s.add_argument(
        "-s",
        "--sparql",
        dest="sparql",
        help="SPARQL query file",
        required=True,
    )

    parent_parser_e = argparse.ArgumentParser(add_help=False)
    parent_parser_e.add_argument(
        "-e",
        "--endpoint",
        dest="endpoint",
        help="Endpoint either 'uniprot' or 'rhea' or 'rhea_endpoint",
        required=True,
    )

    # From_file command argument
    parent_parser_l = argparse.ArgumentParser(add_help=False)
    parent_parser_l.add_argument(
        "-l",
        "--l",
        dest="list",
        help="File containing list of Rhea reacitons IDs or protein IDs",
        required=True,
    )
    parent_parser_r = argparse.ArgumentParser(add_help=False)
    parent_parser_r.add_argument(
        "-r",
        "--reference",
        dest="reference",
        help="Database for reference either 'uniprot' or 'rhea'",
        required=True,
    )

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")

    database_parser = subparsers.add_parser(
        "database",
        help="Create rhea_mapper database [internet connection required]",
        parents=[
            parent_parser_d
        ],
        description=
        "Create rhea_mapper database [internet connection required]"
    )

    genome_parser = subparsers.add_parser(
        "genome",
        help="metabolic network reconstruction from genome",
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_d, parent_parser_c
        ],
        description=
        "metabolic network reconstruction from genome"
    )

    sparql_parser = subparsers.add_parser(
        "sparql",
        help="Using a SPARQL query for reaction ID on Rhea or protein ID on Uniprot, create the corresponding metabolic network",
        parents=[
            parent_parser_s, parent_parser_e, parent_parser_d, parent_parser_o, parent_parser_c
        ],
        description=
        "Using a SPARQL query for reaction ID on Rhea or protein ID on Uniprot, create the corresponding metabolic network"
    )

    list_parser = subparsers.add_parser(
        "from_file",
        help="Using a tsv file with list of Rhea reaction or Uniprot protein, create the corresponding metabolic network",
        parents=[
            parent_parser_l, parent_parser_d, parent_parser_o, parent_parser_c, parent_parser_r
        ],
        description=
        "Using a tsv file with list of Rhea reaction or Uniprot protein, create the corresponding metabolic network"
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
    elif args.cmd == "sparql":
        if args.endpoint in ['uniprot', 'rhea', 'rhea_endpoint']:
            workflow.sparql_query_workflow(args.sparql, args.output, args.database, args.endpoint, args.cpu)
        else:
            print("Invalid endpoint given, either 'uniprot' or 'rhea'")
            sys.exit()
    elif args.cmd == "from_file":
        if args.reference in ['uniprot', 'rhea']:
            workflow.sbml_from_list(args.list, args.database, args.output, args.reference, args.cpu)
        else:
            print("Invalid endpoint given, either 'uniprot' or 'rhea'")
            sys.exit()
    elif args.cmd == "database":
        create_database.download_database(args.database)
    else:
        print("Invalid command given, rhea_mapper commands are: genome, sparql, database")
        sys.exit()
