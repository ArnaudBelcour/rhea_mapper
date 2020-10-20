
from rhea_mapper import workflow

import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='input folder')
    parser.add_argument('-o', help='output folder')
    parser.add_argument('-d', help='database folder')
    parser.add_argument('-c', '--cpu',
                        help='cpu number for multi-process',
                        required=False,
                        type=int,
                        default=1)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    workflow.rhea_mapper_workflow(args.f, args.o, args.d, args.cpu)