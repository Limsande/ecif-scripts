#!/usr/bin/env python

"""
Calculate ECIF descriptor with given cutoff for given receptor-ligand complex.

NOTE: Please ignore any AtomValenceExceptions that RDKit might produce.
"""

import argparse
import csv
import os
import sys

from pandas import DataFrame
from ECIF import ecif


def print_error_and_exit(msg: str):
    sys.exit(f'[ERROR] {os.path.basename(sys.argv[0])}: {msg}')


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    # calculate_ecif.py complx_name receptor cutoff output_file ligand...
    parser.add_argument('complx_name', help='Name of receptor-ligand complex, e.g. 1Q11')
    parser.add_argument('receptor', help='Receptor file in PDB format')
    parser.add_argument('cutoff', type=float, help='Distance cutoff for ECIF calculation')
    parser.add_argument('output_file', help='CSV file write (append) descriptor to')
    parser.add_argument('ligand', help='Ligand file in SD format',  nargs='+')

    return parser.parse_args()


def main(complx_name, receptor_file, ligand_files, cutoff, output_file):

    # Check if all input files exist
    for file in ligand_files + [receptor_file]:
        if file is not None and not os.path.isfile(file):
            print_error_and_exit(f'File not found: {file}')

    # Calculate ECIF::LD descriptors
    ecif_ld = ecif.get_ecif_ld(receptor_files=receptor_file, ligand_files=ligand_files, cutoff=cutoff)

    # Add complex name, and pose number if number of ligands >1
    if len(ligand_files) > 1:
        ecif_ld = DataFrame(
            {'Name': [complx_name] * len(ecif_ld),
             'Ligand_idx': list(range(1, len(ecif_ld) + 1))}
        ).join(ecif_ld)
    else:
        ecif_ld = DataFrame(
            {'Name': [complx_name] * len(ecif_ld)}
        ).join(ecif_ld)

    # Write descriptors to csv file
    if os.path.exists(output_file):
        # Append output to existing file
        write_header = False
    else:
        write_header = True
    ecif_ld.to_csv(output_file, header=write_header, quoting=csv.QUOTE_NONNUMERIC, mode='a', index=False)


if __name__ == '__main__':
    args = parse_args()
    main(complx_name=args.complx_name, receptor_file=args.receptor, ligand_files=args.ligand,
         cutoff=args.cutoff, output_file=args.output_file)
