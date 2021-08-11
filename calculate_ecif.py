"""
Calculate ECIF descriptor with given cutoff for given receptor-ligand complex.

NOTE: Please ignore any AtomValenceExceptions that RDKit might produce.
"""

import argparse
import csv
import os
import sys

from ECIF import ecif


def print_error_and_exit(msg: str):
    sys.exit(f'[ERROR] {sys.argv[0]}: {msg}')


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument('complx_name', help='Name of receptor-ligand complex, e.g. 1Q11')
    parser.add_argument('receptor', help='Receptor file in PDB format')
    parser.add_argument('ligand', help='Ligand file in SDF format')
    parser.add_argument('cutoff', type=float, help='Distance cutoff for ECIF calculation')
    parser.add_argument('output_file', help='CSV file write (append) descriptor to')

    return parser.parse_args()


def main(complx_name, receptor_file, ligand_file, cutoff, output_file):

    for file in [receptor_file, ligand_file]:
        if not os.path.isfile(file):
            print_error_and_exit(f'File not found: {file}')

    ecif_ld = ecif.get_ecif_ld_for_single_complex(
        complx_name=complx_name, receptor_file=receptor_file,ligand_file=ligand_file, cutoff=cutoff)

    if os.path.exists(output_file):
        # Append output to existing file
        write_header = False
    else:
        write_header = True
    ecif_ld.to_csv(output_file, header=write_header, quoting=csv.QUOTE_NONNUMERIC, mode='a', index=False)


if __name__ == '__main__':
    args = parse_args()
    main(complx_name=args.complx_name, receptor_file=args.receptor, ligand_file=args.ligand,
         cutoff=args.cutoff, output_file=args.output_file)
