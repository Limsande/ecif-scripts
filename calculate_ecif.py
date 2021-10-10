#!/usr/bin/env python

"""
Calculate ECIF::LD descriptors for receptor-ligand complexes.

Accepts a receptor and any number of FILE... as ligands. Then computes ECIF::LD descriptors for all receptor-ligand
pairs. The receptor file must be in PDB and the ligand files in SD format. If a file contains several ligands, the
behavior depends on ECIF/ecif.LoadSDFasDF (which uses rdkit.Chem.MolFromMolFile).

ECIF::LD have a parameter, a distance cutoff two atoms must not exceed to be counted as pair. The cutoff to use can be
set via the --cutoff parameter. If --cutoff is omitted, descriptors for all cutoffs from 4.0 to 15.0 A are computed
(in steps of 0.5 A). In this case, OUTPUT is interpreted as a directory, in which one CSV file per cutoff is created.

Besides the descriptors, output has columns "Receptor" and "Ligand". The first contains the given complex name. The
second is to identify each ligand. If a ligand has the "i_glide_XP_PoseRank" or "i_glide_SP_PoseRank" property, this
value is used. If not, it is assigned the index at which it was given on the command line. By default, this issues a
warning about missing pose rank property. Use --no-warn-missing-rank to change this behavior.

NOTE: Please ignore any AtomValenceExceptions that RDKit might produce.
"""

import argparse
import csv
import os
import sys
from typing import Union

from pandas import DataFrame
from ECIF import ecif
from rdkit.Chem import PandasTools


def print_error_and_exit(msg: str):
    sys.exit(f'[ERROR] {os.path.basename(sys.argv[0])}: {msg}')


def print_warning(msg: str):
    print(f'[WARNING] {os.path.basename(sys.argv[0])}: {msg}')


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    # calculate_ecif.py [--no-warn-missing-rank] --complx-name STRING --receptor FILE [--cutoff FLOAT] --output PATH FILE...
    parser.add_argument('--no-warn-missing-rank', required=False, action='store_true',
                        help='Issue a warning, if any ligand has no pose rank property')
    parser.add_argument('--cutoff', required=False, metavar='FLOAT', type=float,
                        help='Distance cutoff for ECIF calculation')
    required = parser.add_argument_group('required arguments')
    required.add_argument('--complx-name', required=True, metavar='STRING',
                          help='Name of receptor-ligand complex, e.g. 1Q11')
    required.add_argument('--receptor', required=True, metavar='FILE', help='Receptor file in PDB format')
    required.add_argument('--output', required=True, metavar='PATH', help='CSV file write (append) descriptor to')
    required.add_argument('ligands', metavar='FILE', help='Ligand file(s) in SD format',  nargs='+')

    return parser.parse_args()


def main(complx_name, receptor_file, ligand_files, cutoff, output, warn_no_rank=False):

    # Check if all input files exist
    for file in ligand_files + [receptor_file]:
        if file is not None and not os.path.isfile(file):
            print_error_and_exit(f'File not found: {file}')

    # Create output directory if needed.
    if not cutoff:
        if not os.path.isdir(args.output):
            try:
                os.makedirs(args.output, exist_ok=True)
            except IOError as e:
                print_error_and_exit(e)

    if not cutoff:
        cutoffs = [round(x / 10, 1) for x in range(40, 155, 5)]  # [4,15] in 0.5 steps
    else:
        cutoffs = [cutoff]

    for cutoff in cutoffs:
        # Calculate ECIF::LD descriptors
        ecif_ld = ecif.get_ecif_ld(receptor_files=receptor_file, ligand_files=ligand_files, cutoff=cutoff)

        # Construct columns "Receptor" and "Ligand".
        # If pose rank is missing, simply enumerate in same order as files given
        # on command line.
        lig_ids = list(range(1, len(ligand_files) + 1))
        for i, lig_file in enumerate(ligand_files):
            rank = get_pose_rank(lig_file, warn_no_rank=warn_no_rank)
            if rank:
                lig_ids[i] = rank
        ecif_ld = DataFrame({'Receptor': [complx_name] * len(lig_ids), 'Ligand': lig_ids}).join(ecif_ld)

        # Treat output as directory, if no cutoff was specified.
        if len(cutoffs) == 1:
            output_file = output
        else:
            output_file = os.path.join(output, f'ECIF_LD_{cutoff}.csv')

        # Write descriptors to csv file
        if os.path.exists(output_file):
            # Append output to existing file
            write_header = False
        else:
            write_header = True

        ecif_ld.to_csv(output_file, header=write_header, quoting=csv.QUOTE_NONNUMERIC, mode='a', index=False)


def get_pose_rank(sdf_file: str, warn_no_rank) -> Union[str, None]:
    mol = PandasTools.LoadSDF(sdf_file, includeFingerprints=False)
    if 'i_glide_XP_PoseRank' in mol.columns:
        return mol.i_glide_XP_PoseRank.iloc[0]
    elif 'i_glide_SP_PoseRank' in mol.columns:
        return mol.i_glide_SP_PoseRank.iloc[0]
    else:
        if warn_no_rank:
            print_warning(f'{sdf_file}: No pose rank property found')
        return None


if __name__ == '__main__':
    args = parse_args()
    if args.no_warn_missing_rank:
        warn_no_rank = False
    else:
        warn_no_rank = True
    main(complx_name=args.complx_name, receptor_file=args.receptor, ligand_files=args.ligands,
         cutoff=args.cutoff, output=args.output, warn_no_rank=warn_no_rank)
