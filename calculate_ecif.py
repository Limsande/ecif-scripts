#!/usr/bin/env python

"""
Calculate ECIF::LD descriptors for receptor-ligand complexes. Accepts as input
two CSV files with receptors and corresponding poses.

NOTE: Please ignore any AtomValenceExceptions that RDKit might produce.
"""

import argparse
import csv
import os
import sys

import pandas as pd
from ECIF import ecif


def print_error_and_exit(msg: str):
    sys.exit(f'[ERROR] {os.path.basename(sys.argv[0])}: {msg}')


def print_warning(msg: str):
    print(f'[WARNING] {os.path.basename(sys.argv[0])}: {msg}')


def parse_args():
    """
    [--cutoff FLOAT] --receptors FILE --poses FILE --output PATH
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--cutoff', required=False, metavar='FLOAT', type=float,
                        help='Distance cutoff for ECIF calculation. two atoms must not exceed to be counted as pair.'
                             'If --cutoff is omitted, descriptors for all cutoffs from 4.0 to 15.0 A are computed (in'
                             'steps of 0.5 A).In this case, --output is interpreted as a directory, in which one CSV'
                             'file per cutoff is created.')

    required = parser.add_argument_group('required arguments')
    required.add_argument('--receptors', required=True, metavar='FILE', type=str,
                          help='CSV file to read receptor files from. It must contain at least the columns'
                               'ID, RECEPTOR, POSE, and POSERANK. Additional columns are ignored. Column RECEPTOR is'
                               'expected to contain file paths of the molecules. These paths are resolved relative to'
                               'the CSV file. Receptors must be in PDB format. Each file is expected to contain only a'
                               'single molecule. If any of these files does not exist, it is skipped with a warning.'
                               'The columns in the CSV file have to be separated by comma (,).')
    required.add_argument('--poses', required=True, metavar='FILE', type=str,
                          help='CSV file to read docking pose filesfrom. It must contain at least the columns ID, POSE,'
                               'and POSERANK. Additional columns are ignored. Column POSE is expected to contain file'
                               'paths of the molecules. These paths are resolved relative to the CSV file. They must be'
                               'in SD format. Each file is expected to contain only a single molecule. If any of these'
                               'files does not exist, it is skipped with a warning. The columns in the CSV file have to'
                               'be separated by comma (,).')
    required.add_argument('--output', required=True, metavar='PATH', type=str,
                          help='File to write the computed descriptors to. This is a CSV file with columns ID,'
                               'POSERANK, and one column for each element of the descriptors. If PATH contains any'
                               'directories which do not exist, they are created. If PATH already exists and is a file,'
                               'it is overwritten. If --cutoff is not given, PATH is interpreted as a directory instead'
                               'of a file name.')

    return parser.parse_args()


def main(receptor_file: str, pose_file: str,  cutoff: float, output_file: str) -> None:

    # Check if input files exists
    for f in [receptor_file, pose_file]:
        if not os.path.isfile(f):
            print_error_and_exit(f'Input file does not exist: {f}')

    # Create output directory if needed.
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except IOError as e:
            print_error_and_exit(e)

    if not cutoff:
        if not os.path.isdir(args.output):
            try:
                os.makedirs(args.output, exist_ok=True)
            except IOError as e:
                print_error_and_exit(e)

    # Load input CSVs into data frames and merge on ID
    try:
        receptors = pd.read_csv(receptor_file, usecols=['ID', 'RECEPTOR'])
        poses = pd.read_csv(pose_file, usecols=['ID', 'POSE', 'POSERANK'])
    except ValueError as e:
        print_error_and_exit(f'Could not load input file: {e}')

    # sanity check
    if not all(receptors.ID.sort_values().eq(poses.ID.sort_values().unique())):
        print_error_and_exit('IDs in receptor and pose input files do not match!')

    receptor_ligand_pairs = receptors.merge(poses, on='ID')

    if not cutoff:
        cutoffs = [round(x / 10, 1) for x in range(40, 155, 5)]  # [4,15] in 0.5 steps
    else:
        cutoffs = [cutoff]

    for cutoff in cutoffs:
        list_of_descriptor_dfs = []
        for pair in receptor_ligand_pairs.itertuples(index=False):
            # Calculate ECIF::LD descriptors
            try:
                descriptors = ecif.get_ecif_ld(receptor_files=pair.RECEPTOR, ligand_files=pair.POSE, cutoff=cutoff)
            except FileNotFoundError as e:
                print_warning(f'{e}. Skipped.')
                continue

            # attach ID and pose rank to descriptors for nicer output
            descriptors = pd.concat(
                [pd.DataFrame({'ID': [pair.ID], 'POSERANK': [pair.POSERANK]}), descriptors],
                axis='columns')
            list_of_descriptor_dfs.append(descriptors)

        result = pd.concat(list_of_descriptor_dfs, axis='index', ignore_index=True)

        # Treat output as directory, if no cutoff was specified.
        if len(cutoffs) > 1:
            output_file = os.path.join(output_file, f'ECIF_LD_{cutoff}.csv')

        result.to_csv(output_file, index=False)


if __name__ == '__main__':
    args = parse_args()
    main(receptor_file=args.receptors, pose_file=args.poses, cutoff=args.cutoff, output_file=args.output)
