#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import tempfile

import calculate_ecif


def print_error_and_exit(msg: str):
    sys.exit(f'[ERROR] {os.path.basename(sys.argv[0])}: {msg}')


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    # ecif_for_docking_poses.py [-- output-dir DIR] complx_name receptor pose_archive
    parser.add_argument('complx_name', help='Name of receptor-ligand complex, e.g. 1Q11')
    parser.add_argument('receptor', help='Receptor file in PDB format')
    parser.add_argument('--output_dir', help='Directory to create all output in', metavar='DIR')
    parser.add_argument('pose_archive', help='File with docking poses in a format understood by structconvert')

    return parser.parse_args()


def main(complx_name, receptor_file, pose_archive, output_dir):
    for file in [receptor_file, pose_archive]:
        if not os.path.isfile(file):
            print_error_and_exit(f'File not found: {file}')

    if output_dir is not None and not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except IOError as e:
            print_error_and_exit(e)
    elif output_dir is None:
        output_dir = os.getcwd()

    # Run the bash script to split pose archive into several SDF files. The ECIF script wants SDF files.
    with tempfile.TemporaryDirectory() as work_dir:
        try:
            p = subprocess.run(
                ['bash', 'split_pose_archive.sh', pose_archive, work_dir],
                capture_output=True, check=True, text=True)
        except subprocess.CalledProcessError as e:
            print_error_and_exit(f'{e.stdout}\n{e.stderr}')
        else:
            pose_files = p.stdout.strip().split()

        cutoffs = [round(x / 10, 1) for x in range(40, 155, 5)]  # [4,15] in 0.5 steps
        for cutoff in cutoffs:
            output_file = os.path.join(output_dir, f'ECIF_LD_{cutoff}.csv')
            calculate_ecif.main(
                complx_name=complx_name, receptor_file=receptor_file, ligand_files=pose_files,
                cutoff=cutoff, output_file=output_file)


if __name__ == '__main__':
    args = parse_args()
    main(complx_name=args.complx_name, receptor_file=args.receptor, pose_archive=args.pose_archive,
         output_dir=args.output_dir)
