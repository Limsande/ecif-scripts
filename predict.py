"""
Predicts binding affinity for given descriptors using the specified trained model.
"""
import os
import pickle
import sys
from argparse import ArgumentParser, Namespace
from typing import Any

import pandas as pd
from pandas import DataFrame


def parse_args() -> Namespace:
    """
    predict.py [-h] --model FILE --descriptors FILE --output FILE
    """
    parser = ArgumentParser(description=__doc__)
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '--model', required=True, type=str, metavar='FILE',
        help='Trained model to make predictions with. Must implement a predict(X) method like sklearn\'s predictors.'
             'Must be saved in binary format as produced by pickle.dump().')
    required.add_argument(
        '--descriptors', required=True, type=str, metavar='FILE',
        help='ECIF::LD descriptors to predict binding affinity for (CSV format). All columns except "Receptor" and'
             '"Ligand" (both optional) are assumed to be part of the descriptors. Descriptors must be of same length as'
             'those used to train the model.')
    required.add_argument(
        '--output', required=True, type=str, metavar='FILE', help='Path to write results to')
    return parser.parse_args()


def load_model(file: str) -> Any:
    with open(file, 'rb') as f:
        model = pickle.load(f)

    # Verify that the model really has a predict method.
    if not callable(getattr(model, 'predict', None)):
        print_error_and_exit('Model has no predict method. Can\'t use it.')

    return model


def load_descriptors(file: str) -> DataFrame:
    return pd.read_csv(file)


def print_error_and_exit(message: str):
    sys.exit(f'[ERROR] {os.path.basename(__file__)}: {message}')


if __name__ == '__main__':
    args = parse_args()

    # Input validation
    for f in [args.model, args.descriptors]:
        if not os.path.isfile(f):
            print_error_and_exit(f'File does not exist or is a directory: {f}')
    if os.path.exists(args.output):
        print_error_and_exit(f'Output file already exists: {args.output}')

    # If output path contains directories, create them.
    if os.path.dirname(args.output) != '':
        try:
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
        except IOError as e:
            print_error_and_exit(e)

    # Load the model
    model = load_model(args.model)

    # Load descriptors
    ecif_ld = load_descriptors(args.descriptors)

    # Predict binding affinities
    predictions = model.predict(ecif_ld.drop(columns=['Receptor', 'Ligand'], errors='ignore'))

    # Construct output
    output = DataFrame()
    if 'Receptor' in ecif_ld.columns:
        output = DataFrame({'Receptor': ecif_ld['Receptor']})
    if 'Ligand' in ecif_ld.columns:
        output = output.join(DataFrame({'Ligand': ecif_ld['Ligand']}))
    output = output.join(DataFrame({'PredictedBindingAffinity': predictions}))
    output.to_csv(args.output, index=False)
