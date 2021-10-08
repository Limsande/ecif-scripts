"""
Trains a model of given TYPE, using ECIF, ligand descriptors, and pK from given input FILEs and writes the trained model
to given output FILE.

Instead of relying on predefined training and validation sets, the training is performed with 10-fold cross validation
across the entire training set.

Training data is accepted in three parts, because unlike ECIF, ligand descriptors and pK are always the same for a
given PDB ID. Merging those with all the different possible ECIF beforehand, would mean a lot of redundancy.
"""
import os
import pickle
import sys
from argparse import ArgumentParser, Namespace
from datetime import datetime
from typing import Union

import pandas as pd
from pandas import DataFrame
from scipy.stats import pearsonr
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import make_scorer, mean_squared_error
from sklearn.model_selection import cross_validate


def parse_args() -> Namespace:
    """
    train.py [-h] --model {rf,gbt} --ecif FILE --ld FILE --pK FILE --output FILE
    """
    parser = ArgumentParser(description=__doc__)
    required = parser.add_argument_group('required arguments')
    required.add_argument('--model', required=True, choices=['rf', 'gbt'], type=str, help='Type of model to train')
    required.add_argument(
        '--ecif', required=True, type=str, metavar='FILE',
        help='ECIF portion of training data (CSV format). PDB column is used to join with LD and pK.')
    required.add_argument(
        '--ld', required=True, type=str, metavar='FILE',
        help='Ligand Descriptor portion of training data (CSV format). PDB column is used to join with ECIF and pK.')
    required.add_argument(
        '--pK', required=True, type=str, metavar='FILE',
        help='pK portion of training data (CSV format). PDB column is used to join with ECIF and LD. Columns except'
             'PDB and pK are ignored.')
    required.add_argument(
        '--output', required=True, type=str, metavar='FILE', help='Path to save trained model to')
    return parser.parse_args()


def print_error_and_exit(message: str):
    sys.exit(f'[ERROR] {os.path.basename(__file__)}: {message}')


def load_data(ecif: str, ld: str, pK: str) -> (DataFrame, DataFrame):
    """
    Loads descriptors from given files in CSV format and returns a tuple of (descriptors, pK).
    """
    ecif = pd.read_csv(ecif)
    ligand_descriptors = pd.read_csv(ld)
    pK = pd.read_csv(pK)[['PDB', 'pK']]

    # Join descriptors to make ECIF:LD. Then join the pK values to make sure
    # that they are assigned to the correct PDB ID (since we throw away their IDs
    # later and use a bare list of pK values).
    descriptors = ecif.merge(ligand_descriptors, left_on="PDB", right_on="PDB")
    descriptors_pK = descriptors.merge(pK, left_on='PDB', right_on='PDB')

    return descriptors_pK.iloc[:, :-1], descriptors_pK.pK


def get_model(model: str) -> Union[GradientBoostingRegressor, RandomForestRegressor]:
    """
    Returns model specified by <model>. Parameters are set as described in the ECIF paper.
    """
    if model == 'gbt':
        return GradientBoostingRegressor(
            random_state=42,
            n_estimators=20000,
            max_features="sqrt",
            max_depth=8,
            min_samples_split=3,
            learning_rate=0.005,
            loss="ls",
            subsample=0.7
        )
    elif model == 'rf':
        # TODO
        return RandomForestRegressor()


def pearsonr_score(y_train, y_test) -> float:
    """Wrapper to be used with cross_validate."""
    return pearsonr(y_train, y_test)[0]


def train(model, descriptors, pK) -> dict:
    scoring_funcs = {
        'mse': make_scorer(mean_squared_error, greater_is_better=False),
        'pearsonr': make_scorer(pearsonr_score)
    }
    start_time = datetime.now()
    scores = cross_validate(model, descriptors, pK, scoring=scoring_funcs, cv=10)
    elapsed_time = str(datetime.now() - start_time).split('.')[0]  # Remove microseconds
    scores['elapsed_time'] = elapsed_time
    return scores


if __name__ == '__main__':
    args = parse_args()

    # Input validation
    for f in [args.ecif, args.ld, args.pK]:
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

    # Load the training data
    descriptors, pK = load_data(ecif=args.ecif, ld=args.ld, pK=args.pK)

    model = get_model(args.model)
    print(f'Training model...')
    scores = train(model, descriptors, pK)
    model.scores_dict = scores
    pearson = scores['test_pearsonr'].mean()
    mse = scores['test_mse'].mean() * (-1)  # sign flipped in cross-val because maximization
    print(f'Done. Took {scores["elapsed_time"]}.')
    print('Scores (mean across all CV splits):')
    print(f'  Pearson correlation coefficient: {pearson}')
    print(f'  RMSE: {mse}')
    print('  Scores can be accessed as model.scores_dict')

    print(f'Saving model to {args.output}.')
    pickle.dump(model, open(args.output, 'wb'))
    print('Finished. Bye.')
