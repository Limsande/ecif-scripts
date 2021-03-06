"""
Trains a model of given TYPE, using ECIF, ligand descriptors, and pK from given input FILEs and writes the trained model
to given output FILE.

In evaluation mode (--evaluate), Pearson correlation coefficient and RMSE are computed between true and predicted pK.
Instead of relying on predefined training and validation sets, the evaluation is performed with 10-fold cross validation
across the entire training set.

Training data is accepted in three parts, because unlike ECIF, ligand descriptors and pK are always the same for a
given PDB ID. Merging those with all the different possible ECIF beforehand, would mean a lot of redundancy.
"""
import os
import pickle
import sys
import warnings
from argparse import ArgumentParser, Namespace
from datetime import datetime
from typing import Union

import pandas as pd
from pandas import DataFrame
from scipy.stats import pearsonr, PearsonRConstantInputWarning
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import make_scorer, mean_squared_error
from sklearn.model_selection import cross_validate


def parse_args() -> Namespace:
    """
    train.py [-h] --model {rf,gbt} --ecif FILE --ld FILE --pK FILE --output FILE
    train.py [-h] --model {rf,gbt} --ecif FILE --ld FILE --pK FILE --evaluate
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
    exclusive = parser.add_mutually_exclusive_group()
    exclusive.add_argument(
        '--output', type=str, metavar='FILE', help='Path to save trained model to')
    exclusive.add_argument('--evaluate', action='store_true', help='Only run evaluation')
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

    # Descriptors span all columns, except first (PDB ID) and last (pK).
    return descriptors_pK.iloc[:, 1:-1], descriptors_pK.pK


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
        return RandomForestRegressor(
            random_state=42,
            n_estimators=500,
            max_features=0.33
        )


def pearsonr_score(y_train, y_test) -> float:
    """
    Wrapper to be used with cross_validate. If any of y_train, y_test has 0 variance,
    i.e. is "constant", pearson_r is undefined, and pearsonr() issues a warning. Can't
    think of a better solution than setting pearson_r=0 in this case.
    """
    warnings.filterwarnings('error')
    try:
        res =  pearsonr(y_train, y_test)[0]
    except PearsonRConstantInputWarning:
        res = 0

    return res


def cv_score(model, descriptors, pK) -> (float, float, float):
    scoring_funcs = {
        'mse': make_scorer(mean_squared_error, greater_is_better=False),
        'pearsonr': make_scorer(pearsonr_score)
    }
    start_time = datetime.now()
    scores = cross_validate(model, descriptors, pK, scoring=scoring_funcs, cv=10)
    elapsed_time = str(datetime.now() - start_time).split('.')[0]  # Remove microseconds

    scores['test_mse'] = scores['test_mse'].mean() * (-1)  # sign flipped in cross-val because maximization
    return scores['test_pearsonr'].mean(), scores['test_mse'], elapsed_time


def train(model, descriptors, pK) -> (Union[GradientBoostingRegressor, RandomForestRegressor], float):
    start_time = datetime.now()
    model = model.fit(descriptors, pK)
    elapsed_time = str(datetime.now() - start_time).split('.')[0]  # Remove microseconds
    return model, elapsed_time


if __name__ == '__main__':
    args = parse_args()

    # Input validation
    for f in [args.ecif, args.ld, args.pK]:
        if not os.path.isfile(f):
            print_error_and_exit(f'File does not exist or is a directory: {f}')
    if args.output and os.path.exists(args.output):
        print_error_and_exit(f'Output file already exists: {args.output}')

    # If output path contains directories, create them.
    if args.output and os.path.dirname(args.output) != '':
        try:
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
        except IOError as e:
            print_error_and_exit(e)

    # Load the training data
    descriptors, pK = load_data(ecif=args.ecif, ld=args.ld, pK=args.pK)

    # Train model
    model = get_model(args.model)

    if args.evaluate:
        print(f'Training model for evaluation with 10-fold CV...')
        pearson, mse, elapsed_time = cv_score(model, descriptors, pK)
        print(f'Done. Took {elapsed_time}.')
        print('Scores (mean across all CV splits):')
        print(f'  Pearson correlation coefficient: {pearson}')
        print(f'  RMSE: {mse}')
    else:
        print(f'Training model...')
        model, elapsed_time = train(model, descriptors, pK)
        print(f'Done. Took {elapsed_time}.')

        # Persist model
        model.input_dict = {'ecif': os.path.abspath(args.ecif), 'ld': os.path.abspath(args.ld), 'pK': os.path.abspath(args.pK)}
        print(f'Saving model to {args.output}.')
        with open(args.output, 'wb') as f:
            pickle.dump(model, f)
        print('Input files can be accessed as model.input_dict')
    print('Finished. Bye.')
