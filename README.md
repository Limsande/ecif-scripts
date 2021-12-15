# ECIF-based protein-ligand binding affinity prediction

Based on the paper [ECIF: improving binding affinity prediction through chemical description](https://doi.org/10.1093/bioinformatics/btaa982)
by Sánchez-Cruz et al. and their implementation of ECIF descriptors.

All scripts support the `-h` or `--help` flag to display detailed usage instructions.

## Setup

Make sure all Python packages listed in [`ECIF/requirements.txt`](./ECIF/requirements.txt) are available on your system. The easiest
way is to use conda and import our predefined environment:

```bash
conda env create -n ecif -f ECIF/conda_env.yaml
```

## Calculate ECIF descriptors for protein-ligand pairs

The script `calculate_ecif.py` offers this command line interface:

```bash
calculate_ecif.py [-h|--help] [--cutoff FLOAT] --receptors FILE --poses FILE --output PATH
```

It calculates ECIF::LD descriptors for receptor-ligand complexes. Accepts as input one CSV file with receptor files and one with pose files.
Run `calculate_ecif.py -h` for more information.

Example invocation:

```bash
cd example_data/docking_poses
calculate_ecif.py --cutoff 6 --input poses.csv --output descriptors.csv
```

The Jupyter notebook `notebooks/validate_ecif_calculation.ipynb` compares descriptors computed with our script to descriptors computed with the code by Sánchez-Cruz et al. to show the correctness of our script.

## Train a model on ECIF descriptors and binding affinity data

`train.py` can train models using ECIF descriptors of protein-ligand complexes and the corresponding binding affinity. It currently
supports the same setup as used by Sánchez-Cruz et al.: Either random forest (rf) or gradient boosting trees (gbt) with a predefined
set of hyperparameters.

In normal model (`--output`), the script trains a model and writes it to the specified output file to be used by other scripts.

In evaluation mode (--evaluate), Pearson correlation coefficient and RMSE are computed between true and predicted pK. Instead of relying
on predefined training and validation sets, the evaluation is performed with 10-fold cross validation across the entire training set.

Training data is accepted in three parts, because unlike ECIF, ligand descriptors and pK are always the same for a given PDB ID. Merging
those with all the different possible ECIF beforehand, would mean a lot of redundancy.

```bash
train.py [-h|--help] --model {rf,gbt} --ecif FILE --ld FILE --pK FILE --output FILE
train.py [-h|--help] --model {rf,gbt} --ecif FILE --ld FILE --pK FILE --evaluate
```

## Predict protein-ligand binding affinity with a trained model

`predict.py` predicts binding affinity for given descriptors using the specified trained model.

```bash
predict.py [-h|--help] --model FILE --descriptors FILE --output FILE
```

## Analyze the rescoring

We prepared a Jupyter notebook (`notebooks/rescoring.ipynb`), that can be used to analyze the rescoring.

## TODO

- [ ] resolve paths in input file relative to input file
