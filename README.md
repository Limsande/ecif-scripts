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
calculate_ecif.py [-h|--help] [--no-warn-missing-rank] --complx-name STRING --receptor FILE [--cutoff FLOAT] --output PATH FILE...
```

It accepts a receptor and any number of FILE... as ligands. Then computes ECIF::LD descriptors for all receptor-ligand
pairs. ECIF::LD have a parameter, a distance cutoff two atoms must not exceed to be counted as pair. The cutoff to use can be
set via the --cutoff parameter. If --cutoff is omitted, descriptors for all cutoffs from 4.0 to 15.0 A are computed
(in steps of 0.5 A). In this case, OUTPUT is interpreted as a directory, in which one CSV file per cutoff is created.

Example invokation:

```bash
calculate_ecif.py \
    --complx-name 1A0Q \
    --receptor ECIF/Example_Structures/1a0q_protein.pdb \
    --cutoff 6.0 \
    --output ecif_ld.csv \
    ECIF/Example_Structures/1a0q_ligandCD1.pdb
```

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
