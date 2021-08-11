# Rescoring Glide Fragment Docking Poses for LEADS-FRAG dataset using Binding Affinity Prediction with ECIF

## Setup

Make sure all Python packages listed in [`ECIF/requirements.txt`](./ECIF/requirements.txt) are available on your system. The easiest way is to use conda and import our predefined environment:

```bash
conda env create -n ecif -f ECIF/conda_env.yaml
```

## Calculate ECIF descriptors for LEADS-FRAG dataset

The script offers this command line interface:

```bash
calculate_ecif.py <receptor> <ligand> <cutoff> <output_file>
```

Example invokation:

```bash
calculate_ecif.py \
    ECIF/Example_Structures/1a0q_protein.pdb \
    ECIF/Example_Structures/1a0q_ligandCD1.pdb \
    6.0 \
    ecif_ld.csv
```
