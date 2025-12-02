__author__ = "Aman Nalakath"

# Under the assumption that gene names in adata.var are HGNC (alphabetical gene names). If Ensembl or Entrez ##TODO: conversion needed

import celltypist
from celltypist import models
from celltypist import annotate
import scanpy as sc
import os
import argparse

p = argparse.ArgumentParser()
p.add_argument("--input")
args = p.parse_args()

if args.input:
    choice = args.input.strip()
else:
    # auto-detect last normalized file
    folder = "data_processed/processed"
    files = [f for f in os.listdir(folder) if f.endswith("_normalized.h5ad")]

    if not files:
        choice = input("No normalized file found. Provide path manually: ").strip()
    else:
        choice = max(files, key=lambda f: os.path.getmtime(os.path.join(folder, f)))
        choice = os.path.join(folder, choice)
        print(f"Using auto-detected file: {choice}")

basename = os.path.basename(choice)
f_prefix = basename.split('.')[0]

adata = sc.read_h5ad(f"{choice}")

print('check appropriate model at https://www.celltypist.org/models')

print(celltypist.models.models_description()) # cache check?

choice = input("\nChoose which model to proceed with: e.g Immune_All_Low.pkl: ")
if not choice.endswith(".pkl"):
    choice += ".pkl"

models.download_models(model=choice)

# load the model
model = models.Model.load(model=choice)

# change to normed lay ## layer name assumed from previous script
if "normed_log" not in adata.layers:
    raise ValueError("Layer 'normed_log' not found in adata. Check previous script")
adata.X = adata.layers["normed_log"].copy()

# get a new annotated adata
predictions = annotate(adata, model=model, majority_voting=True)

## how i did it
collapsed = predictions.predicted_labels.groupby(level=0)['majority_voting'].first() #groups using cell index(first col), followed by taking majority voting col. and selecting the first element(same for all)
  # check -> collapsed.values

# add predicted labels to adata.obs
adata.obs['cell_type'] = adata.obs_names.map(collapsed)

# save
outdir = "data_processed/processed"
os.makedirs(outdir, exist_ok=True)

output_path = os.path.join(outdir, f"{f_prefix}_annotated.h5ad")
adata.write(output_path)

print("Saved annotated file to:", output_path)
print(output_path)