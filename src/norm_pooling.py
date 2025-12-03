__author__ = "Aman Nalakath"

import scanpy as sc
import pandas as pd, numpy as np
import os
import argparse
from rich import print
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt

p = argparse.ArgumentParser()
p.add_argument("--input")
p.add_argument("--celltype")
args = p.parse_args()

console = Console()

# This is to make sure it receives the original input at the beginning from the user. option for interactive input in case using this as a single script
if args.input:
    choice = args.input.strip()
else:
    choice = Prompt.ask("[bold green]\nProvide the complete path to the .h5ad or .loom: ")

# determine .loom or .h5ad
if os.path.splitext(os.path.basename(choice))[-1] == '.loom':
    adata = sc.read_loom(f"{choice}")
else:
    adata = sc.read_h5ad(f"{choice}")

# base name from input file
base = os.path.splitext(os.path.basename(choice))[0]

# find the name of layer
print(adata.layers.keys())

# create raw_counts layer
if adata.layers:
    for k in adata.layers:
        if "count" in k.lower():
            adata.layers["raw_counts"] = adata.layers[k].copy()
            break
    else:
        adata.layers["raw_counts"] = adata.X.copy()
else:
    adata.layers["raw_counts"] = adata.X.copy()

# set .X to copy of raw count matrix
adata.X = adata.layers["raw_counts"].copy()

# check the labels
print(adata.obs) # see if there is cell types

# see all the cell types  # to pool by cell type then
table = Table(title="obs columns", title_style="bold cyan") # How to get the table to the center of tty?

for col in adata.obs.columns:
    table.add_row(col)

console.print(table)

# user input
if args.celltype is not None:
    col = args.celltype.strip()
else:
    col = Prompt.ask("\n[bold green]Enter the column name for pooling-by-celltype -- usually 'cell_type' (press Enter to pool all): ").strip()
# gene names usually stored in .var
genes = adata.var_names

# if no column given, pool all cells
if not col:
    # # sum raw counts and convert to CPM; Should I change to array here?
    x = adata.layers["raw_counts"].toarray().sum(axis=0).flatten()
    cpm = x / x.sum() * 1e6

    os.makedirs("data/data_processed", exist_ok=True)
    output_path = f"data/data_processed/expression_data_{base}.csv"

    pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

else:
    outdir = f"data/data_processed/expression_by_celltype_{base}"
    os.makedirs(outdir, exist_ok=True)
    for t in adata.obs[col].unique():
        # subset for this cell type
        sub = adata[adata.obs[col] == t]
        
        # sum raw counts and convert to CPM
        x = sub.layers["raw_counts"].toarray().sum(axis=0)
        cpm = x / x.sum() * 1e6
        
        # sanitize name
        safe_name = str(t).replace(' ', '_').replace('/', '_').replace('\\', '_')

        # path for save
        output_path = f"{outdir}/{safe_name}.csv"
        
        # save files
        pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(output_path, index=False)
        console.rule("[bold green]File Saved[/bold green]")
        console.print(f"[bold cyan]{output_path}[/bold cyan]")