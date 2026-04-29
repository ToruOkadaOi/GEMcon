__author__ = "Aman Nalakath"

import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
from scipy.sparse import issparse
from rich import print
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt

p = argparse.ArgumentParser()
p.add_argument("--input")
p.add_argument("--celltype")
args = p.parse_args()

console = Console()


def col_sums(matrix):
    """
    Sum along axis=0 (per-gene totals across cells), sparse-safe.
    Returns a 1D numpy array of length n_genes.
    """
    if issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).flatten()
    return np.asarray(matrix).sum(axis=0).flatten()


# Receive original input from user; interactive fallback for standalone use
if args.input:
    choice = args.input.strip()
else:
    choice = Prompt.ask(
        "[bold green]\nProvide the complete path to the .h5ad or .loom: "
    )

# Read .loom or .h5ad
if os.path.splitext(os.path.basename(choice))[-1] == ".loom":
    adata = sc.read_loom(f"{choice}")
else:
    adata = sc.read_h5ad(f"{choice}")

base = os.path.splitext(os.path.basename(choice))[0]

print(adata.layers.keys())

# Locate raw counts: check for an existing count layer, else fall back to .X
# IMPORTANT: do NOT copy here — we only need a reference for read-only summing.
raw_source = None
for k in adata.layers:
    if "count" in k.lower():
        raw_source = adata.layers[k]
        break
if raw_source is None:
    raw_source = adata.X

# Sanity check
print(f"raw matrix type: {type(raw_source).__name__}, shape: {raw_source.shape}")
if issparse(raw_source):
    print(f"sparsity: {1 - raw_source.nnz / (raw_source.shape[0] * raw_source.shape[1]):.4f}")

# Show obs and available columns
print(adata.obs)

table = Table(title="obs columns", title_style="bold cyan")
for col_name in adata.obs.columns:
    table.add_row(col_name)
console.print(table)

# Pick a column for pooling (or pool everything)
if args.celltype is not None:
    col = args.celltype.strip()
else:
    col = Prompt.ask(
        "\n[bold green]Enter the column name for pooling-by-celltype -- usually 'cell_type' (press Enter to pool all): "
    ).strip()

genes = adata.var_names

if not col:
    # Pool all cells: per-gene sum across the full dataset
    x = col_sums(raw_source)
    total = x.sum()
    if total == 0:
        raise ValueError("All counts are zero — nothing to normalize.")
    cpm = x / total * 1e6

    os.makedirs("data/data_processed", exist_ok=True)
    output_path = f"data/data_processed/expression_data_{base}.csv"
    pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

else:
    outdir = f"data/data_processed/expression_by_celltype_{base}"
    os.makedirs(outdir, exist_ok=True)

    if col not in adata.obs.columns:
        raise KeyError(
            f"Column '{col}' not found in adata.obs. Available: {list(adata.obs.columns)}"
        )

    cell_types = adata.obs[col].unique()
    console.print(f"[bold]Pooling by '{col}' across {len(cell_types)} cell types[/bold]")

    # Pre-compute boolean masks once; slice the sparse matrix directly per type
    # (avoids constructing AnnData views, which can copy metadata unnecessarily)
    obs_col = adata.obs[col].values

    for t in cell_types:
        mask = obs_col == t
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue

        # Slice rows of the raw matrix — sparse slicing stays sparse and is cheap
        sub_matrix = raw_source[mask, :]
        x = col_sums(sub_matrix)
        total = x.sum()
        if total == 0:
            console.print(f"[yellow]Skipping '{t}' — all zero counts ({n_cells} cells)[/yellow]")
            continue
        cpm = x / total * 1e6

        safe_name = str(t).replace(" ", "_").replace("/", "_").replace("\\", "_")
        output_path = f"{outdir}/{safe_name}.csv"
        pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(output_path, index=False)
        console.rule("[bold green]File Saved[/bold green]")
        console.print(f"[bold cyan]{output_path}[/bold cyan] [dim]({n_cells} cells)[/dim]")