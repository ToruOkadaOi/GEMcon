import scanpy as sc
import pandas as pd, numpy as np
import os

choice = input("\nProvide the complete path to the .h5ad or .loom: ") # log1p norm'ed

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
print("obs columns:", adata.obs.columns.tolist())

# user input
col = input("\nEnter the column name for pooling-by-celltype (press Enter to pool all): ").strip()

# gene names usually stored in .var
genes = adata.var_names

# if no column given, pool all cells
if not col:
    # # sum raw counts and convert to CPM; Should I change to array here?
    x = adata.layers["raw_counts"].toarray().sum(axis=0).flatten()
    cpm = x / x.sum() * 1e6

    os.makedirs("data_processed", exist_ok=True)
    output_path = f"data_processed/expression_data_{base}.csv"

    pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

else:
    outdir = f"data_processed/expression_by_celltype_{base}"
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
        print(f"Saved: {output_path}")

#--- normalize ---
#import pandas as pd, numpy as np

# summing all the cells to get a global equivalent like the Bulk RNA-Seq
#gene_sums = np.array(adata.layers["raw_counts"].sum(axis=0)).flatten()

# counts per million for the pooled cells ## To get the number of req. cells to pool - Use library(DSAVE)
#cpm = gene_sums / gene_sums.sum() * 1e6

# export
#pd.DataFrame({"gene": genes, "expression": cpm}).to_csv("expression_data.csv", index=False)

# api_hca_userinp.py -> norm_pre.py -> genetoenseml.py -> gimme/tINIT