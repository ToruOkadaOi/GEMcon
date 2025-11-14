import scanpy as sc
import numpy as np

choice = input("\nProvide the complete path to the loom file: ")

adata = sc.read_loom(f"{choice}")

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

# Setting a target sum so as not to default to median based normalization
sc.pp.normalize_total(adata, target_sum=1e4) # or 1e6

# log (1+x) transform
sc.pp.log1p(adata)

# creating a new layer unto the object
adata.layers["normed_log"] = adata.X.copy()

# set .X back to raw_counts
adata.X = adata.layers["raw_counts"]

#check
adata.X

# change to .h5ad format for faster loading
adata.write_h5ad(f"{choice}.h5ad")

print(adata.layers.keys())

#cells
## i([adata.obs])

#genes
## i([adata.var])

# Confirm everything
print('Check before proceeding')

print("\nRaw counts (first 20):")
print(adata.layers["raw_counts"].data[:20])

print("\nNormalized + log1p (first 20):")
print(adata.layers["normed_log"].data[:20])

print("\nCells:", adata.obs_names[:5])
print("\nGenes:", adata.var_names[:5])
print("\nShape:", adata.shape)
