__author__ = "Aman Nalakath"

import scanpy as sc
import numpy as np

choice = input("\nProvide the complete path to the loom file: ")

adata = sc.read_loom(f"{choice}")

# find the name of layers
print('The name/s of the layers in the file downloaded are: ',adata.layers.keys())

# to detect whether values look already normlized
max_val = adata.X.max()

if max_val < 50:  # heuristic threshold
    print("values look already normalized/log-transformed. SELF-CONFIRM")
    needs_norm = False
else:
    print("Raw counts detected. SELF-CONFIRM")
    needs_norm = True

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

if needs_norm:
    # save true raw counts permanently ## Is this nec.?
    adata.raw = adata.copy()

    # Setting a target sum so as not to default to median based normalization
    sc.pp.normalize_total(adata, target_sum=1e4) # or 1e6

    # log (1+x) transform
    sc.pp.log1p(adata)

    # creating a new layer unto the object
    adata.layers["normed_log"] = adata.X.copy()

else:
    print("Normalization skipped.")

    adata.layers["normed_log"] = adata.X.copy()

# set .X back to raw_counts
adata.X = adata.layers["raw_counts"]

#check; #print(adata.X)

# change to .h5ad format for faster loading
output_path = choice.replace('.loom', '.h5ad')
adata.write_h5ad(output_path)

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