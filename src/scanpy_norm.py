__author__ = "Aman Nalakath"
import scanpy as sc
import anndata as ad
import os
import argparse
import pandas as pd
from scipy.sparse import csr_matrix, issparse, vstack

p = argparse.ArgumentParser()
p.add_argument("--input")
p.add_argument("--chunk-size", type=int, default=20000,
               help="Number of cells per chunk when streaming loom files")
p.add_argument("--no-compress", action="store_true",
               help="Disable gzip compression on output h5ad (faster write, much bigger file)")
args = p.parse_args()

# Receive the original input from the user, with interactive fallback
if args.input:
    choice = args.input.strip()
else:
    choice = input("\nProvide the complete path to the loom/h5ad file: ")

ext = os.path.splitext(os.path.basename(choice))[-1]


def read_loom_chunked(path, chunk_size=20000):
    """
    Stream a loom file in cell chunks and build a sparse AnnData.
    Avoids ever holding the full dense matrix in memory.
    """
    import loompy
    print(f"Streaming loom file in chunks of {chunk_size} cells...")
    with loompy.connect(path, mode="r") as ds:
        n_genes, n_cells = ds.shape
        print(f"  shape on disk: {n_genes} genes x {n_cells} cells")

        # Cell metadata (column attrs) and gene metadata (row attrs)
        obs = {k: ds.ca[k][:] for k in ds.ca.keys()}
        var = {k: ds.ra[k][:] for k in ds.ra.keys()}

        # Stream the matrix, sparsifying each chunk before appending
        sparse_chunks = []
        for start in range(0, n_cells, chunk_size):
            end = min(start + chunk_size, n_cells)
            # ds[:, start:end] is genes x cells_chunk; transpose so cells are rows
            chunk = ds[:, start:end].T
            sparse_chunks.append(csr_matrix(chunk))
            del chunk
            print(f"  read {end}/{n_cells} cells")

        X = vstack(sparse_chunks).tocsr()
        del sparse_chunks

    obs_df = pd.DataFrame(obs)
    var_df = pd.DataFrame(var)

    # Standard loom convention: CellID for obs, Gene/Accession for var.
    # Prefer Accession (Ensembl IDs) over Gene (symbols) — fewer collisions.
    if "CellID" in obs_df.columns:
        obs_df.index = obs_df["CellID"].astype(str)
    elif "obs_names" in obs_df.columns:
        obs_df.index = obs_df["obs_names"].astype(str)
    else:
        obs_df.index = obs_df.index.astype(str)

    if "Accession" in var_df.columns:
        var_df.index = var_df["Accession"].astype(str)
    elif "Gene" in var_df.columns:
        var_df.index = var_df["Gene"].astype(str)
    elif "var_names" in var_df.columns:
        var_df.index = var_df["var_names"].astype(str)
    else:
        var_df.index = var_df.index.astype(str)

    obs_df.index = ad.utils.make_index_unique(pd.Index(obs_df.index))
    var_df.index = ad.utils.make_index_unique(pd.Index(var_df.index))

    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    return adata


# Read .loom or .h5ad
if ext == ".loom":
    adata = read_loom_chunked(choice, chunk_size=args.chunk_size)
else:
    adata = sc.read_h5ad(choice)

# Ensure adata.X is sparse CSR
if not issparse(adata.X):
    print(f"adata.X loaded as dense ({type(adata.X).__name__}), converting to CSR sparse...")
    adata.X = csr_matrix(adata.X)
elif adata.X.format != "csr":
    adata.X = adata.X.tocsr()

print(f"adata.X type: {type(adata.X).__name__}, shape: {adata.shape}")
print(f"sparsity: {1 - adata.X.nnz / (adata.shape[0] * adata.shape[1]):.4f}")
print("The name/s of the layers in the file downloaded are:", list(adata.layers.keys()))

# Heuristic: detect whether values look already normalized
max_val = adata.X.max()
if max_val < 50:
    print("values look already normalized/log-transformed. SELF-CONFIRM")
    needs_norm = False
else:
    print("Raw counts detected. SELF-CONFIRM")
    needs_norm = True

# Stash raw counts (we need them after in-place normalization).
# We'll delete this layer before writing — .X is restored to raw at the end,
# so saving the same data twice on disk is wasteful.
adata.layers["raw_counts"] = adata.X.copy()

# Normalize in place on .X
if needs_norm:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["normed_log"] = adata.X.copy()
else:
    print("Normalization skipped.")
    adata.layers["normed_log"] = adata.X.copy()

# Restore .X to raw counts
adata.X = adata.layers["raw_counts"]

# Drop the redundant raw_counts layer — .X already holds raw counts
del adata.layers["raw_counts"]

# Save as compressed .h5ad
outdir = "data/data_processed/processed"
os.makedirs(outdir, exist_ok=True)
base = os.path.splitext(os.path.basename(choice))[0]
output_path = os.path.join(outdir, f"{base}_normalized.h5ad")

if args.no_compress:
    adata.write_h5ad(output_path)
else:
    adata.write_h5ad(output_path, compression="gzip", compression_opts=4)

print(f"Saved normalized file to: {output_path}")
print("Layers:", list(adata.layers.keys()))
print(f"On-disk size: {os.path.getsize(output_path) / 1e9:.2f} GB")

# Sanity prints
print("\nCheck before proceeding")
print("\nRaw counts (first 20, from adata.X):")
print(adata.X.data[:20] if issparse(adata.X) else adata.X.flatten()[:20])
norm = adata.layers["normed_log"]
print("\nNormalized + log1p (first 20):")
print(norm.data[:20] if issparse(norm) else norm.flatten()[:20])
print("\nCells:", adata.obs_names[:5].tolist())
print("\nGenes:", adata.var_names[:5].tolist())
print("\nShape:", adata.shape)