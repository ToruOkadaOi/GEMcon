__author__ = "Aman Nalakath"

import os
import re

import pandas as pd
import scanpy as sc
from rich.console import Console

console = Console()


def normalize_scanpy_data(
    input_file: str,
    output_dir: str = "data/data_processed/processed",
    target_sum: float = 1e4,
) -> str:
    if input_file.endswith(".loom"):
        adata = sc.read_loom(input_file)
    else:
        adata = sc.read_h5ad(input_file)

    console.print(f"[blue]Layers: {list(adata.layers.keys())}")

    # Detect if already normalized
    max_val = adata.X.max()
    needs_norm = max_val >= 50
    console.print(
        f"[yellow]{'Raw counts' if needs_norm else 'Already normalized'} detected"
    )

    # Create raw counts layer
    if adata.layers:
        for k in adata.layers:
            if "count" in k.lower():
                adata.layers["raw_counts"] = adata.layers[k].copy()
                break
        else:
            adata.layers["raw_counts"] = adata.X.copy()
    else:
        adata.layers["raw_counts"] = adata.X.copy()

    # Set .X to raw counts
    adata.X = adata.layers["raw_counts"].copy()

    if needs_norm:
        adata.raw = adata.copy()
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        adata.layers["normed_log"] = adata.X.copy()
    else:
        console.print("[yellow]Normalization skipped")
        adata.layers["normed_log"] = adata.X.copy()

    # Set .X back to raw counts
    adata.X = adata.layers["raw_counts"]

    # Save
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(input_file))[0]
    output_path = os.path.join(output_dir, f"{base}_normalized.h5ad")
    adata.write_h5ad(output_path)

    console.print(f"[green]✓ Normalized: {output_path}")
    console.print(f"  Shape: {adata.shape}")

    return output_path


def annotate_celltypes(
    input_file: str,
    model_name: str = "Immune_All_Low.pkl",
    output_dir: str = "data/data_processed/processed",
) -> str:
    # scoped import for this func.
    from celltypist import annotate, models

    adata = sc.read_h5ad(input_file)

    # Ensure model has .pkl extension
    if not model_name.endswith(".pkl"):
        model_name += ".pkl"

    console.print(f"[blue]Loading model: {model_name}")
    models.download_models(model=model_name)
    model = models.Model.load(model=model_name)

    # Switch to normalized layer
    if "normed_log" not in adata.layers:
        raise ValueError(
            "Layer 'normed_log' not found. Run normalize_scanpy_data first"
        )

    adata.X = adata.layers["normed_log"].copy()

    # Annotate
    console.print("[blue]Running annotation...")
    predictions = annotate(adata, model=model, majority_voting=True)

    # Add celltypist predictions to adata.obs
    collapsed = predictions.predicted_labels.groupby(level=0)[
        "majority_voting"
    ].first()
    adata.obs["cell_type"] = adata.obs_names.map(collapsed)

    # Save
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(input_file))[0]
    base = base.replace("_normalized", "")
    output_path = os.path.join(output_dir, f"{base}_annotated.h5ad")
    adata.write(output_path)

    console.print(f"[green]✓ Annotated: {output_path}")
    console.print(
        f"  Cell types found: {len(adata.obs['cell_type'].unique())}"
    )

    return output_path


def normalize_and_pool(
    input_file: str,
    celltype_column: str | None = None,
    output_dir: str = "data/data_processed",
) -> str:
    # Read input
    if input_file.endswith(".loom"):
        adata = sc.read_loom(input_file)
    else:
        adata = sc.read_h5ad(input_file)

    base = os.path.splitext(os.path.basename(input_file))[0]

    # Create raw counts layer
    if adata.layers:
        for k in adata.layers:
            if "count" in k.lower():
                adata.layers["raw_counts"] = adata.layers[k].copy()
                break
        else:
            adata.layers["raw_counts"] = adata.X.copy()
    else:
        adata.layers["raw_counts"] = adata.X.copy()

    adata.X = adata.layers["raw_counts"].copy()
    genes = adata.var_names

    # Pool all cells
    if not celltype_column:
        x = adata.layers["raw_counts"].toarray().sum(axis=0).flatten()
        cpm = x / x.sum() * 1e6

        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"expression_data_{base}.csv")
        pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(
            output_path, index=False
        )

        console.print(f"[green]✓ Pooled all cells: {output_path}")
        return output_path

    # Pool by cell type
    else:
        if celltype_column not in adata.obs.columns:
            raise ValueError(
                f"Column '{celltype_column}' not found in adata.obs"
            )

        outdir = os.path.join(output_dir, f"expression_by_celltype_{base}")
        os.makedirs(outdir, exist_ok=True)

        for cell_type in adata.obs[celltype_column].unique():
            sub = adata[adata.obs[celltype_column] == cell_type]
            x = sub.layers["raw_counts"].toarray().sum(axis=0)
            cpm = x / x.sum() * 1e6

            safe_name = (
                str(cell_type)
                .replace(" ", "_")
                .replace("/", "_")
                .replace("\\", "_")
            )
            output_path = os.path.join(outdir, f"{safe_name}.csv")
            pd.DataFrame({"gene": genes, "expression": cpm}).to_csv(
                output_path, index=False
            )

        console.print(f"[green]✓ Pooled by cell type: {outdir}/")
        console.print(
            f"  Created {len(adata.obs[celltype_column].unique())} files"
        )
        return outdir


def convert_gene_ids(
    gtf_path: str,
    expr_file: str | None = None,
    output_dir: str = "data/data_processed",
) -> str:
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    # Auto detect
    if not expr_file:
        files = [
            f
            for f in os.listdir(output_dir)
            if f.startswith("expression_data_") and f.endswith(".csv")
        ]
        if not files:
            raise FileNotFoundError(
                f"No expression_data_*.csv found in {output_dir}"
            )

        paths = [os.path.join(output_dir, f) for f in files]
        expr_file = max(paths, key=os.path.getmtime)
        console.print(f"[blue]Auto-detected: {expr_file}")

    if not os.path.exists(expr_file):
        raise FileNotFoundError(f"Expression file not found: {expr_file}")

    # Load expression data
    df = pd.read_csv(expr_file, index_col=0)

    # Gene name to ensembl id
    console.print(f"[blue]Parsing GTF: {gtf_path}")
    mapping = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] != "gene":
                continue
            attrs = parts[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            if gid and gname:
                clean_id = gid.group(1).split(".")[0]
                mapping[gname.group(1)] = clean_id

    console.print(f"[blue]Built mapping for {len(mapping)} genes")

    # Convert
    df["gencode_id"] = df.index.map(mapping)
    mapped = df.dropna(subset=["gencode_id"]).set_index("gencode_id")[
        ["expression"]
    ]

    # Save
    base = os.path.splitext(os.path.basename(expr_file))[0]
    output_path = os.path.join(output_dir, f"{base}_gencode.csv")
    mapped.to_csv(output_path)

    console.print(f"[green]✓ Converted: {output_path}")
    console.print(
        f"  Mapped: {len(mapped)} / {len(df)} genes ({len(mapped) / len(df) * 100:.1f}%)"
    )

    return output_path
