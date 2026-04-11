"""RIPTiDe (Reaction Inclusion by Parsimony and Transcript Distribution) wrapper for GEMcon.

Uses transcriptomic data to contextualize genome-scale metabolic models via
weighted flux minimization (parsimonious FBA weighted by transcript distribution).

Reference:
    Jenior ML, Moutinho Jr TJ, Dougherty BV, & Papin JA. (2020).
    Transcriptome-guided parsimonious flux analysis improves predictions
    with metabolic networks in complex environments. PLOS Comp Biol.
    https://doi.org/10.1371/journal.pcbi.1007099
"""

__author__ = "Aman Nalakath"

import pandas as pd
import numpy as np
import cobra
import os
import re
import argparse
import yaml

import riptide

p = argparse.ArgumentParser(description="RIPTiDe context-specific model extraction")
p.add_argument("--expr", help="Path to expression data csv")
p.add_argument("--model", help="Path to the SBML model e.g. Human-GEM.xml")
args = p.parse_args()

# Load config
cfg = {}
if os.path.exists("config.yaml"):
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f) or {}

model_cfg = cfg.get("model_config", {})

# Transcript pattern (same as gimme.py)
pattern = model_cfg.get("alt_transcript_pattern", "__COBAMPGPRDOT__[0-9]{1}")

# Detect model type
if args.model:
    model_path = args.model.strip()
else:
    model_path = cfg.get("model")

if model_path and (
    "_AT" in pattern or "recon3d" in model_path.lower() or "Recon3D" in model_path
):
    model_type = "recon3d"
    expr_suffix = "_recon3d.csv"
else:
    model_type = "human-gem"
    expr_suffix = "_gencode.csv"

if not model_path:
    model_path = input("\nProvide full path to the metabolic model (.xml): ").strip()

if not os.path.exists(model_path):
    raise FileNotFoundError(f"Model not found: {model_path}")

# Expression file lookup (same logic as gimme.py)
if args.expr:
    expr_path = args.expr.strip()
else:
    files = [
        f
        for f in os.listdir("data/data_processed")
        if f.startswith("expression_data_") and f.endswith(expr_suffix)
    ]

    if not files:
        expr_path = input(
            f"\nNo files found, provide path to expression .csv ({model_type} format): "
        ).strip()
    else:
        paths = [os.path.join("data/data_processed", f) for f in files]
        expr_path = max(paths, key=os.path.getmtime)
        print(f"Using expression file: {expr_path}")

if not os.path.exists(expr_path):
    raise FileNotFoundError(f"Expression file not found: {expr_path}")

base = os.path.splitext(os.path.basename(expr_path))[0]

# ---- Load model ----
print(f"Loading model: {model_path}")
model = cobra.io.read_sbml_model(model_path)

# ---- Load expression data and convert to RIPTiDe transcriptome dict ----
# GEMcon expression CSVs have genes as rows (index) and a value column.
# RIPTiDe expects a dict of {gene_id: [abundance_values]}.
expression_data = pd.read_csv(expr_path, index_col=0)
print(f"Expression data shape: {expression_data.shape}")

# Build transcriptome dict from the DataFrame
# Each gene maps to a list of abundance values (one per column/sample)
transcriptome = {}
for gene_id in expression_data.index:
    values = expression_data.loc[gene_id].values
    if np.isscalar(values):
        transcriptome[str(gene_id)] = [float(values)]
    else:
        transcriptome[str(gene_id)] = [float(v) for v in values]

# RIPTiDe expects a 'replicates' key
n_replicates = len(expression_data.columns)
transcriptome["replicates"] = n_replicates

print(
    f"Transcriptome dict built: {len(transcriptome) - 1} genes, {n_replicates} replicate(s)"
)

# ---- RIPTiDe configuration from config.yaml ----
riptide_cfg = cfg.get("riptide", {})

mode = riptide_cfg.get("mode", "contextualize")  # "contextualize" or "maxfit"
fraction = riptide_cfg.get("fraction", 0.8)
samples = riptide_cfg.get("samples", 500)
prune = riptide_cfg.get("prune", True)
conservative = riptide_cfg.get("conservative", False)
objective = riptide_cfg.get("objective", True)
additive = riptide_cfg.get("additive", False)
set_bounds = riptide_cfg.get("set_bounds", True)
gpr = riptide_cfg.get("gpr", False)
threshold = riptide_cfg.get("threshold", 1e-8)
tasks = riptide_cfg.get("tasks", [])
exclude = riptide_cfg.get("exclude", [])

# maxfit-specific params
frac_min = riptide_cfg.get("frac_min", 0.25)
frac_max = riptide_cfg.get("frac_max", 0.85)
frac_step = riptide_cfg.get("frac_step", 0.1)

# ---- Run RIPTiDe ----
if mode == "maxfit":
    print(
        f"Running RIPTiDe maxfit (fraction range: {frac_min} - {frac_max}, step: {frac_step})"
    )
    riptide_result = riptide.maxfit(
        model=model,
        transcriptome=transcriptome,
        frac_min=frac_min,
        frac_max=frac_max,
        frac_step=frac_step,
        prune=prune,
        samples=samples,
        conservative=conservative,
        objective=objective,
        additive=additive,
        set_bounds=set_bounds,
        tasks=tasks,
        exclude=exclude,
        gpr=gpr,
        threshold=threshold,
    )
else:
    print(f"Running RIPTiDe contextualize (fraction: {fraction})")
    riptide_result = riptide.contextualize(
        model=model,
        transcriptome=transcriptome,
        fraction=fraction,
        samples=samples,
        prune=prune,
        conservative=conservative,
        objective=objective,
        additive=additive,
        set_bounds=set_bounds,
        tasks=tasks,
        exclude=exclude,
        gpr=gpr,
        threshold=threshold,
    )

# ---- Extract context-specific model ----
ctx_model = riptide_result.model

print(
    f"\nThe number of reactions in the base model is: {len(model.reactions)} "
    f"& the number of reactions in the extracted model is: {len(ctx_model.reactions)}"
)

# ---- Save outputs ----
mod_dir = cfg.get("output", {}).get("models_dir")
if not mod_dir:
    mod_dir = (
        input(
            "Enter output directory for saving the models (default './models'): "
        ).strip()
        or "./models"
    )

os.makedirs(mod_dir, exist_ok=True)

# Save the context-specific model as SBML
model_outpath = os.path.join(mod_dir, f"riptide_context_specific_model_{base}.xml")
cobra.io.write_sbml_model(ctx_model, model_outpath)
print(f"Saved context-specific model to: {model_outpath}")

# Save pruned reaction/gene/metabolite info
if hasattr(riptide_result, "pruned") and isinstance(riptide_result.pruned, dict):
    pruned_rxns = list(riptide_result.pruned.get("reactions", []))
    pruned_genes = list(riptide_result.pruned.get("genes", []))
    pruned_mets = list(riptide_result.pruned.get("metabolites", []))

    max_len = max(len(pruned_rxns), len(pruned_genes), len(pruned_mets), 1)
    pruned_df = pd.DataFrame(
        {
            "pruned_reactions": pruned_rxns + [""] * (max_len - len(pruned_rxns)),
            "pruned_genes": pruned_genes + [""] * (max_len - len(pruned_genes)),
            "pruned_metabolites": pruned_mets + [""] * (max_len - len(pruned_mets)),
        }
    )
    pruned_outpath = os.path.join(mod_dir, f"riptide_pruned_{base}.tsv")
    pruned_df.to_csv(pruned_outpath, sep="\t", index=False)
    print(f"Saved pruned components to: {pruned_outpath}")

# Save concordance info
if hasattr(riptide_result, "concordance") and isinstance(
    riptide_result.concordance, dict
):
    r_val = riptide_result.concordance.get("r", "N/A")
    p_val = riptide_result.concordance.get("p", "N/A")
    print(f"Concordance: Spearman r={r_val}, p={p_val}")

# Save flux samples if available
if hasattr(riptide_result, "flux_samples") and isinstance(
    riptide_result.flux_samples, pd.DataFrame
):
    samples_outpath = os.path.join(mod_dir, f"riptide_flux_samples_{base}.tsv")
    riptide_result.flux_samples.to_csv(samples_outpath, sep="\t")
    print(f"Saved flux samples to: {samples_outpath}")

# Print objective value
solution = ctx_model.optimize()
print(f"The objective value is: {solution.objective_value}")
print(f"Fraction of optimum used: {riptide_result.fraction_of_optimum}")
