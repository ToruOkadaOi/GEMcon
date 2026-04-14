"""CORDA wrapper for GEMcon.

CORDA (Cost Optimization Reaction Dependency Assessment) reconstructs
context-specific metabolic networks using confidence levels for reactions.
Unlike other algorithms, CORDA works with free solvers (GLPK) and uses a
5-level confidence scoring system.

Reference:
    Schultz A, et al. Reconstruction of Tissue-Specific Metabolic Networks
    Using CORDA. PLoS Comput Biol. 2016;12(3):e1004808.
"""

__author__ = "Aman Nalakath"

import pandas as pd
import cobra
import os
import argparse
import yaml

from corda import CORDA, reaction_confidence

p = argparse.ArgumentParser(description="CORDA context-specific model extraction")
p.add_argument("--expr", help="Path to expression data csv")
p.add_argument("--model", help="Path to the SBML model e.g. Human-GEM.xml")
args = p.parse_args()

# Load config
cfg = {}
if os.path.exists("config.yaml"):
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f) or {}

model_cfg = cfg.get("model_config", {})

# Detect model type
if args.model:
    model_path = args.model.strip()
else:
    model_path = cfg.get("model")

pattern = model_cfg.get("alt_transcript_pattern", "__COBAMPGPRDOT__[0-9]{1}")

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

# Expression file lookup
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

# Load model
print(f"Loading model: {model_path}")
model = cobra.io.read_sbml_model(model_path)
model.solver = (
    "cplex"  # choosing the faster solver (change to "glpk" if CPLEX unavailable)
)

# Load expression data
expression_data = pd.read_csv(expr_path, index_col=0)
print(f"Expression data shape: {expression_data.shape}")

# Get CORDA config
corda_cfg = cfg.get("corda", {})

# Confidence thresholds (for log-transformed scRNA-seq data)
# These map expression values to confidence levels: -1, 0, 1, 2, 3
# Default quantiles: 5th, 25th, 50th, 75th percentile
quantile_absent = corda_cfg.get("quantile_absent", 0.05)
quantile_unknown = corda_cfg.get("quantile_unknown", 0.25)
quantile_low = corda_cfg.get("quantile_low", 0.50)
quantile_medium = corda_cfg.get("quantile_medium", 0.75)

# Or use absolute thresholds if provided
threshold_absent = corda_cfg.get("threshold_absent", None)
threshold_unknown = corda_cfg.get("threshold_unknown", None)
threshold_low = corda_cfg.get("threshold_low", None)
threshold_medium = corda_cfg.get("threshold_medium", None)

met_prod = corda_cfg.get("met_prod", [])  # Metabolites that must be producible

print("Converting expression to gene confidence levels...")

# Get expression values (single sample or mean across samples)
if expression_data.shape[1] == 1:
    gene_expr = expression_data.iloc[:, 0]
else:
    gene_expr = expression_data.mean(axis=1)

# Calculate thresholds from quantiles if not provided absolutely
if threshold_absent is None:
    threshold_absent = gene_expr.quantile(quantile_absent)
if threshold_unknown is None:
    threshold_unknown = gene_expr.quantile(quantile_unknown)
if threshold_low is None:
    threshold_low = gene_expr.quantile(quantile_low)
if threshold_medium is None:
    threshold_medium = gene_expr.quantile(quantile_medium)

print("Confidence thresholds:")
print(f"  Absent (-1): < {threshold_absent:.4f}")
print(f"  Unknown (0): {threshold_absent:.4f} - {threshold_unknown:.4f}")
print(f"  Low (1): {threshold_unknown:.4f} - {threshold_low:.4f}")
print(f"  Medium (2): {threshold_low:.4f} - {threshold_medium:.4f}")
print(f"  High (3): > {threshold_medium:.4f}")


# Map gene expression to confidence levels
def expr_to_confidence(val):
    if val <= threshold_absent:
        return -1  # absent
    elif val <= threshold_unknown:
        return 0  # unknown
    elif val <= threshold_low:
        return 1  # low confidence
    elif val <= threshold_medium:
        return 2  # medium confidence
    else:
        return 3  # high confidence


# Build gene confidence dictionary
conf_genes = {}
for gene_id, expr_val in gene_expr.items():
    # Strip transcript version if present (e.g., ENSG000001234.5 -> ENSG000001234)
    gene_clean = str(gene_id).split(".")[0]
    conf_genes[gene_clean] = expr_to_confidence(expr_val)

print("Gene confidence distribution:")
print(f"  Absent (-1): {sum(1 for v in conf_genes.values() if v == -1)}")
print(f"  Unknown (0): {sum(1 for v in conf_genes.values() if v == 0)}")
print(f"  Low (1): {sum(1 for v in conf_genes.values() if v == 1)}")
print(f"  Medium (2): {sum(1 for v in conf_genes.values() if v == 2)}")
print(f"  High (3): {sum(1 for v in conf_genes.values() if v == 3)}")

# Map gene confidences to reaction confidences using GPR rules
print("Mapping gene confidences to reaction confidences using GPR rules...")
conf_reactions = {}
for rxn in model.reactions:
    if rxn.gene_reaction_rule:
        conf = reaction_confidence(rxn, conf_genes)
        conf_reactions[rxn.id] = conf
    else:
        # No GPR rule - unknown confidence
        conf_reactions[rxn.id] = 0

print("Reaction confidence distribution:")
print(f"  Absent (-1): {sum(1 for v in conf_reactions.values() if v == -1)}")
print(f"  Unknown (0): {sum(1 for v in conf_reactions.values() if v == 0)}")
print(f"  Low (1): {sum(1 for v in conf_reactions.values() if v == 1)}")
print(f"  Medium (2): {sum(1 for v in conf_reactions.values() if v == 2)}")
print(f"  High (3): {sum(1 for v in conf_reactions.values() if v == 3)}")

# Run CORDA
print("Running CORDA reconstruction...")
corda_obj = CORDA(model, conf_reactions, met_prod=met_prod)
corda_obj.build()

# Get the reconstructed model
print("CORDA reconstruction complete!")
print(f"Included reactions: {corda_obj.included}")

# Build context-specific model
ctx_rxns = [r.id for r in model.reactions if corda_obj.included.get(r.id, False)]
ctx_model = model.copy()
to_remove = [r for r in ctx_model.reactions if r.id not in ctx_rxns]
ctx_model.remove_reactions(to_remove, remove_orphans=True)

print(
    f"\nThe number of reactions in the base model is: {len(model.reactions)} "
    f"& the number of reactions in the extracted model is: {len(ctx_model.reactions)}"
)

# Save output
mod_dir = cfg.get("output", {}).get("models_dir", "./models")
os.makedirs(mod_dir, exist_ok=True)

# Save the context-specific model
output_path = os.path.join(mod_dir, f"corda_context_specific_model_{base}.xml")
cobra.io.write_sbml_model(ctx_model, output_path)
print(f"Saved context-specific model to: {output_path}")

# Save reaction confidences
conf_df = pd.DataFrame(
    [{"reaction_id": k, "confidence": v} for k, v in conf_reactions.items()]
)
conf_path = os.path.join(mod_dir, f"corda_reaction_confidences_{base}.tsv")
# Use csv module to avoid pandas version compatibility issues
import csv

with open(conf_path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["reaction_id", "confidence"])
    for k, v in conf_reactions.items():
        writer.writerow([k, v])
print(f"Saved reaction confidences to: {conf_path}")

# Print objective value
solution = ctx_model.optimize()
print(f"The objective value is: {solution.objective_value}")
