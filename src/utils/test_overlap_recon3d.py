import pandas as pd
import cobra
import argparse

p = argparse.ArgumentParser()
p.add_argument("--model", help="Path to Recon3D model")
p.add_argument("--expr", help="Path to expression CSV (with Entrez IDs)")
args = p.parse_args()

if args.model:
    model_path = args.model.strip()
else:
    model_path = input("model path: ").strip()

if args.expr:
    expr_path = args.expr.strip()
else:
    expr_path = input("expression path: ").strip()

model = cobra.io.read_sbml_model(model_path)
expr = pd.read_csv(expr_path, index_col=0)

print(f"model: {len(model.genes)} genes")
print(f"expression csv: {len(expr)} genes")

model_genes = set(g.id for g in model.genes)
expr_genes = set(expr.index)

direct = model_genes & expr_genes
model_no_at = set(g.id.split('_AT')[0] for g in model.genes)
expr_no_at = set(str(g).split('_AT')[0] for g in expr.index)
no_at = model_no_at & expr_no_at

print(f"overlap: {len(no_at)} genes")