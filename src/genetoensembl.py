import pandas as pd
import re
import argparse
import os

p = argparse.ArgumentParser()
p.add_argument("--gtf")
p.add_argument("--expr")
args = p.parse_args()

if args.gtf:
    gtf_path = args.gtf.strip()
else:
    gtf_path = input("\nProvide path to gencode .gtf: ").strip()

if not os.path.exists(gtf_path):
    raise FileNotFoundError(gtf_path)

if args.expr:
    expr_path = args.expr.strip()
else:
    files = [f for f in os.listdir("data/data_processed") if f.startswith("expression_data_") and f.endswith(".csv")]

    if not files:
        expr_path = input("No files found, please provide the abs. path to an expression .csv: ").strip()
    else:
        paths = [os.path.join("data/data_processed", f) for f in files]
        expr_path = max(paths, key=os.path.getmtime)   # the maximum time will be for the last created file ## TODO: verify twice
        print(f"Using the last made file: {expr_path}")

if not os.path.exists(expr_path):
    raise FileNotFoundError(expr_path)

# basename
base = os.path.splitext(os.path.basename(expr_path))[0]

df = pd.read_csv(expr_path, index_col=0)

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
            clean_id = gid.group(1).split('.')[0]
            mapping[gname.group(1)] = clean_id

df["gencode_id"] = df.index.map(mapping)
mapped = df.dropna(subset=["gencode_id"]).set_index("gencode_id")[["expression"]]

# save
out_path = f"data/data_processed/{base}_gencode.csv"
mapped.to_csv(out_path)
print(f"Saved: {out_path}")