import pandas as pd
import mygene
import argparse
import os

p = argparse.ArgumentParser()
p.add_argument("--expr", help="Path to expression CSV with Ensembl IDs")
p.add_argument("--output", help="Output path (optional)")
args = p.parse_args()

if args.expr:
    expr_path = args.expr.strip()
else:
    files = [f for f in os.listdir("data/data_processed") 
             if f.endswith("_gencode.csv")]
    expr_path = max([os.path.join("data/data_processed", f) for f in files], 
                    key=os.path.getmtime) if files else input("Path: ").strip()

if not os.path.exists(expr_path):
    raise FileNotFoundError(expr_path)

df = pd.read_csv(expr_path, index_col=0)
print(f"loaded {len(df)} genes")

mg = mygene.MyGeneInfo()
clean_ids = [str(eid).split('.')[0] for eid in df.index]

results = mg.querymany(clean_ids, scopes='ensembl.gene', 
                       fields='entrezgene', species='human', returnall=True)

mapping = {r['query']: str(int(r['entrezgene'])) 
           for r in results['out'] if 'entrezgene' in r and r['entrezgene']}

print(f"mapped {len(mapping)}/{len(clean_ids)}")

df.index = [mapping.get(cid, cid) for cid in clean_ids]
df = df[df.index.str.isdigit()]

if args.output:
    out_path = args.output
else:
    base = os.path.basename(expr_path).replace('_gencode.csv', '')
    out_path = f"data/data_processed/{base}_recon3d.csv"

os.makedirs(os.path.dirname(out_path), exist_ok=True)
df.to_csv(out_path)
print(f"wrote {len(df)} genes to {out_path}")