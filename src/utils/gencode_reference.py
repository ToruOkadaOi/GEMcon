import pandas as pd
import os
import re
import argparse

p = argparse.ArgumentParser()
p.add_argument("--gtf", help="Path to gencode GTF file")
args = p.parse_args()

if args.gtf:
    gtf_path = args.gtf.strip()
else:
    gtf_path = input("\nProvide the complete path to the gencode gtf file: ").strip()

if not os.path.exists(gtf_path):
    raise FileNotFoundError(f"File not found: {gtf_path}")

regex = re.search(r"gencode\.v\d+", gtf_path)
f_prefix = regex.group(0) if regex else 'gencode'

# I am using gencode.v49 gtf
gtf_data = pd.read_csv(gtf_path, sep="\t", header=6, comment="#") # test header again
gtf_data.columns = ["chr","source","feature","start","end","score","strand","frame","attribute"]

# extract gene_id and gene_name
gtf_data["gene_id"] = gtf_data["attribute"].str.extract('gene_id "([^"]+)"')
gtf_data["gene_name"] = gtf_data["attribute"].str.extract('gene_name "([^"]+)"')

gtf_data = gtf_data[["gene_id", "gene_name"]].drop_duplicates()

print(gtf_data[["gene_id", "gene_name"]].head())

# outfile = f"{f_prefix}_gencode_gene_id.tsv"
# gtf_data.to_csv(outfile, sep="\t", index=False)

outdir = "data_processed/reference"
os.makedirs(outdir, exist_ok=True)

outfile = os.path.join(outdir, f"{f_prefix}_gene_id_map.tsv")
gtf_data.to_csv(outfile, sep="\t", index=False)

print("Saved gene-id map to:", outfile)