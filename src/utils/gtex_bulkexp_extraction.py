import os
import requests
import pandas as pd
import sys
from pathlib import Path

prj_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(prj_root))

from src.api_hca_userinp import download_file
import asyncio
import argparse

p = argparse.ArgumentParser()
p.add_argument("--tissue", required=True)
args = p.parse_args()
tissue = args.tissue

def walk(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from walk(v)
    elif isinstance(obj, list):
        for item in obj:
            yield from walk(item)
    yield obj

def get_urls(project_id="adult-gtex", keyword="bulk-gex/v10"):
    url = f"https://gtexportal.org/api/v2/dataset/openAccessFilesMetadata?project_id={project_id}"
    data = requests.get(url).json()
    urls = []

    for node in walk(data):
        if isinstance(node, dict) and node.get("url"):
            if keyword in node["url"]:
                urls.append(node["url"])

    return urls

# Fetch only the 54 V10 files
urls = get_urls()

# Build DataFrame
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)
df = pd.DataFrame({
    "name": [os.path.basename(u) for u in urls],
    "url": urls
})

# print(df)

tpm_df = df[df['name'].str.startswith("gene_tpm_v10")].copy()
tpm_df['tissue'] = tpm_df['name'].str.replace('gene_tpm_v10_', '').str.replace('.gct.gz', '')
tissues = tpm_df['tissue'].tolist()
# print(tissues) # STDOUT and seek u.input

# ------------
async def download_tissue(tissue_name, output_path="."):
    row = tpm_df.loc[tpm_df["tissue"] == tissue_name]

    if row.empty:
        print(f"Tissue '{tissue_name}' not found")
        return None

    url = row["url"].iloc[0]
    filename = row["name"].iloc[0]

    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    out_file = output_path / filename

    await download_file(url, str(out_file))

    return str(out_file)


out_file = asyncio.run(download_tissue(tissue)) # argparse

if out_file is None:
    sys.exit(1)

# expression csv extraction ### Need to fix this part

tpm_median = pd.read_csv(out_file, sep="\t", skiprows=2)

# Select Name + first sample # TODO: Multi-samples
df = tpm_median.iloc[:, [0, 2]].copy()

df.columns = ["Name", "expression"]

# Strip Ensembl version numbers
df["Name"] = df["Name"].astype(str).str.split(".").str[0]

bulkcsv_path = f"{prj_root}/data/data_processed/bulk_expr_{tissue}_single_sample.csv"
df.to_csv(bulkcsv_path, index=False)

print(f"Saved file at {bulkcsv_path}")
print(df.head())

# plug to flow.py