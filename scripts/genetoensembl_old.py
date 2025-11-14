import pandas as pd
import re

# load the expression data from before
df = pd.read_csv("expression_data.csv", index_col=0)

# extract mapping from Gencode GTF
mapping = {}
with open("gencode.v43.annotation.gtf") as f:
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
            mapping[gname.group(1)] = gid.group(1)

print(f"Loaded {len(mapping)} gene symbol mappings from Gencode.")

# apply mapping
df["ensembl_id"] = df.index.map(mapping)
mapped = df.dropna(subset=["ensembl_id"]).set_index("ensembl_id")[["expression"]]
print(f"Mapped {len(mapped)}/{len(df)} genes using Gencode reference.")

mapped.to_csv("expression_data_ensembl.csv")
print("Saved: expression_data_ensembl.csv")
