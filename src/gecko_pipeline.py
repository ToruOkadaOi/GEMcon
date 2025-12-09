__author__ = "Aman Nalakath"

import geckopy
from geckopy.io import read_sbml_ec_model
import argparse, os, requests, time, pandas as pd, yaml

p = argparse.ArgumentParser()
p.add_argument("--expr", help="Path to the protein abundance csv file")
p.add_argument("--model", help="Path to the sbml model eg Human-GEM.xml")
args = p.parse_args()

cfg = {}
if os.path.exists("config.yaml"):
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f) or {}
# geckopy part
geckopy_cfg = cfg.get('geckopy', {})

cell_params = geckopy_cfg.get("cell_params", {})
vol = cell_params.get("vol", 2.3e-12)
dens = cell_params.get("dens", 1.05)
water = cell_params.get("water", 0.7)

## Model
# if model (eg. Human-GEM or Recon3D) is specified as cli arg
if args.model:
    model_path = args.model.strip()
else:
    model_path = geckopy_cfg.get("ec_model")

if not model_path:
    model_path = input("\nProvide full path to the expanded/normal metabolic model (.xml): ").strip()

if not os.path.exists(model_path):
    raise FileNotFoundError(model_path)

## Abundance file
if args.expr:
    expr_path = args.expr.strip()

else:
    expr_path = geckopy_cfg.get("paxdb_data")

if not expr_path:
    expr_path = input("\nProvide path to PaxDB protein abundance file: ").strip()

if not os.path.exists(expr_path):
    raise FileNotFoundError(expr_path)

# Using the expanded model here to save on the loading time; from geckopy.io.write_sbml_ec_model()
model = read_sbml_ec_model(model_path)

df = pd.read_csv(expr_path, sep='\t', comment="#", header=None, names=["internal_id", "string_external_id", "abundance"], usecols = [0, 1, 2])

df["ENSP"] = df["string_external_id"].str.replace("9606.", "", regex=False) # this strips down to just the enseml tag

ensps = df["ENSP"].drop_duplicates().tolist()

all_results = []

try:
    job = requests.post(
        "https://rest.uniprot.org/idmapping/run",
        data={
            "from": "Ensembl_Protein",
            "to": "UniProtKB",
            "ids": ensps
        }
    ).json()

    job_id = job["jobId"]

except Exception as e:
    raise RuntimeError(f"Failed to post: {e}") from e

time.sleep(5)

try:
    # Get results
    result = requests.get(
        f"https://rest.uniprot.org/idmapping/stream/{job_id}"
    ).json()

except Exception as e:
    raise RuntimeError(f"Failed to get: {e}") from e

all_results.extend(result.get("results", []))

map_df = pd.DataFrame(all_results)

map_df = map_df.rename(columns={"from": "ENSP", "to": "UniProt"})

merged = df.merge(map_df, on="ENSP", how="left") # keeps all rows

merged["protein_gecko_id"] = merged["UniProt"].apply(lambda x: f"prot_{x}[c]" if pd.notna(x) else None) # to match the print(model.proteins)

#print(merged.head(20))
print("Total rows mapped: ", merged["UniProt"].notna().sum())

# switching to a format geckopy expects
prot = merged[["UniProt", "abundance"]].copy()

prot = prot.dropna(subset=["UniProt"])

prot = prot.rename(columns={"abundance": "copies_per_cell"})

prot["protein_gecko_id"] = prot["UniProt"].apply(lambda x: f"prot_{x}[c]")

# standard deviation set to zero; Sample infos not available
prot["stdev"] = 0

print(f"The proteins in the model are like this: {model.proteins}\n")

# prot # 39mins!!!
ec_model_exp = geckopy.experimental.from_copy_number(
    model,
    index=prot["protein_gecko_id"],
    cell_copies = prot["copies_per_cell"],
    stdev = prot["stdev"],
    vol=vol, dens=dens, water=water # are these close to bio.reality? # test values only
) # pipe output to additional log.f?

print("The extracted model is: ", ec_model_exp.summary())

# optimize
print(f"The objective is: {ec_model_exp.slim_optimize()}")

# SAVE

base = os.path.splitext(os.path.basename(expr_path))[0]

# get the save path from config
mod_dir = cfg.get('output', {}).get('models_dir')
os.makedirs(mod_dir, exist_ok=True)
if not mod_dir:
    mod_dir = input("Enter output directory for saving the models (default './models'): ").strip() or './models'

# export the model
geckopy.io.write_sbml_ec_model(ec_model_exp, f"{mod_dir}/geckopy_context_specific_model_{base}.xml")