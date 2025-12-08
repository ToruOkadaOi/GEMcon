import requests, time, pandas as pd

import geckopy
model = geckopy.io.read_sbml_ec_model('/home/biodata/aman/notebooks/proteome/ecHumanGEM.json', hardcoded_rev_reactions=False)  # how is this made actually? principle? # parameter for faster loading
print('Human Model is imported ... ') # 8+ mins # why ec model taking more time than the normal one?ODD # 68 min is tooflong

df = pd.read_csv('/home/biodata/aman/data/data_raw/paxdb_human_data/9606-Adult_Lung_Kim_2014_SEQUEST.txt', sep='\t', comment="#", header=None, names=["internal_id", "string_external_id", "abundance"], usecols = [0, 1, 2])
# df = df[:19000] # testing

df["ENSP"] = df["string_external_id"].str.replace("9606.", "", regex=False) # this strips down to just the enseml tag

ensps = df["ENSP"].drop_duplicates().tolist()

all_results = []

job = requests.post(
    "https://rest.uniprot.org/idmapping/run",
    data={
        "from": "Ensembl_Protein",
        "to": "UniProtKB",
        "ids": ensps
    }
).json()

job_id = job["jobId"]

time.sleep(3)

# Get results
result = requests.get(
    f"https://rest.uniprot.org/idmapping/stream/{job_id}"
).json()

all_results.extend(result.get("results", []))

map_df = pd.DataFrame(all_results)
map_df = map_df.rename(columns={"from": "ENSP", "to": "UniProt"})

merged = df.merge(map_df, on="ENSP", how="left")

merged["protein_gecko_id"] = merged["UniProt"].apply(lambda x: f"prot_{x}[c]" if pd.notna(x) else None)

print(merged.head(20))
print("Total rows mapped:", merged["UniProt"].notna().sum())

prot = merged[["UniProt", "abundance"]].copy()
prot = prot.dropna(subset=["UniProt"])
prot = prot.rename(columns={"abundance": "copies_per_cell"})

prot["protein_gecko_id"] = prot["UniProt"].apply(lambda x: f"prot_{x}[c]")

prot["stdev"] = 0

# prot
# prot
ec_model_exp = geckopy.experimental.from_copy_number(
    model,
    index=prot["protein_gecko_id"],
    cell_copies = prot["copies_per_cell"],
    stdev = prot["stdev"],
    vol=2.3e-12, dens=1.05, water=0.7
) # pipe output to additional log.f?

#ec_model_exp.slim_optimize()