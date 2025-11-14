#--- code written with .loom files in mind i.e DCP processed matrix with just 1 file in the project. Many projects have multiple count matrix files. # TODO: think of scrapping all this together if user is going to use those instead ## add before pd.DataFrame printing ---#
#TODO: make changes to download locations

from pydantic import BaseModel, Field
from typing import List, Optional
import requests
from requests.exceptions import HTTPError
import json, time
import aiohttp
import asyncio
from aiohttp import ClientTimeout
import pandas as pd
import os
import ast
import logging

class File(BaseModel):
    name: str
    format: Optional[str] = None # Maybe not optional??
    url: Optional[str] = None

class Project(BaseModel):
    projectTitle: Optional[List[Optional[str]]] = None
    laboratory: Optional[List[Optional[str]]] = None

class Sample(BaseModel):
    organ: Optional[List[Optional[str]]] = None
    disease: Optional[List[Optional[str]]] = None

class Hit(BaseModel):
    projects: List[Project]
    samples: Optional[List[Sample]] = None
    files: List[File]

#dict
# filters available on HCA
supported_fields = [
    "accessions", "aggregateLastModifiedDate", "aggregateSubmissionDate", "aggregateUpdateDate",
    "assayType", "biologicalSex", "bionetworkName", "bundleUuid", "bundleVersion", "cellCount",
    "cellLineType", "contactName", "contentDescription", "dataUseRestriction", "developmentStage",
    "donorCount", "donorDisease", "duosId", "effectiveCellCount", "effectiveOrgan", "entryId",
    "fileFormat", "fileId", "fileName", "fileSize", "fileSource", "fileVersion", "genusSpecies",
    "institution", "instrumentManufacturerModel", "isIntermediate", "isTissueAtlasProject",
    "laboratory", "lastModifiedDate", "libraryConstructionApproach", "matrixCellCount",
    "modelOrgan", "modelOrganPart", "nucleicAcidSource", "organ", "organPart", "organismAge",
    "organismAgeRange", "pairedEnd", "preservationMethod", "project", "projectDescription",
    "projectEstimatedCellCount", "projectId", "projectTitle", "publicationTitle", "sampleDisease",
    "sampleEntityType", "sampleId", "selectedCellType", "sourceId", "sourceSpec", "specimenDisease",
    "specimenOrgan", "specimenOrganPart", "submissionDate", "tissueAtlas", "updateDate",
    "workflow", "accessible"
]

# ask user
path = input("\nPath to a config file for filtering (press Enter for manual input): ").strip()

# init dict for filter
filters = {}

if path and os.path.isfile(path):
    with open(path) as f:
        # set json as filter
        filters = json.load(f)
    print(f"\nLoaded filters from {path}")
else:
    # print available fields
    print("\nSupported fields:\n" + ", ".join(supported_fields))

    # get input for which fields to select
    chosen = input("\nEnter fields to filter (comma separated, e.g. fileFormat,genusSpecies,isIntermediate, fileSource): ").strip()
    # sanitize
    selected = [f.strip() for f in chosen.split(",") if f.strip()]

    for f in selected:
        # ensure the provided field is present in HCA api
        if f not in supported_fields:
            print(f"Skipping unknown field: {f}")
            continue

        # get input for getting the values for the fields selected in the step b4
        val = input(f"Enter value(s) for '{f}' (comma separated, True/False if applicable): ").strip()

        if not val:
            continue

        try: 
            # ast function used here for passing in the User input directly as a bool. ≠ string. 
            filters[f] = {"is": [ast.literal_eval(val)]}
        except Exception: 
            filters[f] = {"is": [v.strip() for v in val.split(",") if v.strip()]}

# logging.info
print("\nUsing filters:")
print(json.dumps(filters, indent=2))

size=100
catalog = "dcp54"

base_url = "https://service.azul.data.humancellatlas.org/index/files"

def fetch_all_pages(catalog=catalog, filters=filters, size=50):
    params = {"catalog": catalog, "filters": json.dumps(filters), "size": size}
    hits, url, page = [], base_url, 1 # empty list to store evrytg

    while url:
        print(f"page {page}")
        r = requests.get(url, params=params if url == base_url else None)
        # stat at this point
        ## logging.error here
        if r.status_code != 200:
            print("error", r.status_code, r.text[:200])
            break

        data = r.json()
        if page == 1:
            facets = data.get("termFacets")
        hits += data.get("hits", [])
        url = data.get("pagination", {}).get("next")
        params = None
        page += 1

        time.sleep(0.2)

    return hits, facets

raw_hits, facets = fetch_all_pages()
#raw_hits = fetch_all_pages()
hits = [Hit(**h) for h in raw_hits]
# logging.info
print("pydantic validated:", len(hits))

formats = facets["fileFormat"]["terms"]
formats_pd = pd.DataFrame(formats)

# accounting for multiple samples in hits

rows = []
for hit in hits:
    for project in hit.projects or []:
        project_title = ", ".join(project.projectTitle or [])
        labs = project.laboratory or []
        lab_names = "; ".join(str(lab) for lab in labs if lab) or "N/A"

        # collect all organs/diseases across samples
        samples = hit.samples or []
        if samples:
            organs = set()
            diseases = set()
            for s in samples:
                if getattr(s, "organ", None):
                    organs.update([str(o) for o in s.organ if o])
                if getattr(s, "disease", None):
                    diseases.update([str(d) for d in s.disease if d])
            organ = ", ".join(organs) if organs else "N/A"
            disease = ", ".join(diseases) if diseases else "N/A"
        else:
            organ = "N/A"
            disease = "N/A"

# Sasb4; in future add other metadata here
        for file in hit.files or []:
            rows.append({
                "Project": project_title,
                "Project lab": lab_names,
                "Organ": organ,
                "Disease": disease,
                "File": file.name,
                "Url": file.url
            })

df = pd.DataFrame(rows)
#print(df)

## Extract Url
urls = df["Url"].tolist()

# test with metadata loom file
from aiohttp import ClientTimeout
timeout = ClientTimeout(total=None)
headers = {"User-Agent": "Mozilla/5.0"}

# download func; In 1 mb chunks
async def download_file(url, filename, chunk_size=1024*1024, retries=3):
    for i in range(retries):
        try:
            connector = aiohttp.TCPConnector(limit=10)
            async with aiohttp.ClientSession(timeout=timeout, headers=headers, connector=connector) as session:
                async with session.get(url) as response:
                    if response.status != 200:
                        raise Exception(f"Failed: {response.status}")
                    
                    with open(filename, "wb") as f:
                        async for chunk in response.content.iter_chunked(chunk_size):
                            f.write(chunk)
            # logging.info
            print("Download complete:", filename)
            return
        except Exception as e:
            # logging.warning
            print(f"Retry {i+1} failed:", e)
            time.sleep(2 ** i) # adjust

# choose the save location
save_dir = input("\nDirectory to save downloads (press Enter for default: data_raw/HCA_downloads): ").strip()
if not save_dir:
    save_dir = "data_raw/HCA_downloads"

os.makedirs(save_dir, exist_ok=True)

# allow user to specify
pd.set_option('display.max_rows', None)
print(df[["File", "Organ", "Disease"]].reset_index())
#pd.reset_option("display.max_rows")
choice = int(input("\nWhich file to download? Enter the index: "))

url = df.loc[choice, "Url"]
filename = os.path.join(save_dir, df.loc[choice, 'File'] or f"file_{choice+1}")

df.to_csv(os.path.join(save_dir, "metadata.csv"), index=False)

# logging.info
print(f"\nDownloading: {filename}")

asyncio.run(download_file(url, filename))
# logging.info
print("\nSuccess.")

# To play nice with workflow.py
with open("data_raw/_last_downloaded.txt", "w") as f:
    f.write(filename)

#add logging??