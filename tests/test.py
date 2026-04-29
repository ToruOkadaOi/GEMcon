import subprocess
import os
import pytest
import shutil
import requests
import json
import time

# import aiohttp
# import asyncio
from aiohttp import ClientTimeout
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import yaml


# ---
## --help
def test_help_runs():
    result = subprocess.run(
        ["python", "-m", "scripts.flow", "--help"], capture_output=True
    )
    assert result.returncode == 0


# ---
## config.yaml
### reads config f # Remove maybe
def test_reads_config_yaml():
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f)
    assert isinstance(cfg, dict)


# ---
## input
@pytest.mark.skipif(
    os.getenv("CI") == "true", reason="CI skip"
)  # this test fails on git ci # TODO: workarounds?
def test_fetches_data_when_no_input():
    # Temporarily move config.yaml if it exists
    config_exists = os.path.exists("config.yaml")
    if config_exists:
        shutil.move("config.yaml", "config.yaml.backup")
    try:
        result = subprocess.run(
            [
                "python",
                "-m", "scripts.flow",
                "--branch",
                "transcriptomic",
                "--task",
                "annotate",
            ],
            capture_output=True,
            text=True,
            timeout=20,
        )
        output = result.stdout + result.stderr
        # Should try to run the api script
        assert (
            "api_hca_userinp.py" in output or "run_in_scanpy.sh" in output
        )  # indicate that it is running/on the api call code

    finally:
        # Restore config
        if config_exists:
            shutil.move("config.yaml.backup", "config.yaml")


# ---
## HCA api ### test reachability
base_url = "https://service.azul.data.humancellatlas.org/index/files"


def get_latest_catalog():
    """Fetch the latest DCP catalog from HCA API."""
    import re

    try:
        r = requests.get(
            "https://service.azul.data.humancellatlas.org/index/catalogs", timeout=10
        )
        r.raise_for_status()
        catalogs = r.json().get("catalogs", {})
        dcp_catalogs = [k for k in catalogs if k.startswith("dcp")]
        # Prefer non-it catalogs (dcp57 over dcp57-it) - they have more file formats
        base_catalogs = [k for k in dcp_catalogs if "-it" not in k]
        if base_catalogs:

            def get_num(name):
                m = re.search(r"dcp(\d+)", name)
                return int(m.group(1)) if m else 0

            base_catalogs.sort(key=get_num)
            return base_catalogs[-1]

        def get_num(name):
            m = re.search(r"dcp(\d+)", name)
            return int(m.group(1)) if m else 0

        dcp_catalogs.sort(key=get_num)
        return dcp_catalogs[-1] if dcp_catalogs else "dcp57"
    except Exception:
        return "dcp57"


def test_api_reachable():
    filters = {"fileFormat": {"is": ["loom"]}}
    catalog = get_latest_catalog()
    params = {"catalog": catalog, "filters": json.dumps(filters), "size": 1}
    r = requests.get(base_url, params=params)
    assert r.status_code == 200


## uniprot api ### geckopy #### test reachability

base_post = "https://rest.uniprot.org/idmapping/run"
base_get = "https://rest.uniprot.org/idmapping/results/"


def test_uniprot():
    data = {"from": "Ensembl_Protein", "to": "UniProtKB", "ids": "ENSP00000370010"}
    r = requests.post(base_post, data=data)
    assert r.status_code == 200
    assert "jobId" in r.json()
    job_id = r.json()["jobId"]

    # ---
    time.sleep(3)

    result = requests.get(f"{base_get}{job_id}").json()
    assert "results" in result
    assert len(result["results"]) > 0


# ---
## PaxDB ### test reachability # write better test
timeout = ClientTimeout(total=None)
headers = {"User-Agent": "Mozilla/5.0"}


def get_links():
    base = "https://pax-db.org"
    listing_url = f"{base}/downloads/3.0/datasets/9606/"

    r = requests.get(listing_url)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    links = []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if href.endswith(".txt"):
            links.append(urljoin(listing_url, href))

    return links


def test_links():
    links = get_links()
    assert len(links) > 0, "No links found"

    r = requests.head(links[0], headers=headers, timeout=10)
    assert r.status_code == 200


# Test that the reactions are getting a score according to the exp.csv
# Test extracted model # CI ignore or dummy
# Test loom and h5ad produce the same result
# Any data losses in the branches
# Algo specific tests
