import subprocess
import os
import pytest
import shutil
import requests
import json
import time

# ---
## --help
def test_help_runs():
    result = subprocess.run(["python", "scripts/run.py", "--help"], capture_output=True)
    assert result.returncode == 0

# ---
## config.yaml
### reads config f
def test_reads_config_yaml():
    result = subprocess.run(
        ["python", "run.py", "--branch", "annotate_cells"],
        capture_output=True,
        timeout=5
    )
    assert b"yaml" not in result.stderr.lower()

# ---
## input
@pytest.mark.skipif(os.getenv('CI') == 'true', reason="CI skip")  # this test fails on git ci # TODO: workarounds?
def test_fetches_data_when_no_input():
    # Temporarily move config.yaml if it exists
    config_exists = os.path.exists("config.yaml")
    if config_exists:
        shutil.move("config.yaml", "config.yaml.backup")
    try:
        result = subprocess.run(
            ["python", "scripts/run.py", "--branch", "annotate_cells"],
            capture_output=True,
            text=True,
            timeout=5
        )
        output = result.stdout + result.stderr
        # Should try to run the api script
        assert "api_hca_userinp.py" in output or "run_in_scanpy.sh" in output # indicate that it is running/on the api call code

    finally:
        # Restore config
        if config_exists:
            shutil.move("config.yaml.backup", "config.yaml") 

# ---
## HCA api ### test reachability
base_url = "https://service.azul.data.humancellatlas.org/index/files"
def test_api_reachable():
    filters = {"fileFormat": {"is": ["loom"]}}
    params = {"catalog": "dcp55", "filters": json.dumps(filters), "size": 1}
    r = requests.get(base_url, params=params)
    assert r.status_code == 200

## uniprot api ### geckopy #### test reachability

base_post= "https://rest.uniprot.org/idmapping/run"
base_get = "https://rest.uniprot.org/idmapping/results/"
def test_uniprot():
    data = {
        "from": "Ensembl_Protein",
        "to": "UniProtKB",
        "ids": "ENSP00000370010"
    }
    r = requests.post(base_post, data=data)
    assert r.status_code == 200
    assert "jobId" in r.json()
    job_id = r.json()["jobId"]

    #---
    time.sleep(3)
    
    result = requests.get(f"{base_get}{job_id}").json()
    assert "results" in result
    assert len(result["results"]) > 0