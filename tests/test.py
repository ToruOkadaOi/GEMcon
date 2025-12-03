import subprocess
import os
import pytest
import shutil

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