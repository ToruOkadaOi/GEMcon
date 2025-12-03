import subprocess

def test_help_runs():
    result = subprocess.run(["python", "scripts/run.py", "--help"], capture_output=True)
    assert result.returncode == 0