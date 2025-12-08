from prefect import flow, task
import subprocess
import os
import yaml
from typing import Optional
from rich import print
from rich.console import Console

console = Console()

@task
def run_cmd(cmd: list):
    # print("running:", cmd)
    console.rule(f"[bold cyan]running: {cmd}")
    subprocess.run(cmd, check=True)

@task
def load_config():
    try:
        with open("config.yaml") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        return {}

@task
def resolve_input_file(args_input, cfg):
    # user-provided
    if args_input:
        if not os.path.exists(args_input):
            raise FileNotFoundError(args_input)
        return args_input

    # config-provided
    if cfg.get("input_file"):
        if not os.path.exists(cfg["input_file"]):
            raise FileNotFoundError(cfg["input_file"])
        return cfg["input_file"]

    # fetch
    run_cmd(["bash", "src/run_in_scanpy.sh", "src/api_hca_userinp.py"])
    with open("data_raw/_last_downloaded.txt") as f:
        return f.read().strip()

@task
def branch_annotate(input_file):
    run_cmd(["bash", "src/run_in_scanpy.sh", "src/scanpy_norm.py", "--input", input_file])
    run_cmd(["bash", "src/run_in_scanpy.sh", "src/annotate_celltypes.py"])

@task
def branch_metabolic(input_file, cfg):
    celltype = cfg.get("celltype")
    gtf      = cfg.get("gtf")

    cmd = ["bash", "src/run_in_scanpy.sh", "src/norm_pooling.py", "--input", input_file]
    if celltype is not None:
        cmd += ["--celltype", celltype]
    run_cmd(cmd)

    cmd = ["bash", "src/run_in_scanpy.sh", "src/genetoensembl.py"]
    if gtf:
        cmd += ["--gtf", gtf]
    run_cmd(cmd)

    cmd = ["bash", "src/run_in_cplex.sh", "src/gimme.py"]
    model = cfg.get("model")
    if model:
        cmd += ["--model", model]
    run_cmd(cmd)

@flow
def main_flow(branch: str, input_file: Optional[str] = None):
    cfg = load_config()
    resolved_input = resolve_input_file(input_file, cfg)

    if branch == "annotate":
        branch_annotate(resolved_input)
    elif branch == "metabolic":
        branch_metabolic(resolved_input, cfg)

if __name__ == "__main__":
    # inside the if block, as this file would be imported later d.line
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--branch", choices=["annotate", "metabolic"], required=True)
    p.add_argument("--input")
    args = p.parse_args()
    main_flow(args.branch, args.input)