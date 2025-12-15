# same as run.py but with prefect orchestrator
from prefect import flow, task
import subprocess
import os
import yaml
from typing import Optional
from rich import print
from rich.console import Console

console = Console()

# usuable context algo.s
algorithms = {
    "gimme":       ("gimme.py", "cplex"),
    "tinit":       ("tinit.py", "cplex"),
    "fastcore":    ("fastcore_beta.py", "cplex"),
    "geckopy":     ("gecko_pipeline.py", "gecko"),
    "imat":        ("imat_beta.py", "cplex"),
#    "cadre":       ("pymCADRE.py", "cplex"), # sep. env also needed
}

@task
def run_cmd(cmd: list):
    # print("running:", cmd)
    console.rule(f"[bold cyan]running: {cmd}")
    subprocess.run(cmd, check=True)

def run_algo(script_name: str, env: str, args: list = None):
    cmd = ["bash", f"src/run_in_{env}.sh", f"src/{script_name}"] # TODO: better way to handle conda envs?? single bash file???
    if args:
        cmd += args
    run_cmd(cmd)

@task
def load_config():
    try:
        with open("config.yaml") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        return {}

@task
def resolve_input_file(args_input, cfg):
    # cli
    if args_input:
        if not os.path.exists(args_input):
            raise FileNotFoundError(args_input)
        return args_input

    # config
    if cfg.get("input_file"):
        if not os.path.exists(cfg["input_file"]):
            raise FileNotFoundError(cfg["input_file"])
        return cfg["input_file"]

    # fetch
    run_cmd(["bash", "src/run_in_scanpy.sh", "src/api_hca_userinp.py"])
    with open("data/data_raw/_last_downloaded.txt") as f:
        return f.read().strip()

@task
def branch_transcriptomic(input_file, cfg, task_type: str, algo: str = "gimme"):
    if task_type == "annotate":     # /home/biodata/aman/docs/explanation.md

        run_cmd(["bash", "src/run_in_scanpy.sh", "src/scanpy_norm.py", "--input", input_file])
        run_cmd(["bash", "src/run_in_scanpy.sh", "src/annotate_celltypes.py"])

    elif task_type == "metabolic":      # /home/biodata/aman/docs/explanation.md

        celltype = cfg.get("celltype")
        gtf      = cfg.get("gtf")
        model    = cfg.get("model")

        # 1. normalizing
        cmd = ["bash", "src/run_in_scanpy.sh", "src/norm_pooling.py", "--input", input_file]
        if celltype is not None:
            cmd += ["--celltype", celltype]
        run_cmd(cmd)

        # 2. id converting
        cmd = ["bash", "src/run_in_scanpy.sh", "src/genetoensembl.py"]
        if gtf:
            cmd += ["--gtf", gtf]
        run_cmd(cmd)

        # 3. algorithm
        script, env = algorithms[algo]
        args = ["--model", model] if model else [] # internal resolving instead??
        run_algo(script, env, args=args)

@task
def branch_proteomic(cfg, algo: str = "geckopy"):
    script, env = algorithms[algo]
    run_algo(script, env, args=None)       


@flow
def main_flow(branch: str, input_file: Optional[str] = None, task: Optional[str] = None, algo: Optional[str] = None):
    cfg = load_config()
    
    if branch == "transcriptomic":
        
        algo = algo or "gimme" # defau.

        resolved_input = resolve_input_file(input_file, cfg)
        branch_transcriptomic(resolved_input, cfg, task, algo)

    elif branch == "proteomic":

        algo = algo or "geckopy" # defau.

        branch_proteomic(cfg, algo)

if __name__ == "__main__":
    # inside the if block, as this file would be imported later d.line
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--branch", choices=["transcriptomic", "proteomic"], required=True)
    p.add_argument("--task", choices=["annotate", "metabolic"])
    p.add_argument("--input")
    p.add_argument("--algo", default=None, help="Algorithm for metabolic modeling (default - gimme for transcriptomic, geckopy for proteomic)")
    args = p.parse_args()

    # check
    if args.branch == "transcriptomic" and not args.task:
        p.error("--task is required when --branch is transcriptomic")

    main_flow(args.branch, args.input, args.task, args.algo)