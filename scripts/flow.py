import subprocess
import os
import yaml
from typing import Optional
from rich import print
from rich.console import Console

from src.algorithms import (
    GimmeAlgorithm,
    TinitAlgorithm,
    CordaAlgorithm,
    RiptideAlgorithm,
    FastcoreAlgorithm,
    ImatAlgorithm,
)

console = Console()

# algorithm name -> class
algorithms = {
    "gimme": GimmeAlgorithm,
    "tinit": TinitAlgorithm,
    "corda": CordaAlgorithm,
    "fastcore": FastcoreAlgorithm,
    "imat": ImatAlgorithm,
    "riptide": RiptideAlgorithm,
}


def run_cmd(cmd: list):
    console.rule(f"[bold cyan]running: {cmd}")
    subprocess.run(cmd, check=True)


def load_config():
    try:
        with open("config.yaml") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        return {}


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


def branch_transcriptomic(input_file, cfg, task_type: str, algo: str = "gimme"):
    if task_type == "annotate":
        run_cmd(
            ["bash", "src/run_in_scanpy.sh", "src/scanpy_norm.py", "--input", input_file]
        )
        run_cmd(["bash", "src/run_in_scanpy.sh", "src/annotate_celltypes.py"])

    elif task_type == "metabolic":
        celltype = cfg.get("celltype")
        gtf = cfg.get("gtf")
        model = cfg.get("model")

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
        AlgoClass = algorithms[algo]
        args = type("Args", (), {"model": model, "expr": None})()
        algo_instance = AlgoClass(args=args)
        algo_instance.run()


def branch_proteomic(cfg, algo: str = "geckopy"):
    # geckopy still needs its own env
    run_cmd(["bash", "src/run_in_gecko.sh", "src/gecko_pipeline.py"])


def main_flow(
    branch: str,
    input_file: Optional[str] = None,
    task: Optional[str] = None,
    algo: Optional[str] = None,
):
    cfg = load_config()

    if branch == "transcriptomic":
        algo = algo or "gimme"
        resolved_input = resolve_input_file(input_file, cfg)
        branch_transcriptomic(resolved_input, cfg, task, algo)

    elif branch == "proteomic":
        algo = algo or "geckopy"
        branch_proteomic(cfg, algo)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--branch", choices=["transcriptomic", "proteomic"], required=True)
    p.add_argument("--task", choices=["annotate", "metabolic"])
    p.add_argument("--input")
    p.add_argument("--algo", default=None)
    args = p.parse_args()

    if args.branch == "transcriptomic" and not args.task:
        p.error("--task is required when --branch is transcriptomic")

    main_flow(args.branch, args.input, args.task, args.algo)