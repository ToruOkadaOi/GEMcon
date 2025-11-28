#!/usr/bin/env python3
__author__ = "Aman Nalakath"

import argparse
import subprocess
import os
import yaml
from rich import print
from rich.console import Console
from rich.table import Table

console = Console()

def run(cmd):
    print("\n>>>", " ".join(cmd))
    subprocess.run(cmd, check=True)

# idea is to have seperate conda envs for each modules
# def run_api(cmd):
#     subprocess.run(["conda", "run", "-n", "hca_api"] + cmd, check=True)
# def run_scanpy(cmd):
#     subprocess.run(["conda", "run", "-n", "scanpy_legacy"] + cmd, check=True)
# def run_cplex(cmd):
#     subprocess.run(["conda", "run", "-n", "cplex_aman"] + cmd, check=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--branch", choices=["celltype_annotated", "metabolic"], required=True)
    p.add_argument("--input", help="loom/h5ad file. If not given, it will be fetched.")
    args = p.parse_args()

    try:
        with open("config.yaml") as f:
            cfg = yaml.safe_load(f) # config var
    except FileNotFoundError:
        cfg = {}

    # parameters in .yaml
    input_file = args.input or cfg.get("input_file")
    celltype   = cfg.get("celltype")
    gtf        = cfg.get("gtf")


    #1 - count matrices already available
    if args.input: # here passed as --input para
        input_file = args.input
        if not os.path.exists(input_file):
            raise FileNotFoundError(input_file)
        #print(f"\n Using user file: {input_file}")
        console.rule(f"[bold cyan]Using user file: {input_file}")

    elif input_file:   # loaded from config instead
        if not os.path.exists(input_file):
            raise FileNotFoundError(input_file)
        #print(f"\n Using file specified in config: {input_file}")
        console.rule(f"[bold cyan]Using file specified in config: {input_file}")

    #2 no matrices i.e to fetch from databases
    else:
        # run(["python", "scripts/api_hca_userinp.py"])
        run(["bash", "scripts/run_in_scanpy.sh", "scripts/api_hca_userinp.py"])
        # adjust this to the fetch script outputs
        with open("data_raw/_last_downloaded.txt") as f: # using _last_downloaded.txt as identifier
            input_file = f.read().strip()
        #print(f"\n Using fetched file: {input_file}")
        console.rule(f"[bold cyan]Using fetched file: {input_file}")

    # branch A
    if args.branch == "celltype_annotated":
        run(["bash", "scripts/run_in_scanpy.sh", "scripts/scanpy_norm.py", input_file])
        run(["bash", "scripts/run_in_scanpy.sh", "scripts/annotate_celltypes.py"])
        return

    # branch B
    if args.branch == "metabolic":
        # 1. norm_pre.py = Normalize by pooling all celltypes or extract cell types
        cmd = ["bash", "scripts/run_in_scanpy.sh", "scripts/norm_pre.py"]
        if input_file: cmd += ["--input", input_file]
        if celltype is not None:
            cmd += ["--celltype", celltype]
        run(cmd)

        # 2. genetoensembl.py = To change the HGNC symbol to Ensembl ids (compatible with HumanGEM). #TODO : EntrezId for Recon3D
        cmd = ["bash", "scripts/run_in_scanpy.sh", "scripts/genetoensembl.py"] # I need to run this before context algorithms
        if gtf: cmd += ["--gtf", gtf]
        run(cmd) 

        # 3. Any of the algorithm. GIMME initially
        cmd = ["bash", "scripts/run_in_cplex.sh", "scripts/gimme.py"]
        model_path = cfg.get("model")
        if model_path:
            cmd += ["--model", model_path]
        run(cmd)
        #run(["bash", "scripts/run_in_cplex.sh", "scripts/tinit.py"])
        return

if __name__ == "__main__":
    main()