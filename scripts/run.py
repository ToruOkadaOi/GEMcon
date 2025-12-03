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

# idea is to have seperate conda envs

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--branch", choices=["annotate_cells", "metabolic"], required=True) # simpler names in CLI? -a -b indicative of seq. to follow?
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
        # run(["python", "src/api_hca_userinp.py"])
        run(["bash", "src/run_in_scanpy.sh", "src/api_hca_userinp.py"])
        # adjust this to the fetch script outputs
        with open("data_raw/_last_downloaded.txt") as f: # using _last_downloaded.txt as identifier
            input_file = f.read().strip()
        #print(f"\n Using fetched file: {input_file}")
        console.rule(f"[bold cyan]Using fetched file: {input_file}")

    # branch A
    if args.branch == "annotate_cells":
        # normalize per cell, apply log transform to pass onto celltypist. Clustering could be added, but CellTypist can do clustering on its own
        run(["bash", "src/run_in_scanpy.sh", "src/scanpy_norm.py", "--input", input_file])
        run(["bash", "src/run_in_scanpy.sh", "src/annotate_celltypes.py"])
        return

    # branch B
    if args.branch == "metabolic":
        # 1. norm_pooling.py = Normalize by pooling all celltypes or extract cell types # psuedo-bulk
        cmd = ["bash", "src/run_in_scanpy.sh", "src/norm_pooling.py"]
        if input_file: cmd += ["--input", input_file]
        if celltype is not None:
            cmd += ["--celltype", celltype]
        run(cmd)

        # 2. genetoensembl.py = To change the HGNC symbol to Ensembl ids (compatible with HumanGEM). #TODO : EntrezId for Recon3D?
        cmd = ["bash", "src/run_in_scanpy.sh", "src/genetoensembl.py"] # I need to run this before context algorithms
        if gtf: cmd += ["--gtf", gtf]
        run(cmd) 

        # 3. Any of the algorithm. GIMME initially #TODO: specifiable in config
        cmd = ["bash", "src/run_in_cplex.sh", "src/gimme.py"]
        model_path = cfg.get("model")
        if model_path:
            cmd += ["--model", model_path]
        run(cmd)
        #run(["bash", "src/run_in_cplex.sh", "src/tinit.py"])
        return

if __name__ == "__main__":
    main()