#!/usr/bin/env python3
import argparse
import subprocess
import os

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

    #1 - count matrices already available
    if args.input:
        input_file = args.input
        if not os.path.exists(input_file):
            raise FileNotFoundError(input_file)
        print(f"\n>>> using user file: {input_file}")

    #2 no matrices i.e to fetch from databases
    else:
        run(["python", "scripts/api_hca_userinp.py"])
        # adjust this to whatever your fetch script outputs
        with open("data_raw/_last_downloaded.txt") as f:
            input_file = f.read().strip()
        print(f"\n>>> using fetched file: {input_file}")

    # branch A
    if args.branch == "celltype_annotated":
        run(["python", "scripts/scanpy_norm.py", input_file])
        run(["python", "scripts/annotate_celltypes.py"])
        return

    # branch B
    if args.branch == "metabolic":
        run(["bash", "scripts/run_in_scanpy.sh", "scripts/norm_pre.py", input_file])
        run(["bash", "scripts/run_in_cplex.sh", "scripts/gimme.py"])
        return

if __name__ == "__main__":
    main()