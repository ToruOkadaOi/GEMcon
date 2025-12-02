__author__ = "Aman Nalakath"
__description__ = ""

#TODO: logging + testing; prints --> rich;

import pandas as pd
import numpy as np
import cobra
import re
import os
import argparse
import yaml

# Troppo imports
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

p = argparse.ArgumentParser()
p.add_argument("--expr", help="Path to expression data csv") # or just -e instead?? # go cobra-cli subcommands
p.add_argument("--model", help="Path to the sbml model eg Human-GEM.xml")
args = p.parse_args()

# Load the config if existing
cfg = {}
if os.path.exists("config.yaml"):
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f) or {}

model_cfg = cfg.get('model_config', {})

# get the pattern from the config
pattern = model_cfg.get('alt_transcript_pattern')
pattern = model_cfg.get('alt_transcript_pattern', '__COBAMPGPRDOT__[0-9]{1}') # keep def. here?? or ask inter.?
patt = re.compile(pattern)
replace_alt_transcripts = lambda x: patt.sub('', x)

# for recon3d; # TODO: decide on the proper pattern # patt = re.compile(r'_AT[0-9]+$') ## replace_alt_transcripts = lambda x: patt.sub('', x)

# if model (eg. Human-GEM or Recon3D) is specified as cli arg
if args.model:
    model_path = args.model.strip()
else:
    model_path = cfg.get("model")

if not model_path:
    model_path = input("\nProvide full path to the metabolic model (.xml): ").strip()

# model_path = input("\nProvide full path to the SBML model (.xml): ").strip()
if not os.path.exists(model_path):
    raise FileNotFoundError(model_path)

if args.expr:   # add argparse in here too
    expr_path = args.expr.strip()
else:
    files = [f for f in os.listdir("data_processed") if f.startswith("expression_data_") and f.endswith("_gencode.csv")]

    if not files:
        expr_path = input("\nNo files found, please provide the abs. path to an expression .csv with gencode symbols: ").strip()
    else:
        paths = [os.path.join("data_processed", f) for f in files]
        expr_path = max(paths, key=os.path.getmtime)
        print(f"Using the last made file: {expr_path}")

if not os.path.exists(expr_path):
    raise FileNotFoundError(expr_path)

# basename
base = os.path.splitext(os.path.basename(expr_path))[0]

model = cobra.io.read_sbml_model(model_path)

expression_data = pd.read_csv(expr_path, index_col=0)
print(expression_data.head())

# Transpose to pass to Troppo TabularReader()
expression_data_transposed = expression_data.T

# Check
print("\nTransposed shape:", expression_data_transposed.shape)
print(expression_data_transposed.head())

omics_container = TabularReader(path_or_df=expression_data_transposed, 
                                nomenclature='gene',
                                omics_type='transcriptomics').to_containers()

# most single cell data are one sample only. TODO: Have to fig. out cases for multi samples esp. in case of Bulk RNASeq
single_sample = omics_container[0]

model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999,
                                      gpr_gene_parse_function=replace_alt_transcripts)

data_map = single_sample.get_integrated_data_map(model_reader=model_wrapper.model_reader,
                                                 and_func=min, or_func=sum)

print(list(data_map.get_scores().items())[:10])

def score_apply(reaction_map_scores):
    return {k:0  if v is None else v for k, v in reaction_map_scores.items()}

continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
scores = continuous_integration.integrate(data_map=data_map)
model = model_wrapper.model_reader.model

objective_rxn = model_cfg.get('objective_reaction')

# if reaction.id (eg. MAR02388) is not specified in config, ask user interactively
if not objective_rxn:
    objective_rxn = input("\nEnter objective reaction ID (e.g., MAR02388): ").strip()

idx_objective = model_wrapper.model_reader.r_ids.index(objective_rxn)

# get the params from config or ask interactively; should I bother with hvn cli args?
gimme_cfg = cfg.get('gimme', {})

# Get parameters from config or ask user # TODO: revisit the fallback defaults
obj_frac = gimme_cfg.get('obj_frac')
if obj_frac is None:
    obj_frac = float(input("\nEnter obj_frac (default 0.1): ").strip() or 0.1)

flux_threshold = gimme_cfg.get('flux_threshold')
if flux_threshold is None:
    flux_threshold = float(input("Enter flux_threshold (default 0.25): ").strip() or 0.25)

solver = gimme_cfg.get('solver')
if not solver:
    solver = input("Enter solver (default CPLEX; make sure it is installed properly!): ").strip() or 'CPLEX' 

properties = GIMMEProperties(
    exp_vector=[v for k, v in scores.items()],
    obj_frac=obj_frac,
    objectives=[{idx_objective: 1}],
    preprocess=True,
    flux_threshold=flux_threshold,
    solver=solver,
    reaction_ids=model_wrapper.model_reader.r_ids,
    metabolite_ids=model_wrapper.model_reader.m_ids
)

# Run the GIMME algorithm.
gimme = GIMME(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

gimme_run = gimme.run()

# SAVE
# get the save path from config
mod_dir = cfg.get('output', {}).get('models_dir')
if not mod_dir:
    mod_dir = input("Enter output directory for saving the models (default './models'): ").strip() or './models'

os.makedirs(mod_dir, exist_ok=True) # raise err?

pd.DataFrame(gimme_run).to_csv(f'{mod_dir}/gimme_output_{base}.tsv', sep='\t')

print(pd.DataFrame(gimme_run).head()) # TODO:rich + context

selected_reactions = [model.reactions[i] for i in gimme_run]
print(pd.DataFrame(selected_reactions).head())

ctx_model = model.copy()

selected_ids = [r.id for r in selected_reactions]
to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]

ctx_model.remove_reactions(to_remove, remove_orphans=True)

print(f"The number of reactions in the base model is: {len(model.reactions)} "
      f"& the number of reactions in the extracted model is: {len(ctx_model.reactions)}")

# export as .xml. # TODO:db integration?
cobra.io.write_sbml_model(ctx_model, f"{mod_dir}/gimme_context_specific_model_{base}.xml")

# printing the objective val.
solution = ctx_model.optimize()
print('The objective value is: ', solution.objective_value)