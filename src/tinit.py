import pandas as pd
import numpy as np
import cobra
import re
import os
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
import argparse
import yaml

p = argparse.ArgumentParser()
p.add_argument("--expr", help="Path to expression data csv")
p.add_argument("--model", help="Path to the sbml model eg Human-GEM.xml")
args = p.parse_args()

# Load config if existing
cfg = {}

if os.path.exists("config.yaml"):
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f) or {}

model_cfg = cfg.get('model_config', {})
# tinit config part
tinit_cfg = cfg.get('tinit', {})

# get the pattern from the config
pattern = model_cfg.get('alt_transcript_pattern')
pattern = model_cfg.get('alt_transcript_pattern', '__COBAMPGPRDOT__[0-9]{1}') # keep def. here?? or ask inter.? ## update: fallback to Human-GEM pattern
patt = re.compile(pattern)
replace_alt_transcripts = lambda x: patt.sub('', x)

# detect model type
if args.model:
    model_path = args.model.strip()
else:
    model_path = cfg.get("model")

if model_path and ('_AT' in pattern or 'recon3d' in model_path.lower() or 'Recon3D' in model_path):
    model_type = 'recon3d'
    expr_suffix = '_recon3d.csv'
else:
    model_type = 'human-gem'
    expr_suffix = '_gencode.csv'

if not model_path:
    model_path = input("\nProvide full path to the metabolic model (.xml): ").strip()

if not os.path.exists(model_path):
    raise FileNotFoundError(model_path)

# expression file lookup
if args.expr:
    expr_path = args.expr.strip()
else:
    files = [f for f in os.listdir("data/data_processed") 
             if f.startswith("expression_data_") and f.endswith(expr_suffix)]
    
    if not files:
        expr_path = input(f"\nNo files found, provide path to expression .csv ({model_type} format): ").strip()
    else:
        paths = [os.path.join("data/data_processed", f) for f in files]
        expr_path = max(paths, key=os.path.getmtime)
        print(f"using {expr_path}")

if not os.path.exists(expr_path):
    raise FileNotFoundError(expr_path)

# Get basename for output files
base = os.path.splitext(os.path.basename(expr_path))[0]

model = cobra.io.read_sbml_model(model_path)
expression_data = pd.read_csv(expr_path, index_col=0)
print(expression_data.head())

expression_data_transposed = expression_data.T

nomenclature = tinit_cfg.get('nomenclature', 'ensembl_gene_id') # default to ensembl here
omics_container = TabularReader(path_or_df=expression_data_transposed, 
                                nomenclature=nomenclature,
                                omics_type='transcriptomics').to_containers()

single_sample = omics_container[0]

model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999,
                                      gpr_gene_parse_function=replace_alt_transcripts)

data_map = single_sample.get_integrated_data_map(model_reader=model_wrapper.model_reader,
                                                 and_func=min, or_func=sum)

def score_apply(reaction_map_scores):
    return {k:0  if v is None else v for k, v in reaction_map_scores.items()}

continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
scores = continuous_integration.integrate(data_map=data_map)
print(f'\nScores are:\n {scores.items()}')

essential_rxn_ids = tinit_cfg.get('essential_reactions')
if not essential_rxn_ids:
    user_input = input("\nEnter essential reaction ids (eg,MAR13082 in Human-GEM): ").strip()
    essential_rxn_ids = [r.strip() for r in user_input.split(',') if r.strip()]

essential_reactions = [model_wrapper.model_reader.r_ids.index(rid) for rid in essential_rxn_ids]

solver = tinit_cfg.get('solver')
if not solver:
    solver = input("\nEnter solver (default CPLEX): ").strip() or 'CPLEX'

properties = tINITProperties(
    reactions_scores=[v for k, v in scores.items()], 
    solver=solver, 
    essential_reactions=essential_reactions
)

tinit = tINIT(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)
# run
model_tinit = tinit.run()

mod_dir = cfg.get('output', {}).get('models_dir')
if not mod_dir:
    mod_dir = input("Enter output directory for saving the models (default './models'): ").strip() or './models'

os.makedirs(mod_dir, exist_ok=True)

selected_reactions = [model.reactions[i] for i in model_tinit.flatten().tolist()]
ctx_model = model.copy()

selected_ids = [r.id for r in selected_reactions]
to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]

ctx_model.remove_reactions(to_remove, remove_orphans=True) # TODO: altenative to removing?

print(f"The number of reactions in the base model is: {len(model.reactions)} "
      f"& the number of reactions in the extracted model is: {len(ctx_model.reactions)}")

# export as .xml. # TODO:db integration?
cobra.io.write_sbml_model(ctx_model, f"{mod_dir}/gimme_context_specific_model_{base}.xml")
print('model exported')

solution = ctx_model.optimize()
print('The objective value is: ', solution.objective_value)

# TODO: output the selcted reaction indices as csv (like in gimme)
#get the model summary
print(f"Summary of extracted model: {ctx_model.summary()}")