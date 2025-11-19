import pandas as pd
import numpy as np
import cobra
import re
import os
import argparse
import yaml
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

p = argparse.ArgumentParser()
p.add_argument("--expr", help="Path to expression CSV") # or just e instead?? # go cobra-cli subcommands
p.add_argument("--model", help="Path to the SBML model file")
args = p.parse_args()


patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
replace_alt_transcripts = lambda x: patt.sub('', x)

# for recon3d
# patt = re.compile(r'_AT[0-9]+$')
# replace_alt_transcripts = lambda x: patt.sub('', x)

if args.model:
    model_path = args.model.strip()
else:
    if os.path.exists("config.yaml"):
        with open("config.yaml") as f:
            cfg = yaml.safe_load(f)
        model_path = cfg.get("model")
    else:
        model_path = None

if not model_path:
    model_path = input("\nProvide full path to the SBML model (.xml): ").strip()

# model_path = input("\nProvide full path to the SBML model (.xml): ").strip()
if not os.path.exists(model_path):
    raise FileNotFoundError(model_path)

# expr_path = input("Provide full path to the expression CSV file: ").strip()
# if not os.path.exists(expr_path):
#     raise FileNotFoundError(expr_path)

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

expression_data_transposed = expression_data.T

print("\nTransposed shape:", expression_data_transposed.shape)
print(expression_data_transposed.head())

omics_container = TabularReader(path_or_df=expression_data_transposed, 
                                nomenclature='gene',
                                omics_type='transcriptomics').to_containers()

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

# search for biomass in reaction names
# for r in model.reactions:
#     if 'biomass' in r.name.lower() or 'growth' in r.name.lower():
#         print(r.id, r.name)


# print(model_wrapper.model_reader.r_ids)
# Get the index of the biomass reaction in the model. This will be used as objective for the GIMME algorithm. # TODO: important
idx_objective = model_wrapper.model_reader.r_ids.index('MAR02388') # chaneg each time maybe
# idx_objective = model_wrapper.model_reader.r_ids.index('BIOMASS_reaction')
# Create the properties for the GIMME algorithm.
properties = GIMMEProperties(exp_vector=[v for k, v in scores.items()], obj_frac=0.1, objectives=[{idx_objective: 1}],
                             preprocess=True, flux_threshold=0.25, solver='CPLEX',
                             reaction_ids= model_wrapper.model_reader.r_ids, metabolite_ids=model_wrapper.model_reader.m_ids)

# Run the GIMME algorithm.
gimme = GIMME(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

gimme_run = gimme.run()
# print(type(gimme_run))
# print(gimme_run)
pd.DataFrame(gimme_run).to_csv(f'gimme_output_{base}.tsv', sep='\t')
print(pd.DataFrame(gimme_run).head())

#---
# model.reactions
# len(model.reactions)
selected_reactions = [model.reactions[i] for i in gimme_run]
print(pd.DataFrame(selected_reactions).head())
#---
ctx_model = model.copy()

selected_ids = [r.id for r in selected_reactions]
to_remove = [r for r in ctx_model.reactions if r.id not in selected_ids]

ctx_model.remove_reactions(to_remove, remove_orphans=True)

#print(ctx_model)
print(len(model.reactions), len(ctx_model.reactions))
# export
cobra.io.write_sbml_model(ctx_model, f"gimme_context_specific_model_{base}.xml")
solution = ctx_model.optimize()
print(solution.objective_value)
# print(ctx_model.summary())
# fluxes = gimme.sol.__dict__['_Solution__value_map']
# print(fluxes).head()