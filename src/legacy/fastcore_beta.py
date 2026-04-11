import pandas as pd
import numpy as np
import cobra
import re
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper
from troppo.omics.integration import CustomSelectionIntegrationStrategy
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.gapfill.fastcc import FastCC, FastCCProperties
from cobra.io import read_sbml_model
import numpy as np

# parsing rule 
patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
replace_alt_transcripts = lambda x: patt.sub('', x)

# load model and expression data
model = cobra.io.read_sbml_model('/home/biodata/aman/Human-GEM/model/Human-GEM.xml')

expression_data = pd.read_csv('/home/biodata/aman/data/data_processed/expression_data_hippocampus-development-human-brain-10XV2_annotated_normalized_annotated_gencode.csv', index_col=0)

expression_data_T = expression_data.T # transpose

omics_container = TabularReader(path_or_df=expression_data_T, nomenclature='gene',
                                omics_type='transcriptomics').to_containers()

single_sample = omics_container[0]

model_wrapper = ReconstructionWrapper(model=model, ttg_ratio=9999,
                                      gpr_gene_parse_function=replace_alt_transcripts)

data_map = single_sample.get_integrated_data_map(model_reader=model_wrapper.model_reader,
                                                 and_func=min, or_func=sum)
from math import log
threshold =  (50 * log(2))

# keywords = [
#     "oxidative", "phosphorylation",
#     "taurine", 
#     "glycerophospholipid",
#     "glutamate", 
#     "glutamine", "aspartate", 
#     "malate", 
#     "glycine", "serine", "threonine",
#     "cholesterol", 
#     "sphingolipid"
# ] # for brain hippocampus dev. ## doi: 10.1021/acschemneuro.5c00006

rid_lists = []

# for kw in keywords:
#     for r in model.reactions:
#         if r.subsystem and kw in r.subsystem.lower():
#             rid_lists.append(r.id)

# print(rid_lists)

if rid_lists:
    protected_reactions = rid_lists # switch to the real rid. 
else:
    protected_reactions = []

def integration_fx(reaction_map_scores):
    return [[k for k, v in reaction_map_scores.get_scores().items() if (v is not None and v > threshold) or k in protected_reactions]]

threshold_integration = CustomSelectionIntegrationStrategy(group_functions=[integration_fx])
threshold_scores = threshold_integration.integrate(data_map=data_map)

print(threshold_scores)
# Get the index of the reaction of the CORE reaction set
ordered_ids = {r:i for i,r in enumerate(model_wrapper.model_reader.r_ids)}
core_idx = [[ordered_ids[k] for k in l] for l in threshold_scores]

# Define the FastCORE properties
properties = FastcoreProperties(core=core_idx, solver='CPLEX')

# instantiate the FastCORE class
fastcore = FASTcore(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

# Run the algorithm
model_fastcore = fastcore.run()