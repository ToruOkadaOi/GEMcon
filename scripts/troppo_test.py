import pandas as pd
import numpy as np
import cobra
import re

from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

# parsing rule 
patt = re.compile('__COBAMPGPRDOT__[0-9]{1}') # e.g __COBAMPGPRDOT__2
replace_alt_transcripts = lambda x: patt.sub('', x) #  empty string (prune)

# load model and expression data
model = cobra.io.read_sbml_model('/home/biodata/aman/Human-GEM/model/Human-GEM.xml')
expression_data = pd.read_csv('expression_data.csv', index_col=0)

omics_container = TabularReader(path_or_df=expression_data, nomenclature='symbol', # mine is not entrezid is it?
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

properties = tINITProperties(
    reactions_scores=[v for k, v in scores.items()],
    solver='CPLEX'
) # works

tinit = tINIT(
    S=model_wrapper.S,
    lb=model_wrapper.lb,
    ub=model_wrapper.ub,
    properties=properties
)

model_tinit = tinit.run()

cobra.io.write_sbml_model(model_tinit, "tinit_context_specific_model.xml")
model_tinit.to_json("tinit_context_specific_model.json")
