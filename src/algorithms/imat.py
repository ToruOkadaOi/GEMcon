"""
Status: Testing
"""
import re
import os
import argparse
import numpy as np
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.imat import IMAT, IMATProperties

from .base import BaseAlgorithm


class ImatAlgorithm(BaseAlgorithm):
    """iMAT algorithm for context-specific metabolic model extraction"""
    
    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to expression data csv")
        p.add_argument("--model", help="Path to SBML model")
        return p.parse_args()
    
    @property
    def algorithm_name(self):
        return "imat"
    
    def run(self):
        # Load model and expression
        self.model = self.load_model()
        self.expression_data = self.load_expression()
        
        # Get gene pattern for GPR parsing
        pattern = self.model_cfg.get('alt_transcript_pattern', '__COBAMPGPRDOT__[0-9]{1}')
        patt = re.compile(pattern)
        replace_alt_transcripts = lambda x: patt.sub('', x)
        
        # Prepare expression data for Troppo
        expr_transposed = self.expression_data.T
        
        omics_container = TabularReader(
            path_or_df=expr_transposed,
            nomenclature='gene',
            omics_type='transcriptomics'
        ).to_containers()
        
        single_sample = omics_container[0]
        
        model_wrapper = ReconstructionWrapper(
            model=self.model,
            ttg_ratio=9999,
            gpr_gene_parse_function=replace_alt_transcripts
        )
        
        data_map = single_sample.get_integrated_data_map(
            model_reader=model_wrapper.model_reader,
            and_func=min,
            or_func=sum
        )
        
        def score_apply(reaction_map_scores):
            return {k: 0 if v is None else v for k, v in reaction_map_scores.items()}
        
        continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
        scores = continuous_integration.integrate(data_map=data_map)
        
        # Build expression vector in reaction order
        reaction_ids = model_wrapper.model_reader.r_ids
        exp_vector = np.array([scores[rid] for rid in reaction_ids])
        
        # Get iMAT-specific parameters
        imat_cfg = self.cfg.get('imat', {})
        exp_thresholds = tuple(imat_cfg.get('exp_thresholds', [25, 75]))
        solver = imat_cfg.get('solver', 'CPLEX')
        
        # Run iMAT
        properties = IMATProperties(
            exp_vector=exp_vector,
            exp_thresholds=exp_thresholds
        )
        
        imat = IMAT(
            S=model_wrapper.S,
            lb=model_wrapper.lb,
            ub=model_wrapper.ub,
            properties=properties
        )
        
        imat_result = imat.run()
        
        # Build context-specific model
        base = os.path.splitext(os.path.basename(self.expr_path))[0]
        
        selected_reactions = [self.model.reactions[i] for i in imat_result.flatten().tolist()]
        selected_ids = [r.id for r in selected_reactions]
        
        ctx_model = self.create_context_model(self.model, selected_ids)
        
        print(f"Base model reactions: {len(self.model.reactions)}")
        print(f"Context model reactions: {len(ctx_model.reactions)}")
        
        # Export model
        output_path = self.export_model(ctx_model, base)
        print(f"Saved to: {output_path}")
        
        solution = ctx_model.optimize()
        print(f"Objective value: {solution.objective_value}")
        
        return ctx_model