import re
import os
import argparse
from math import log
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import CustomSelectionIntegrationStrategy
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties

from .base import BaseAlgorithm


class FastcoreAlgorithm(BaseAlgorithm):
    """FastCORE algorithm for context-specific metabolic model extraction"""

    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to expression data csv")
        p.add_argument("--model", help="Path to SBML model")
        return p.parse_args()

    @property
    def algorithm_name(self):
        return "fastcore"

    def run(self):
        # Load model and expression
        self.model = self.load_model()
        self.expression_data = self.load_expression()

        # Get gene pattern for GPR parsing
        pattern = self.model_cfg.get(
            "alt_transcript_pattern", "__COBAMPGPRDOT__[0-9]{1}"
        )
        patt = re.compile(pattern)
        replace_alt_transcripts = lambda x: patt.sub("", x)

        # Prepare expression data for Troppo
        expr_transposed = self.expression_data.T

        omics_container = TabularReader(
            path_or_df=expr_transposed,
            nomenclature="gene",
            omics_type="transcriptomics",
        ).to_containers()

        single_sample = omics_container[0]

        model_wrapper = ReconstructionWrapper(
            model=self.model,
            ttg_ratio=9999,
            gpr_gene_parse_function=replace_alt_transcripts,
        )

        data_map = single_sample.get_integrated_data_map(
            model_reader=model_wrapper.model_reader, and_func=min, or_func=sum
        )

        # Expression threshold for core reaction selection
        fastcore_cfg = self.cfg.get("fastcore", {})
        threshold = fastcore_cfg.get("threshold", 50 * log(2))
        solver = fastcore_cfg.get("solver", "CPLEX")
        protected_reactions = fastcore_cfg.get("protected_reactions", [])

        # Select core reactions — those above threshold or protected
        def integration_fx(reaction_map_scores):
            return [
                [
                    k
                    for k, v in reaction_map_scores.get_scores().items()
                    if (v is not None and v > threshold) or k in protected_reactions
                ]
            ]

        threshold_integration = CustomSelectionIntegrationStrategy(
            group_functions=[integration_fx]
        )
        threshold_scores = threshold_integration.integrate(data_map=data_map)

        # Get indices of core reactions
        ordered_ids = {r: i for i, r in enumerate(model_wrapper.model_reader.r_ids)}
        core_idx = [[ordered_ids[k] for k in l] for l in threshold_scores]

        # Run FastCORE
        properties = FastcoreProperties(core=core_idx, solver=solver)

        fastcore = FASTcore(
            S=model_wrapper.S,
            lb=model_wrapper.lb,
            ub=model_wrapper.ub,
            properties=properties,
        )

        fastcore_result = fastcore.run()

        # Build context-specific model
        base = os.path.splitext(os.path.basename(self.expr_path))[0]

        selected_reactions = [self.model.reactions[i] for i in fastcore_result]
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