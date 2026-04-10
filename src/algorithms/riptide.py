import os
import argparse
import numpy as np
import pandas as pd
import riptide

from .base import BaseAlgorithm


class RiptideAlgorithm(BaseAlgorithm):
    """RIPTiDe algorithm for context-specific metabolic model extraction"""
    
    def _parse_args(self):
        p = argparse.ArgumentParser()
        p.add_argument("--expr", help="Path to expression data csv")
        p.add_argument("--model", help="Path to SBML model")
        return p.parse_args()
    
    @property
    def algorithm_name(self):
        return "riptide"
    
    def run(self):
        # Load model and expression
        self.model = self.load_model()
        self.expression_data = self.load_expression()
        
        # Build transcriptome dict from the DataFrame
        # RIPTiDe expects {gene_id: [abundance_values]} not a Troppo container
        transcriptome = {}
        for gene_id in self.expression_data.index:
            values = self.expression_data.loc[gene_id].values
            if np.isscalar(values):
                transcriptome[str(gene_id)] = [float(values)]
            else:
                transcriptome[str(gene_id)] = [float(v) for v in values]
        
        # RIPTiDe expects a 'replicates' key
        n_replicates = len(self.expression_data.columns)
        transcriptome["replicates"] = n_replicates
        
        print(f"Transcriptome dict built: {len(transcriptome) - 1} genes, {n_replicates} replicate(s)")
        
        # Get RIPTiDe config
        riptide_cfg = self.cfg.get("riptide", {})
        
        mode = riptide_cfg.get("mode", "contextualize")
        fraction = riptide_cfg.get("fraction", 0.8)
        samples = riptide_cfg.get("samples", 500)
        prune = riptide_cfg.get("prune", True)
        conservative = riptide_cfg.get("conservative", False)
        objective = riptide_cfg.get("objective", True)
        additive = riptide_cfg.get("additive", False)
        set_bounds = riptide_cfg.get("set_bounds", True)
        gpr = riptide_cfg.get("gpr", False)
        threshold = riptide_cfg.get("threshold", 1e-8)
        tasks = riptide_cfg.get("tasks", [])
        exclude = riptide_cfg.get("exclude", [])
        
        # maxfit-specific params
        frac_min = riptide_cfg.get("frac_min", 0.25)
        frac_max = riptide_cfg.get("frac_max", 0.85)
        frac_step = riptide_cfg.get("frac_step", 0.1)
        
        # Run RIPTiDe — contextualize or maxfit mode
        if mode == "maxfit":
            print(f"Running RIPTiDe maxfit (fraction range: {frac_min} - {frac_max}, step: {frac_step})")
            riptide_result = riptide.maxfit(
                model=self.model,
                transcriptome=transcriptome,
                frac_min=frac_min,
                frac_max=frac_max,
                frac_step=frac_step,
                prune=prune,
                samples=samples,
                conservative=conservative,
                objective=objective,
                additive=additive,
                set_bounds=set_bounds,
                tasks=tasks,
                exclude=exclude,
                gpr=gpr,
                threshold=threshold,
            )
        else:
            print(f"Running RIPTiDe contextualize (fraction: {fraction})")
            riptide_result = riptide.contextualize(
                model=self.model,
                transcriptome=transcriptome,
                fraction=fraction,
                samples=samples,
                prune=prune,
                conservative=conservative,
                objective=objective,
                additive=additive,
                set_bounds=set_bounds,
                tasks=tasks,
                exclude=exclude,
                gpr=gpr,
                threshold=threshold,
            )
        
        # RIPTiDe returns the model directly — no need for create_context_model()
        ctx_model = riptide_result.model
        
        print(f"Base model reactions: {len(self.model.reactions)}")
        print(f"Context model reactions: {len(ctx_model.reactions)}")
        
        # Export model
        base = os.path.splitext(os.path.basename(self.expr_path))[0]
        output_path = self.export_model(ctx_model, base)
        print(f"Saved to: {output_path}")
        
        # Save pruned reaction/gene/metabolite info
        if hasattr(riptide_result, "pruned") and isinstance(riptide_result.pruned, dict):
            pruned_rxns = list(riptide_result.pruned.get("reactions", []))
            pruned_genes = list(riptide_result.pruned.get("genes", []))
            pruned_mets = list(riptide_result.pruned.get("metabolites", []))
            
            max_len = max(len(pruned_rxns), len(pruned_genes), len(pruned_mets), 1)
            pruned_df = pd.DataFrame({
                "pruned_reactions": pruned_rxns + [""] * (max_len - len(pruned_rxns)),
                "pruned_genes": pruned_genes + [""] * (max_len - len(pruned_genes)),
                "pruned_metabolites": pruned_mets + [""] * (max_len - len(pruned_mets)),
            })
            pruned_outpath = os.path.join(self.output_dir, f"riptide_pruned_{base}.tsv")
            pruned_df.to_csv(pruned_outpath, sep="\t", index=False)
            print(f"Saved pruned components to: {pruned_outpath}")
        
        # Save concordance info
        if hasattr(riptide_result, "concordance") and isinstance(riptide_result.concordance, dict):
            r_val = riptide_result.concordance.get("r", "N/A")
            p_val = riptide_result.concordance.get("p", "N/A")
            print(f"Concordance: Spearman r={r_val}, p={p_val}")
        
        # Save flux samples if available
        if hasattr(riptide_result, "flux_samples") and isinstance(riptide_result.flux_samples, pd.DataFrame):
            samples_outpath = os.path.join(self.output_dir, f"riptide_flux_samples_{base}.tsv")
            riptide_result.flux_samples.to_csv(samples_outpath, sep="\t")
            print(f"Saved flux samples to: {samples_outpath}")
        
        solution = ctx_model.optimize()
        print(f"Objective value: {solution.objective_value}")
        print(f"Fraction of optimum used: {riptide_result.fraction_of_optimum}")
        
        return ctx_model