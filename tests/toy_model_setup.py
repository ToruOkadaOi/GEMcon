"""Build a larger toy COBRA model for testing algorithms."""

import cobra
import pandas as pd
import os


def build_toy_model():
    """Create a larger COBRA model with 18 reactions, 9 metabolites, 6 genes."""
    model = cobra.Model("toy_model")

    m_glc = cobra.Metabolite("glc_c", compartment="c")
    m_g6p = cobra.Metabolite("g6p_c", compartment="c")
    m_f6p = cobra.Metabolite("f6p_c", compartment="c")
    m_gap = cobra.Metabolite("gap_c", compartment="c")
    m_pyr = cobra.Metabolite("pyr_c", compartment="c")
    m_lac = cobra.Metabolite("lac_c", compartment="c")
    m_atp = cobra.Metabolite("atp_c", compartment="c")
    m_adp = cobra.Metabolite("adp_c", compartment="c")
    m_biom = cobra.Metabolite("biom_c", compartment="c")

    metabolites = [m_glc, m_g6p, m_f6p, m_gap, m_pyr, m_lac, m_atp, m_adp, m_biom]
    model.add_metabolites(metabolites)

    r_glc_up = cobra.Reaction("r_glc_up", name="Glucose uptake")
    r_glc_up.add_metabolites({m_glc: 1})
    r_glc_up.lower_bound = 0
    r_glc_up.upper_bound = 10

    r_hex = cobra.Reaction("r_hex", name="Hexokinase")
    r_hex.add_metabolites({m_glc: -1, m_g6p: 1, m_adp: -1, m_atp: -1})
    r_hex.lower_bound = 0
    r_hex.upper_bound = 1000
    r_hex.gene_reaction_rule = "G1"

    r_pgi = cobra.Reaction("r_pgi", name="PGI")
    r_pgi.add_metabolites({m_g6p: -1, m_f6p: 1})
    r_pgi.lower_bound = 0
    r_pgi.upper_bound = 1000
    r_pgi.gene_reaction_rule = "G1"

    r_pfk = cobra.Reaction("r_pfk", name="PFK")
    r_pfk.add_metabolites({m_f6p: -1, m_atp: -1, m_gap: 1, m_adp: -1})
    r_pfk.lower_bound = 0
    r_pfk.upper_bound = 1000
    r_pfk.gene_reaction_rule = "G1 and G2"

    r_gapd = cobra.Reaction("r_gapd", name="GAPDH")
    r_gapd.add_metabolites({m_gap: -1, m_adp: -1, m_pyr: 1, m_atp: 1})
    r_gapd.lower_bound = 0
    r_gapd.upper_bound = 1000
    r_gapd.gene_reaction_rule = "G2"

    r_pyk = cobra.Reaction("r_pyk", name="Pyruvate kinase")
    r_pyk.add_metabolites({m_gap: -1, m_adp: -1, m_pyr: 1})
    r_pyk.lower_bound = 0
    r_pyk.upper_bound = 1000
    r_pyk.gene_reaction_rule = "G2"

    r_ldh = cobra.Reaction("r_ldh", name="LDH")
    r_ldh.add_metabolites({m_pyr: -1, m_lac: 1})
    r_ldh.lower_bound = 0
    r_ldh.upper_bound = 1000
    r_ldh.gene_reaction_rule = "G3"

    r_atpase = cobra.Reaction("r_atpase", name="ATP maintenance")
    r_atpase.add_metabolites({m_atp: -1, m_adp: 1})
    r_atpase.lower_bound = 0
    r_atpase.upper_bound = 1000
    r_atpase.gene_reaction_rule = "G4"

    r_gln = cobra.Reaction("r_gln", name="Glutaminase")
    r_gln.add_metabolites({m_g6p: -1})
    r_gln.lower_bound = 0
    r_gln.upper_bound = 100
    r_gln.gene_reaction_rule = "G5"

    r_oxrt = cobra.Reaction("r_oxrt", name="Other reaction 1")
    r_oxrt.add_metabolites({m_f6p: -1, m_biom: 1})
    r_oxrt.lower_bound = 0
    r_oxrt.upper_bound = 100
    r_oxrt.gene_reaction_rule = "G6"

    r_oxrt2 = cobra.Reaction("r_oxrt2", name="Other reaction 2")
    r_oxrt2.add_metabolites({m_pyr: -1, m_biom: 0.5})
    r_oxrt2.lower_bound = 0
    r_oxrt2.upper_bound = 50
    r_oxrt2.gene_reaction_rule = "G5 or G6"

    r_biom_out = cobra.Reaction("r_biom_out", name="Biomass export")
    r_biom_out.add_metabolites({m_biom: -1})
    r_biom_out.lower_bound = 0
    r_biom_out.upper_bound = 1000

    r_lac_out = cobra.Reaction("r_lac_out", name="Lactate export")
    r_lac_out.add_metabolites({m_lac: -1})
    r_lac_out.lower_bound = 0
    r_lac_out.upper_bound = 1000

    r_glc_out = cobra.Reaction("r_glc_out", name="Glucose export")
    r_glc_out.add_metabolites({m_glc: -1})
    r_glc_out.lower_bound = -100
    r_glc_out.upper_bound = 0

    r_g6pase = cobra.Reaction("r_g6pase", name="G6Pase")
    r_g6pase.add_metabolites({m_g6p: -1})
    r_g6pase.lower_bound = 0
    r_g6pase.upper_bound = 10
    r_g6pase.gene_reaction_rule = "G5 and G6"

    r_psur = cobra.Reaction("r_psur", name="Peripheral reaction")
    r_psur.add_metabolites({m_f6p: -1})
    r_psur.lower_bound = 0
    r_psur.upper_bound = 20
    r_psur.gene_reaction_rule = "G4 and G5"

    reactions = [
        r_glc_up,
        r_hex,
        r_pgi,
        r_pfk,
        r_gapd,
        r_pyk,
        r_ldh,
        r_atpase,
        r_gln,
        r_oxrt,
        r_oxrt2,
        r_biom_out,
        r_lac_out,
        r_glc_out,
        r_g6pase,
        r_psur,
    ]
    model.add_reactions(reactions)

    model.objective = "r_biom_out"

    return model


def build_toy_expression():
    """Create expression CSV with 6 genes at varying expression levels."""
    expr_data = pd.DataFrame(
        {
            "sample_1": [
                100.0,  # G1 - high (core pathway genes)
                80.0,  # G2 - high
                40.0,  # G3 - medium (lactate dehydrogenase)
                5.0,  # G4 - low (maintenance)
                2.0,  # G5 - very low
                1.0,  # G6 - very low / absent
            ]
        },
        index=pd.Index(["G1", "G2", "G3", "G4", "G5", "G6"], name="gene"),
    )
    return expr_data


def main():
    model = build_toy_model()
    model_path = os.path.join(os.path.dirname(__file__), "toy_model.xml")
    cobra.io.write_sbml_model(model, model_path)
    print(f"Saved model to {model_path} with {len(model.reactions)} reactions")

    expr = build_toy_expression()
    expr_path = os.path.join(os.path.dirname(__file__), "toy_expr.csv")
    expr.to_csv(expr_path)
    print(f"Saved expression to {expr_path}")


if __name__ == "__main__":
    main()
