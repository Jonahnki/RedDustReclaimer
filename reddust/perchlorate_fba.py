"""
Perchlorate detoxification pathway integration for the Bacillus subtilis
genome-scale model iYO844, and analysis of its effect on growth under
oxygen limitation.

The pathway adds two reactions:

    PCR  (perchlorate reductase, pcrAB):
         ClO4- + 2 FADH2  ->  ClO2- + 2 H2O + 2 FAD
    CLD  (chlorite dismutase, cld):
         ClO2-            ->  Cl- + O2

Both reactions are elementally and charge balanced (verified at build time).
The key biochemical feature is that chlorite dismutase releases molecular
oxygen, which the host can use as a terminal electron acceptor.
"""

from __future__ import annotations

import cobra
from cobra import Metabolite, Reaction


# Correctly balanced pathway reactions. Formulae/charges follow the
# conventions already used inside iYO844 for fadh2_c / fad_c / o2_c / h2o_c.
_NEW_METABOLITES = {
    "clo4_c": ("perchlorate", "ClO4", -1, "c"),
    "clo4_e": ("perchlorate", "ClO4", -1, "e"),
    "clo2_c": ("chlorite", "ClO2", -1, "c"),
    "cl_c": ("chloride", "Cl", -1, "c"),
    "cl_e": ("chloride", "Cl", -1, "e"),
}


def load_model(path: str = "models/iYO844.xml") -> cobra.Model:
    """Load the base iYO844 model from an SBML file."""
    return cobra.io.read_sbml_model(path)


def add_perchlorate_pathway(model: cobra.Model) -> cobra.Model:
    """Return a copy of `model` with the perchlorate detoxification pathway
    integrated. Raises AssertionError if either new reaction is not mass and
    charge balanced, so an unbalanced edit can never pass silently."""
    m = model.copy()

    for mid, (name, formula, charge, comp) in _NEW_METABOLITES.items():
        if mid not in m.metabolites:
            m.add_metabolites(
                [Metabolite(mid, name=name, formula=formula, charge=charge,
                            compartment=comp)]
            )

    pcr = Reaction("PCR", name="Perchlorate reductase (pcrAB)")
    pcr.bounds = (0, 1000)
    pcr.add_metabolites({
        m.metabolites.clo4_c: -1,
        m.metabolites.fadh2_c: -2,
        m.metabolites.clo2_c: 1,
        m.metabolites.h2o_c: 2,
        m.metabolites.fad_c: 2,
    })

    cld = Reaction("CLD", name="Chlorite dismutase (cld)")
    cld.bounds = (0, 1000)
    cld.add_metabolites({
        m.metabolites.clo2_c: -1,
        m.metabolites.cl_c: 1,
        m.metabolites.o2_c: 1,
    })

    clo4t = Reaction("CLO4t", name="Perchlorate uptake")
    clo4t.bounds = (0, 1000)
    clo4t.add_metabolites({m.metabolites.clo4_e: -1, m.metabolites.clo4_c: 1})

    clt = Reaction("CLt", name="Chloride export")
    clt.bounds = (-1000, 1000)
    clt.add_metabolites({m.metabolites.cl_c: -1, m.metabolites.cl_e: 1})

    ex_clo4 = Reaction("EX_clo4_e", name="Perchlorate exchange")
    ex_clo4.add_metabolites({m.metabolites.clo4_e: -1})
    ex_clo4.bounds = (0, 1000)  # set uptake per scenario via lower_bound

    ex_cl = Reaction("EX_cl_e", name="Chloride exchange")
    ex_cl.add_metabolites({m.metabolites.cl_e: -1})
    ex_cl.bounds = (0, 1000)

    m.add_reactions([pcr, cld, clo4t, clt, ex_clo4, ex_cl])

    assert not pcr.check_mass_balance(), f"PCR unbalanced: {pcr.check_mass_balance()}"
    assert not cld.check_mass_balance(), f"CLD unbalanced: {cld.check_mass_balance()}"
    return m


def growth(model: cobra.Model, o2_uptake_cap: float | None = None,
           perchlorate_uptake_cap: float = 0.0) -> dict:
    """Optimize growth under a given external-O2 cap and perchlorate supply.

    Parameters
    ----------
    o2_uptake_cap : maximum external O2 uptake (positive number); None leaves
        the model default. Internally applied as EX_o2_e.lower_bound = -cap.
    perchlorate_uptake_cap : maximum perchlorate uptake (positive number);
        0 means no perchlorate available.

    Returns a dict with growth rate and pathway fluxes, or status='infeasible'.
    """
    with model as m:
        if o2_uptake_cap is not None:
            m.reactions.EX_o2_e.lower_bound = -abs(o2_uptake_cap)
        if "EX_clo4_e" in m.reactions:
            m.reactions.EX_clo4_e.lower_bound = -abs(perchlorate_uptake_cap)
        sol = m.optimize()
        if sol.status != "optimal":
            return {"status": sol.status, "growth": None}
        out = {"status": "optimal", "growth": float(sol.objective_value),
               "o2_exchange": float(sol.fluxes.get("EX_o2_e", 0.0))}
        for rxn in ("PCR", "CLD"):
            if rxn in m.reactions:
                out[rxn.lower() + "_flux"] = float(sol.fluxes[rxn])
        return out
