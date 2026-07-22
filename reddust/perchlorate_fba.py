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
from cobra.exceptions import Infeasible
from cobra.flux_analysis import pfba


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
        _apply_scenario(m, o2_uptake_cap, perchlorate_uptake_cap)
        sol = m.optimize()
        if sol.status != "optimal":
            return {"status": sol.status, "growth": None}
        return _result_from_solution(m, sol, float(sol.objective_value))


def _apply_scenario(m: cobra.Model, o2_uptake_cap: float | None,
                    perchlorate_uptake_cap: float) -> None:
    """Apply the external-O2 cap and perchlorate supply to an (already
    context-managed) model. Shared by growth() and growth_pfba() so their
    scenario setup can never drift apart."""
    if o2_uptake_cap is not None:
        m.reactions.EX_o2_e.lower_bound = -abs(o2_uptake_cap)
    if "EX_clo4_e" in m.reactions:
        m.reactions.EX_clo4_e.lower_bound = -abs(perchlorate_uptake_cap)


def _result_from_solution(m: cobra.Model, sol, growth_value: float) -> dict:
    """Build the standard result dict from a solution object. `growth_value`
    is passed explicitly because plain FBA reports growth as the objective
    value while pFBA's objective value is the total flux, so the growth rate
    must be read from the biomass flux instead."""
    out = {"status": "optimal", "growth": float(growth_value),
           "o2_exchange": float(sol.fluxes.get("EX_o2_e", 0.0))}
    for rxn in ("PCR", "CLD"):
        if rxn in m.reactions:
            out[rxn.lower() + "_flux"] = float(sol.fluxes[rxn])
    return out


def growth_pfba(model: cobra.Model, o2_uptake_cap: float | None = None,
                perchlorate_uptake_cap: float = 0.0) -> dict:
    """Same scenario setup as growth(), but solved with parsimonious FBA
    (minimum total flux among optimal-growth solutions) as a cross-check
    against the plain FBA result. Returns growth rate and pathway fluxes
    from the pFBA solution, or status='infeasible'.

    Note: pFBA's returned objective_value is the minimised total flux, not
    the growth rate, so growth is read from the biomass reaction flux (which
    equals the FBA optimum by construction of pFBA).
    """
    with model as m:
        _apply_scenario(m, o2_uptake_cap, perchlorate_uptake_cap)
        try:
            sol = pfba(m)
        except Infeasible:
            return {"status": "infeasible", "growth": None}
        biomass_id = next(r.id for r in m.reactions
                          if r.objective_coefficient != 0)
        return _result_from_solution(m, sol, float(sol.fluxes[biomass_id]))


def ngam_sensitivity(model: cobra.Model, ngam_reaction_id: str,
                     o2_uptake_cap: float, perchlorate_uptake_cap: float,
                     deltas: tuple[float, ...] = (-0.20, -0.10, 0.0, 0.10, 0.20)
                     ) -> list[dict]:
    """Re-run growth() at the given O2/perchlorate scenario across a range of
    NGAM flux values, each scaled by (1 + delta) from the model's default
    NGAM bound. Returns a list of dicts, one per delta, each containing:
    delta, ngam_flux_used, and the full growth() output dict for that run.

    Uses the reaction ID confirmed by inspection (parameterised here so it is
    not hardcoded). Assumes the NGAM reaction has equal lower and upper bounds
    (a fixed-flux maintenance demand); this is asserted up front because the
    scaling logic assumes a single fixed value, not a range.
    """
    ngam = model.reactions.get_by_id(ngam_reaction_id)
    if ngam.lower_bound != ngam.upper_bound:
        raise ValueError(
            f"NGAM reaction {ngam_reaction_id} is not fixed-flux "
            f"(bounds {ngam.bounds}); the scaling logic assumes lb == ub."
        )
    original_value = ngam.lower_bound

    results = []
    for delta in deltas:
        scaled = original_value * (1 + delta)
        with model as m:
            r = m.reactions.get_by_id(ngam_reaction_id)
            r.bounds = (scaled, scaled)
            outcome = growth(m, o2_uptake_cap=o2_uptake_cap,
                             perchlorate_uptake_cap=perchlorate_uptake_cap)
        results.append({"delta": delta, "ngam_flux_used": scaled,
                        "result": outcome})
    return results
