#!/usr/bin/env python3
"""Core analysis: does a perchlorate detoxification pathway benefit growth
under oxygen limitation?

Integrates the PCR/CLD pathway into iYO844 and compares growth with and
without perchlorate available across a range of external-O2 caps. Writes
results/perchlorate_o2_benefit.csv and prints the mechanism at one condition.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reddust import perchlorate_fba as pf


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    model = pf.load_model(os.path.join(here, "models", "iYO844.xml"))
    print(f"Loaded iYO844: {len(model.genes)} genes, "
          f"{len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    integrated = pf.add_perchlorate_pathway(model)
    print("Pathway integrated; PCR and CLD verified mass/charge balanced.")

    # External-O2 caps to sweep (positive = uptake capacity).
    o2_caps = [None, 6, 4, 2, 1.5, 1.3, 1.25, 1.2]
    rows = []
    for cap in o2_caps:
        no_p = pf.growth(integrated, o2_uptake_cap=cap, perchlorate_uptake_cap=0)
        with_p = pf.growth(integrated, o2_uptake_cap=cap, perchlorate_uptake_cap=20)

        def g(d):
            return None if d["growth"] is None else round(d["growth"], 5)

        benefit = "n/a"
        a, b = g(no_p), g(with_p)
        if a is None and b is not None:
            benefit = "rescue"
        elif a and b is not None:
            benefit = f"{100 * (b - a) / a:+.0f}%"
        rows.append({
            "o2_uptake_cap": "unlimited" if cap is None else cap,
            "growth_no_perchlorate": "infeasible" if a is None else a,
            "growth_with_perchlorate": "infeasible" if b is None else b,
            "pcr_flux": None if with_p["growth"] is None else round(with_p.get("pcr_flux", 0), 4),
            "cld_flux": None if with_p["growth"] is None else round(with_p.get("cld_flux", 0), 4),
            "benefit": benefit,
        })

    os.makedirs(os.path.join(here, "results"), exist_ok=True)
    out_path = os.path.join(here, "results", "perchlorate_o2_benefit.csv")
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_path}")
    for r in rows:
        print(f"  O2cap={str(r['o2_uptake_cap']):>9}  no_ClO4={str(r['growth_no_perchlorate']):>10}"
              f"  with_ClO4={str(r['growth_with_perchlorate']):>10}  benefit={r['benefit']}")

    # Mechanism at one microaerophilic condition.
    d = pf.growth(integrated, o2_uptake_cap=1.3, perchlorate_uptake_cap=20)
    ext = -d["o2_exchange"]
    internal = d["cld_flux"]
    print(f"\nMechanism at external-O2 cap = 1.3 mmol/gDW/h:")
    print(f"  growth                 : {d['growth']:.5f} 1/h")
    print(f"  external O2 uptake      : {ext:.4f} mmol/gDW/h")
    print(f"  CLD-generated O2        : {internal:.4f} mmol/gDW/h")
    print(f"  total O2 to cell        : {ext + internal:.4f} mmol/gDW/h "
          f"(+{100 * internal / ext:.0f}% over external alone)")


if __name__ == "__main__":
    main()
