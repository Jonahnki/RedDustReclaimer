#!/usr/bin/env python3
"""Robustness checks for the perchlorate O2-rescue result:

(A) parsimonious FBA (pFBA) cross-check alongside plain FBA at the O2 caps
    used in the main result, confirming both agree on growth rate; and
(B) an NGAM (non-growth maintenance) sensitivity sweep at the headline
    rescue scenario, confirming the qualitative rescue effect holds across
    a +/-20% range of the maintenance demand.

Writes results/sensitivity_checks.csv (tidy long format) and prints a
plain-English summary.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reddust import perchlorate_fba as pf

# Parameters copied from scripts/02_perchlorate_o2_analysis.py (the script that
# produced the paper's headline table) so nothing is re-derived or guessed.
O2_CAPS = [None, 6, 4, 2, 1.5, 1.3, 1.25, 1.2]   # matches script 02
PERCHLORATE_CAP = 20                              # matches script 02
RESCUE_O2_CAP = 1.2                               # infeasible w/o perchlorate
NGAM_REACTION_ID = "ATPM"                         # confirmed by inspection


def _g(d, key):
    return None if d.get("growth") is None else round(d.get(key, 0.0), 5)


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    model = pf.load_model(os.path.join(here, "models", "iYO844.xml"))
    integrated = pf.add_perchlorate_pathway(model)
    ngam = integrated.reactions.get_by_id(NGAM_REACTION_ID)
    print(f"Model: {len(integrated.genes)} genes; "
          f"NGAM reaction {NGAM_REACTION_ID} fixed at {ngam.lower_bound} "
          f"(bounds {ngam.bounds})\n")

    rows = []  # tidy long format

    # ---- (A) pFBA cross-check at every main-result O2 cap ------------------
    print("=" * 74)
    print("(A) FBA vs pFBA growth-rate cross-check")
    print("=" * 74)
    print(f"{'O2 cap':>9} | {'ClO4 cap':>8} | {'FBA growth':>11} | "
          f"{'pFBA growth':>11} | {'match':>6}")
    all_match = True
    for cap in O2_CAPS:
        for pcap in (0, PERCHLORATE_CAP):
            fba = pf.growth(integrated, o2_uptake_cap=cap, perchlorate_uptake_cap=pcap)
            pf_ = pf.growth_pfba(integrated, o2_uptake_cap=cap, perchlorate_uptake_cap=pcap)
            gf, gp = fba["growth"], pf_["growth"]
            if gf is None or gp is None:
                match = (gf is None and gp is None)
            else:
                match = abs(gf - gp) < 1e-6
            all_match = all_match and match
            label = "unlimited" if cap is None else cap
            print(f"{str(label):>9} | {pcap:>8} | "
                  f"{str('infeasible' if gf is None else round(gf,5)):>11} | "
                  f"{str('infeasible' if gp is None else round(gp,5)):>11} | "
                  f"{str(match):>6}")
            for method, d in (("fba", fba), ("pfba", pf_)):
                rows.append({
                    "analysis": "pfba_crosscheck",
                    "o2_cap": "unlimited" if cap is None else cap,
                    "perchlorate_cap": pcap,
                    "ngam_delta": "",
                    "ngam_flux": ngam.lower_bound,
                    "method": method,
                    "growth": "infeasible" if d["growth"] is None else round(d["growth"], 6),
                    "pcr_flux": _g(d, "pcr_flux"),
                    "cld_flux": _g(d, "cld_flux"),
                    "status": d["status"],
                })

    # ---- (B) NGAM sensitivity at the headline rescue scenario -------------
    print("\n" + "=" * 74)
    print(f"(B) NGAM +/-20% sensitivity at rescue scenario "
          f"(O2 cap={RESCUE_O2_CAP})")
    print("=" * 74)
    print(f"{'delta':>6} | {'NGAM flux':>9} | {'ClO4':>5} | {'growth':>10} | "
          f"{'PCR':>7} | {'CLD':>7} | {'status':>10}")
    # Run the sweep both without and with perchlorate so the rescue contrast
    # (infeasible w/o perchlorate, feasible with) is visible at each NGAM level.
    sweeps = {}
    for pcap in (0, PERCHLORATE_CAP):
        sweeps[pcap] = pf.ngam_sensitivity(
            integrated, NGAM_REACTION_ID,
            o2_uptake_cap=RESCUE_O2_CAP, perchlorate_uptake_cap=pcap)
    # print interleaved by delta for readability
    for i, delta in enumerate([r["delta"] for r in sweeps[PERCHLORATE_CAP]]):
        for pcap in (0, PERCHLORATE_CAP):
            entry = sweeps[pcap][i]
            d = entry["result"]
            print(f"{entry['delta']:>+6.2f} | {entry['ngam_flux_used']:>9.3f} | "
                  f"{pcap:>5} | "
                  f"{str('infeasible' if d['growth'] is None else round(d['growth'],5)):>10} | "
                  f"{str(_g(d,'pcr_flux')):>7} | {str(_g(d,'cld_flux')):>7} | "
                  f"{d['status']:>10}")
            rows.append({
                "analysis": "ngam_sensitivity",
                "o2_cap": RESCUE_O2_CAP,
                "perchlorate_cap": pcap,
                "ngam_delta": entry["delta"],
                "ngam_flux": round(entry["ngam_flux_used"], 4),
                "method": "fba",
                "growth": "infeasible" if d["growth"] is None else round(d["growth"], 6),
                "pcr_flux": _g(d, "pcr_flux"),
                "cld_flux": _g(d, "cld_flux"),
                "status": d["status"],
            })

    # ---- write CSV --------------------------------------------------------
    os.makedirs(os.path.join(here, "results"), exist_ok=True)
    out_path = os.path.join(here, "results", "sensitivity_checks.csv")
    fields = ["analysis", "o2_cap", "perchlorate_cap", "ngam_delta",
              "ngam_flux", "method", "growth", "pcr_flux", "cld_flux", "status"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # ---- plain-English summary -------------------------------------------
    rescue_held = all(
        sweeps[0][i]["result"]["growth"] is None
        and sweeps[PERCHLORATE_CAP][i]["result"]["growth"] is not None
        for i in range(len(sweeps[0]))
    )
    print("\n" + "=" * 74)
    print("SUMMARY")
    print("=" * 74)
    print(
        f"pFBA vs FBA growth-rate agreement: "
        f"{'ALL MATCH' if all_match else 'MISMATCH DETECTED - INVESTIGATE'} "
        f"(across {len(O2_CAPS)} O2 caps x 2 perchlorate conditions). By "
        f"construction pFBA holds the biomass flux at the FBA optimum while "
        f"minimising total flux, so growth should match; internal flux "
        f"distributions may differ.\n")
    if rescue_held:
        print(
            f"Rescue effect: HELD at all NGAM deltas ({[r['delta'] for r in sweeps[0]]}). "
            f"At O2 cap {RESCUE_O2_CAP}, growth was infeasible without perchlorate "
            f"and feasible with perchlorate available at every maintenance level "
            f"from -20% to +20% of the default NGAM.")
    else:
        print(
            f"Rescue effect: DID NOT hold uniformly across the NGAM range. See "
            f"table above and CSV for the exact delta(s) where the "
            f"infeasible-without / feasible-with contrast breaks down.")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
