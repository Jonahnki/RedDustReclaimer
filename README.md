# RedDustReclaimer

A small, reproducible computational study of whether an engineered perchlorate
detoxification pathway could benefit *Bacillus subtilis* under the
oxygen-limited conditions relevant to Martian regolith.

Perchlorate salts are abundant in Martian soil. The bacterial enzymes that
break them down — perchlorate reductase (PcrAB) and chlorite dismutase (Cld) —
are notable because the second step releases molecular oxygen. This repository
integrates that pathway into the published genome-scale metabolic model of
*B. subtilis* (iYO844) and asks a focused question: does the pathway change
predicted growth, and if so, under what conditions?

The short answer, reproduced by the scripts here: the pathway is metabolically
neutral when oxygen is plentiful, but under oxygen limitation the oxygen
released by chlorite dismutase supplements respiration and improves predicted
growth — from a modest gain in mild limitation to full growth rescue below the
model's aerobic feasibility threshold.

## What's here

```
reddust/
  perchlorate_fba.py   pathway integration + growth analysis for iYO844
  codon_usage.py       codon-usage scoring and synonymous recoding
scripts/
  01_codon_analysis.py            -> results/codon_usage.csv
  02_perchlorate_o2_analysis.py   -> results/perchlorate_o2_benefit.csv
models/iYO844.xml      base genome-scale model (Oh et al., 2007)
data/                  candidate coding sequences and enzyme references
tests/                 unit tests
paper/manuscript.md    the write-up
```

## Reproducing the results

```bash
# environment (conda/micromamba recommended)
micromamba create -y -f environment.yml
micromamba activate reddust
# or: pip install -r requirements.txt

python scripts/01_codon_analysis.py
python scripts/02_perchlorate_o2_analysis.py
pytest -q
```

Both scripts write CSV files to `results/` and print their key numbers. The
pathway reactions are checked for elemental and charge balance at build time,
so an unbalanced edit cannot pass silently.

## Key result

Predicted specific growth rate of pathway-integrated *B. subtilis* as the
external oxygen supply is tightened, with and without perchlorate available:

| External O₂ cap (mmol gDW⁻¹ h⁻¹) | Growth, no perchlorate | Growth, perchlorate available | Change |
|---|---|---|---|
| unlimited | 0.118 | 0.118 | none |
| 4 | 0.080 | 0.106 | +32% |
| 2 | 0.025 | 0.032 | +25% |
| 1.5 | 0.0085 | 0.013 | +53% |
| 1.3 | 0.0019 | 0.0056 | +202% |
| 1.2 | infeasible | 0.0019 | growth rescued |

At an external-O₂ cap of 1.3, chlorite dismutase supplies an additional
0.32 mmol gDW⁻¹ h⁻¹ of oxygen internally — about a 25% increase over the
external supply.

## Scope and limitations

This is a constraint-based (flux balance) modeling study. It predicts what is
metabolically feasible and favorable under stated assumptions; it does not
measure enzyme kinetics, perchlorate toxicity, or in-vivo behavior. The
codon-usage tool is a frequency-based heuristic against a documented reference
table, not a species-calibrated Codon Adaptation Index. Results are a
hypothesis to be tested experimentally, not a validated strain design.

## License

MIT — see `LICENSE`.
