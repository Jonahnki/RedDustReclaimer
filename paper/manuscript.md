# Oxygen released by chlorite dismutase can support growth of *Bacillus subtilis* under oxygen limitation: a genome-scale modeling study of a Mars-relevant perchlorate detoxification pathway

**John Adedeji**
Osun State University, Osogbo, Nigeria
Correspondence: johnnaadedeji2018@gmail.com

---

## Abstract

Perchlorate salts are abundant in Martian regolith and are toxic to most
organisms, making their microbial removal a recurring theme in proposals for
in-situ resource use and eventual soil remediation. The bacterial perchlorate
respiration pathway is unusual in that its second enzyme, chlorite dismutase,
releases molecular oxygen. We tested, using constraint-based metabolic
modeling, whether adding this pathway to the model bacterium *Bacillus
subtilis* would change predicted growth, and under what conditions. We
integrated correctly balanced perchlorate reductase (PcrAB) and chlorite
dismutase (Cld) reactions into the published genome-scale model iYO844 and
computed growth by flux balance analysis across a range of external oxygen
availabilities. When oxygen is plentiful the pathway is metabolically neutral
and carries no flux. As oxygen becomes limiting, the oxygen liberated by
chlorite dismutase supplements respiration: predicted growth rises by 25–200%
relative to cells without perchlorate, and below the model's aerobic
feasibility threshold the pathway rescues growth that is otherwise impossible.
At an external oxygen uptake cap of 1.3 mmol gDW⁻¹ h⁻¹, chlorite dismutase
supplies an additional 0.32 mmol gDW⁻¹ h⁻¹ of oxygen internally, a 25%
increment over the external supply. These results identify oxygen
self-supply — not detoxification per se — as the mechanism by which a
perchlorate pathway could benefit an engineered host under the low-oxygen
conditions characteristic of Mars, and they give a specific, testable
hypothesis for experimental work. All code, data, and the model are provided
and reproduce the reported numbers.

**Keywords:** perchlorate, chlorite dismutase, flux balance analysis, iYO844,
*Bacillus subtilis*, Mars, oxygen limitation

---

## 1. Introduction

The Martian surface contains perchlorate (ClO₄⁻) at concentrations of roughly
0.5–1% by weight, first measured directly by the Phoenix lander and since
corroborated at other sites [1]. Perchlorate is chemically stable, water
soluble, and toxic to most terrestrial life at elevated concentrations, and it
is a practical obstacle to both human presence and any attempt to establish
biological activity in Martian soil. A microorganism able to reduce perchlorate
would therefore be of interest both as a soil-conditioning agent and as a
model system for studying life-support-relevant metabolism under Martian
constraints.

A well-characterized bacterial pathway already performs this reduction. In
dissimilatory perchlorate-reducing bacteria such as *Dechloromonas* and
*Azospira* species, perchlorate reductase (PcrAB) reduces perchlorate to
chlorite, and chlorite dismutase (Cld) then converts chlorite to chloride and
molecular oxygen [2,3]. The oxygen-generating step is biochemically remarkable:
chlorite dismutase is one of very few enzymes that produce O₂, and in
perchlorate reducers the liberated oxygen is used by the same cell as a
terminal electron acceptor, effectively allowing anaerobic access to an
oxidant [3].

This oxygen byproduct is the feature we focus on here. Mars has an atmosphere
that is only about 0.13% oxygen and roughly 0.6% of Earth's surface pressure,
so any heterotroph engineered for Martian conditions would face severe oxygen
limitation. If a perchlorate pathway can furnish oxygen internally, its value
to an engineered host may lie less in detoxification and more in respiratory
support. Whether this is metabolically significant, and under what conditions,
is a quantitative question that can be examined before any strain is built.

Genome-scale metabolic models make such questions tractable. For the
industrially and genetically well-developed chassis *Bacillus subtilis*, the
manually curated model iYO844 reconstructs 844 genes and the associated
reaction network and has been validated against growth and gene-essentiality
data [4]. By adding a small number of reactions to such a model and applying
flux balance analysis (FBA), one can predict whether and when an engineered
pathway changes growth, without committing to laboratory work whose outcome is
uncertain.

Here we integrate the PcrAB/Cld pathway into iYO844 and ask a deliberately
narrow question: does the pathway change predicted growth of *B. subtilis*, and
if so, is the effect attributable to the oxygen released by chlorite
dismutase? We report that the pathway is neutral under oxygen-replete
conditions but becomes beneficial, and eventually essential for growth, as
oxygen is restricted — a result that reframes the pathway's Mars relevance
around oxygen self-supply and yields a concrete experimental prediction.

## 2. Related work

**Microbial perchlorate reduction.** The biochemistry and ecology of
dissimilatory perchlorate reduction are well established [2]. Perchlorate
reductase is a periplasmic molybdoenzyme; chlorite dismutase is a heme enzyme
whose O₂-forming mechanism has been studied structurally and kinetically
[3]. These organisms couple perchlorate reduction to growth, using the oxygen
produced from chlorite as an electron acceptor. This body of work establishes
that the pathway is real, genetically transferable in principle, and
oxygen-generating — the premises on which the present study rests.

**Perchlorate on Mars and its biological implications.** Direct detection of
perchlorate in Martian soil [1] prompted a substantial literature on its
implications for habitability, for the interpretation of past life-detection
experiments, and for in-situ resource utilization. Proposals to use
perchlorate-reducing microbes for oxygen generation or soil treatment have been
advanced conceptually; what has been missing is quantitative, model-based
assessment of whether a common laboratory chassis would actually benefit from
carrying the pathway under Mars-like constraints.

**Constraint-based modeling of engineered pathways.** Flux balance analysis is
a standard method for predicting the growth consequences of adding or removing
reactions in genome-scale models [5]. Its application to *B. subtilis* rests on
the iYO844 reconstruction [4]. FBA does not require kinetic parameters and is
therefore well suited to an early, hypothesis-generating assessment of a
proposed engineering strategy, with the understanding that its predictions are
statements about metabolic feasibility rather than about rates observed in
vivo.

The present work sits at the intersection of these three threads: it takes an
established, oxygen-generating pathway, places it in a validated chassis model,
and evaluates its effect specifically under the oxygen limitation that
distinguishes the Martian context from ordinary aerobic culture.

## 3. Methods

### 3.1 Base model

We used iYO844, the genome-scale metabolic reconstruction of *Bacillus
subtilis* 168 [4], in SBML level 3 (fbc version 2) format. On loading with
COBRApy the model contained 844 genes, 1250 reactions, and 990 metabolites,
with biomass reaction `BIOMASS_BS_10` as the objective. Unless otherwise
stated, exchange bounds were left at their file defaults (glucose uptake
limited to 1.7 mmol gDW⁻¹ h⁻¹; a fixed ATP maintenance requirement of
9 mmol gDW⁻¹ h⁻¹).

### 3.2 Pathway integration

Two reactions were added, both balanced for all elements and charge:

- **Perchlorate reductase (PCR):**
  ClO₄⁻ + 2 FADH₂ → ClO₂⁻ + 2 H₂O + 2 FAD
- **Chlorite dismutase (CLD):**
  ClO₂⁻ → Cl⁻ + O₂

together with intracellular transport of perchlorate in, chloride out, and
exchange reactions for extracellular perchlorate and chloride. FADH₂/FAD, H₂O
and O₂ were the model's existing cytosolic metabolites, so the O₂ produced by
CLD enters the host's oxygen pool directly. Elemental and charge balance of
PCR and CLD was verified programmatically at construction; an unbalanced
reaction raises an error rather than being admitted to the model.

### 3.3 Growth analysis

Growth was computed by flux balance analysis, maximizing the biomass reaction,
using the GLPK linear-programming solver through COBRApy. To probe the effect
of oxygen limitation we varied the external oxygen uptake capacity
(`EX_o2_e` lower bound) from unlimited down to the point of infeasibility,
each time comparing two conditions: perchlorate unavailable (perchlorate
exchange closed) and perchlorate available (uptake permitted up to
20 mmol gDW⁻¹ h⁻¹). Solutions were used only when the solver reported an
optimal status; infeasible conditions are reported as such rather than by
their (meaningless) flux values. Flux variability analysis was used to confirm
whether the pathway reactions were merely permitted or actually required at a
given growth demand.

### 3.4 Codon usage analysis

As a supporting utility for eventual heterologous expression, coding sequences
of candidate genes were scored by a codon-frequency heuristic against a
documented reference table and recoded to the highest-frequency synonym for
each residue, preserving the encoded protein (verified by back-translation).
This is a frequency-based adaptation score, not a species-calibrated Codon
Adaptation Index, and is reported as such.

### 3.5 Availability

All code, the iYO844 model file, the input sequences, and the analysis scripts
are provided in the accompanying repository. The two scripts regenerate the
result tables below and print their key values; the unit tests check model
dimensions, reaction balance, and the qualitative neutral-vs-beneficial
behavior.

## 4. Results

### 4.1 The pathway is mass balanced and, when oxygen is plentiful, inert

Integration added the PCR and CLD reactions to iYO844; both passed elemental
and charge-balance checks. Under oxygen-replete conditions (external oxygen
uptake at or above the model's unconstrained optimum), pathway-integrated
*B. subtilis* grew at 0.118 h⁻¹, identical to the unintegrated model, and both
PCR and CLD carried zero flux whether or not perchlorate was supplied
(Table 1, top rows). In other words, when oxygen is freely available the cell
has no metabolic incentive to run the pathway, and FBA correctly leaves it
idle. This also confirms that the pathway introduces no spurious growth benefit
of its own — an important negative control, since an incorrectly coupled or
unbalanced pathway can inflate predicted growth artifactually.

### 4.2 Under oxygen limitation the pathway supplements respiration

As the external oxygen supply was tightened, the two conditions diverged
(Table 1). With perchlorate available, PCR and CLD began to carry flux, and
predicted growth exceeded that of cells without perchlorate by an increasing
margin: +32% at an external-O₂ cap of 4, +25% at 2, +53% at 1.5, and +202% at
1.3 mmol gDW⁻¹ h⁻¹. Below an external-O₂ cap of about 1.24 mmol gDW⁻¹ h⁻¹ the
model without perchlorate was infeasible — no growth was possible — whereas
with perchlorate available the cell still grew (0.0019 h⁻¹ at a cap of 1.2),
i.e. the pathway rescued growth outright.

**Table 1.** Predicted specific growth rate (h⁻¹) of pathway-integrated
*B. subtilis* as external oxygen supply is reduced, with and without
perchlorate available. PCR and CLD flux are shown for the perchlorate-available
condition (mmol gDW⁻¹ h⁻¹).

| External O₂ cap | Growth, no ClO₄⁻ | Growth, ClO₄⁻ available | PCR = CLD flux | Change |
|---|---|---|---|---|
| unlimited | 0.118 | 0.118 | 0.00 | none |
| 6 | 0.118 | 0.118 | 0.00 | none |
| 4 | 0.080 | 0.106 | 0.97 | +32% |
| 2 | 0.025 | 0.032 | 0.49 | +25% |
| 1.5 | 0.0085 | 0.013 | 0.37 | +53% |
| 1.3 | 0.0019 | 0.0056 | 0.32 | +202% |
| 1.25 | 0.0002 | 0.0038 | 0.31 | +1879% |
| 1.2 | infeasible | 0.0019 | 0.30 | growth rescued |

### 4.3 The benefit is oxygen self-supply by chlorite dismutase

The mechanism is explicit in the flux solution. At an external-O₂ cap of
1.3 mmol gDW⁻¹ h⁻¹ with perchlorate available, the cell drew its full external
allowance of 1.30 mmol gDW⁻¹ h⁻¹ of oxygen and, in addition, generated
0.32 mmol gDW⁻¹ h⁻¹ of oxygen internally through chlorite dismutase — a total
oxygen supply of 1.62 mmol gDW⁻¹ h⁻¹, 25% above what the environment alone
provided. The equal PCR and CLD fluxes across all conditions (Table 1) confirm
that flux runs through the complete pathway to the oxygen-releasing step, not
through perchlorate reduction alone. Flux variability analysis under oxygen
limitation showed nonzero minimum flux for both reactions at high fractions of
the growth optimum, meaning the pathway is not merely permitted but required to
achieve maximal growth once oxygen is scarce — the opposite of its idle status
under oxygen-replete conditions.

### 4.4 Codon usage of candidate genes

For the six candidate coding sequences examined, the frequency-based
adaptation score improved modestly on synonymous recoding — from 0.44–0.47 to
0.46–0.48, gains of 1–4% — while preserving each encoded protein exactly
(Supplementary data, `results/codon_usage.csv`). These are small effects, as
expected from a single-reference frequency heuristic, and are reported only to
document the candidate sequences and the recoding utility; they are not a claim
of optimized expression.

## 5. Discussion

The central finding is mechanistic and, in retrospect, intuitive: a pathway
whose terminal step produces oxygen is worthless to a cell that already has
oxygen, and valuable to one that does not. Flux balance analysis of
*B. subtilis* carrying the perchlorate pathway makes this quantitative. The
pathway is inert under aerobic conditions and progressively beneficial as
oxygen becomes limiting, up to rescuing growth that is otherwise infeasible.
Because Mars imposes exactly the oxygen-poor regime in which the pathway helps,
the result reframes the appeal of engineering perchlorate reduction into a
Martian chassis: the near-term metabolic payoff is respiratory oxygen supply,
with perchlorate detoxification as a secondary consequence.

Several caveats bound this claim. FBA predicts metabolic feasibility and
optimal flux allocation under stated bounds; it does not model enzyme kinetics,
perchlorate or chlorite toxicity, regulatory constraints, or the metabolic cost
of expressing the introduced enzymes. The oxygen benefit we report is an
upper-bound, stoichiometric statement — the cell *can* use CLD-derived oxygen
to grow more — not a guarantee that a real strain would. The magnitude of the
benefit also depends on the assumed perchlorate and glucose uptake capacities,
which we varied only coarsely. And the model represents a single reference
condition; a full treatment would couple oxygen limitation to the temperature,
pressure, and nutrient constraints of an actual Martian scenario. These are the
natural next steps, and the provided code is structured to support them.

The value of the study is therefore twofold. Substantively, it isolates a
specific, testable prediction — that perchlorate availability should improve
growth of a PcrAB/Cld-expressing *B. subtilis* under microaerophilic but not
aerobic culture, in proportion to the severity of oxygen limitation — which can
be examined directly in a controlled-atmosphere growth experiment.
Methodologically, it illustrates that the benefit of an oxygen-generating
pathway is invisible under the aerobic conditions in which such models are
usually evaluated, and appears only when the relevant environmental stress is
imposed; the same caution applies to other engineered pathways whose value is
conditional on the environment.

## 6. Conclusion

Adding a correctly balanced perchlorate reductase and chlorite dismutase
pathway to the *Bacillus subtilis* model iYO844 produces no growth change under
oxygen-replete conditions but a substantial, oxygen-limitation-dependent
benefit — from a 25–200% growth increase to outright growth rescue — driven by
the oxygen that chlorite dismutase releases. This positions oxygen self-supply,
rather than detoxification alone, as the mechanistic basis for engineering
perchlorate metabolism into a Mars-relevant chassis, and provides a concrete
hypothesis for experimental testing. The complete, reproducible analysis is
provided.

## Data and code availability

All code, the iYO844 model, input sequences, analysis scripts, and unit tests
are available in the accompanying repository. The scripts regenerate every
number reported here.

## References

[1] Hecht MH, Kounaves SP, Quinn RC, et al. Detection of perchlorate and the
soluble chemistry of Martian soil at the Phoenix lander site. *Science*.
2009;325(5936):64–67.

[2] Coates JD, Achenbach LA. Microbial perchlorate reduction: rocket-fuelled
metabolism. *Nature Reviews Microbiology*. 2004;2(7):569–580.

[3] Lee AQ, Streit BR, Zdilla MJ, Abu-Omar MM, DuBois JL. Mechanism of and
exquisite selectivity for O–O bond formation by the heme-dependent chlorite
dismutase. *Proceedings of the National Academy of Sciences*.
2008;105(41):15654–15659.

[4] Oh YK, Palsson BO, Park SM, Schilling CH, Mahadevan R. Genome-scale
reconstruction of metabolic network in *Bacillus subtilis* based on
high-throughput phenotyping and gene essentiality data. *Journal of Biological
Chemistry*. 2007;282(39):28791–28799.

[5] Orth JD, Thiele I, Palsson BO. What is flux balance analysis?
*Nature Biotechnology*. 2010;28(3):245–248.
