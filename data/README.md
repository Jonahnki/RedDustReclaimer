# Data

## extremophile_sequences.fasta
Six coding DNA sequences from extremophile organisms, used for the codon
usage analysis (`scripts/01_codon_analysis.py`).

## candidate_enzymes/
Protein (amino acid) FASTA files for the enzymes relevant to perchlorate
detoxification and stress tolerance, retrieved from NCBI (accession in each
header):

- `perchlorate_reductase_a_azospira_pmj.fasta` (ALS20366.1)
- `perchlorate_reductase_b_azospira_kj.fasta` (ACB69918.1)
- `chlorite_dismutase_dechloromonas_agitata.fasta` (AAM92878.1)
- `reca_family_radA_ecoli.fasta` (WOY73680.1)
- `uvrabc_excinuclease_atpase_ecoli.fasta` (BAE78060.1)
- `atp_synthase_subunit_b_ecoli.fasta` (PPA54414.1)
- `rusticyanin_acidithiobacillus_ferridurans.fasta` (CAA07038.1)

## models/iYO844.xml
The published *Bacillus subtilis* 168 genome-scale metabolic model iYO844
(Oh et al., 2007), in SBML level 3 (fbc) format. Used as the base model for
the perchlorate pathway analysis.
