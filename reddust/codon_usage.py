"""
Codon usage analysis for candidate genes.

Given a coding DNA sequence, compute a codon-usage adaptation score against a
reference codon-frequency table and produce a synonymous recoding that raises
that score while preserving the encoded protein. This is a lightweight
frequency-based heuristic, not a species-calibrated Codon Adaptation Index.
"""

from __future__ import annotations

# Standard genetic code
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reference relative synonymous codon frequencies (fraction within each amino
# acid group). A representative B. subtilis-leaning table; documented so the
# score is fully reproducible and the reference is explicit.
_CODON_FREQ = {
    "A": {"GCT": 0.35, "GCC": 0.40, "GCA": 0.15, "GCG": 0.10},
    "R": {"CGT": 0.25, "CGC": 0.30, "CGA": 0.10, "CGG": 0.15, "AGA": 0.10, "AGG": 0.10},
    "N": {"AAT": 0.45, "AAC": 0.55},
    "D": {"GAT": 0.40, "GAC": 0.60},
    "C": {"TGT": 0.40, "TGC": 0.60},
    "Q": {"CAA": 0.35, "CAG": 0.65},
    "E": {"GAA": 0.45, "GAG": 0.55},
    "G": {"GGT": 0.30, "GGC": 0.35, "GGA": 0.20, "GGG": 0.15},
    "H": {"CAT": 0.45, "CAC": 0.55},
    "I": {"ATT": 0.40, "ATC": 0.50, "ATA": 0.10},
    "L": {"TTA": 0.08, "TTG": 0.12, "CTT": 0.15, "CTC": 0.20, "CTA": 0.10, "CTG": 0.35},
    "K": {"AAA": 0.40, "AAG": 0.60},
    "M": {"ATG": 1.00},
    "F": {"TTT": 0.45, "TTC": 0.55},
    "P": {"CCT": 0.25, "CCC": 0.30, "CCA": 0.25, "CCG": 0.20},
    "S": {"TCT": 0.15, "TCC": 0.20, "TCA": 0.15, "TCG": 0.15, "AGT": 0.15, "AGC": 0.20},
    "T": {"ACT": 0.25, "ACC": 0.40, "ACA": 0.20, "ACG": 0.15},
    "W": {"TGG": 1.00},
    "Y": {"TAT": 0.45, "TAC": 0.55},
    "V": {"GTT": 0.25, "GTC": 0.30, "GTA": 0.15, "GTG": 0.30},
    "*": {"TAA": 0.50, "TAG": 0.20, "TGA": 0.30},
}


def translate(dna: str) -> str:
    dna = dna.upper().replace("U", "T")
    return "".join(_CODON_TABLE.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def adaptation_score(dna: str) -> float:
    """Mean reference frequency of the codons used (0-1; higher = better
    adapted to the reference table)."""
    dna = dna.upper().replace("U", "T")
    vals = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = _CODON_TABLE.get(codon)
        if aa and aa in _CODON_FREQ and codon in _CODON_FREQ[aa]:
            vals.append(_CODON_FREQ[aa][codon])
    return sum(vals) / len(vals) if vals else 0.0


def optimize(dna: str) -> str:
    """Recode each codon to the highest-frequency synonym for its amino acid,
    preserving the encoded protein."""
    if not dna or len(dna) % 3 != 0:
        raise ValueError("Sequence must be non-empty and a multiple of 3.")
    dna = dna.upper().replace("U", "T")
    out = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        aa = _CODON_TABLE.get(codon)
        if aa and aa in _CODON_FREQ:
            best = max(_CODON_FREQ[aa], key=_CODON_FREQ[aa].get)
            out.append(best)
        else:
            out.append(codon)
    return "".join(out)


def gc_content(dna: str) -> float:
    dna = dna.upper()
    return 100.0 * (dna.count("G") + dna.count("C")) / len(dna) if dna else 0.0
