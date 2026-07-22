#!/usr/bin/env python3
"""Codon usage analysis of candidate coding sequences.

Reads data/extremophile_sequences.fasta, scores each sequence, produces a
synonymous recoding, and writes results/codon_usage.csv.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reddust import codon_usage as cu


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    return seqs


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    seqs = read_fasta(os.path.join(here, "data", "extremophile_sequences.fasta"))
    os.makedirs(os.path.join(here, "results"), exist_ok=True)
    out_path = os.path.join(here, "results", "codon_usage.csv")

    rows = []
    for name, seq in seqs.items():
        if len(seq) % 3 != 0:
            continue
        base = cu.adaptation_score(seq)
        opt_seq = cu.optimize(seq)
        opt = cu.adaptation_score(opt_seq)
        assert cu.translate(seq) == cu.translate(opt_seq), "protein not preserved"
        rows.append({
            "sequence": name,
            "length_bp": len(seq),
            "gc_percent": round(cu.gc_content(seq), 2),
            "score_baseline": round(base, 3),
            "score_optimized": round(opt, 3),
            "change_percent": round(100 * (opt - base) / base, 1) if base else 0.0,
        })

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {out_path} ({len(rows)} sequences)")
    for r in rows:
        print(f"  {r['sequence'][:45]:45} {r['score_baseline']} -> "
              f"{r['score_optimized']}  ({r['change_percent']:+}%)")


if __name__ == "__main__":
    main()
