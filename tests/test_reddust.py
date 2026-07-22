"""Tests for RedDustReclaimer. Metabolic tests are skipped if the iYO844
model file or cobra is unavailable, so the codon tests always run."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from reddust import codon_usage as cu

MODEL_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "models", "iYO844.xml",
)


# --- codon usage -----------------------------------------------------------

def test_translate_basic():
    assert cu.translate("ATGAAATTTGGGTAG") == "MKFG*"


def test_optimize_preserves_protein():
    seq = "ATGAAATTTGGGTGGCTGACC"
    assert cu.translate(seq) == cu.translate(cu.optimize(seq))


def test_optimize_rejects_bad_length():
    with pytest.raises(ValueError):
        cu.optimize("ATGC")
    with pytest.raises(ValueError):
        cu.optimize("")


def test_score_in_range():
    s = cu.adaptation_score("ATGAAATTTGGGTAG")
    assert 0.0 <= s <= 1.0


def test_optimize_does_not_reduce_score():
    seq = "TTATTATTACTACTA"  # low-frequency leucine codons
    assert cu.adaptation_score(cu.optimize(seq)) >= cu.adaptation_score(seq)


# --- metabolic model -------------------------------------------------------

cobra = pytest.importorskip("cobra")
pytestmark = pytest.mark.skipif(
    not os.path.exists(MODEL_PATH), reason="iYO844.xml not present"
)


def _load():
    from reddust import perchlorate_fba as pf
    return pf, pf.load_model(MODEL_PATH)


def test_model_dimensions():
    _, m = _load()
    assert len(m.genes) == 844


def test_pathway_is_balanced():
    pf, m = _load()
    integrated = pf.add_perchlorate_pathway(m)  # asserts balance internally
    assert "PCR" in integrated.reactions
    assert "CLD" in integrated.reactions
    assert not integrated.reactions.PCR.check_mass_balance()
    assert not integrated.reactions.CLD.check_mass_balance()


def test_neutral_when_oxygen_replete():
    pf, m = _load()
    integrated = pf.add_perchlorate_pathway(m)
    a = pf.growth(integrated, o2_uptake_cap=6, perchlorate_uptake_cap=0)
    b = pf.growth(integrated, o2_uptake_cap=6, perchlorate_uptake_cap=20)
    assert abs(a["growth"] - b["growth"]) < 1e-6


def test_benefit_under_oxygen_limitation():
    pf, m = _load()
    integrated = pf.add_perchlorate_pathway(m)
    a = pf.growth(integrated, o2_uptake_cap=1.3, perchlorate_uptake_cap=0)
    b = pf.growth(integrated, o2_uptake_cap=1.3, perchlorate_uptake_cap=20)
    assert b["growth"] > a["growth"]
    assert b["cld_flux"] > 0
