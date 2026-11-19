"""Trait vocabulary: the biotic stress entries recognise the diseases, pests and pathogens by their common and Latin names."""
import pytest

from tools.gene_prioritization.backend import traits

CASES = {
    "Brown spot disease": ("brown_spot", "Bipolaris oryzae"),
    "Cochliobolus miyabeanus": ("brown_spot", "brown spot"),
    "Rice blast": ("blast", "Magnaporthe"),
    "Bacterial leaf blight": ("bacterial_blight", "Xanthomonas"),
    "Sheath blight": ("sheath_blight", "Rhizoctonia"),
    "False smut": ("false_smut", "Ustilaginoidea"),
    "Karnal bunt": ("smut_bunt", "Tilletia"),
    "Northern corn leaf blight": ("leaf_spot_blotch", "Exserohilum turcicum"),
    "Fusarium head blight": ("fusarium", "Fusarium"),
    "Bakanae": ("fusarium", "Fusarium fujikuroi"),
    "Aflatoxin": ("mycotoxin_aspergillus", "Aspergillus flavus"),
    "Wheat stripe rust": ("rust_wheat", "stripe rust"),
    "Powdery mildew": ("mildew", "Blumeria"),
    "Bacterial wilt": ("wilt", "Ralstonia"),
    "Rice tungro virus": ("virus", "tungro"),
    "Brown planthopper resistance": ("planthopper_leafhopper", "Nilaparvata lugens"),
    "Yellow stem borer": ("stem_borer", "Scirpophaga incertulas"),
    "Fall armyworm": ("lepidoptera_defoliator", "Spodoptera frugiperda"),
    "Russian wheat aphid": ("sucking_pest", "Diuraphis noxia"),
    "Asian rice gall midge": ("midge_fly", "Orseolia oryzae"),
    "Rice weevil": ("storage_pest", "Sitophilus oryzae"),
    "Western corn rootworm": ("root_pest", "Diabrotica"),
    "Root-knot nematode": ("nematode", "Meloidogyne"),
    "Striga resistance": ("parasitic_weed", "Striga hermonthica"),
    "Phytoplasma": ("phytoplasma", "phytoplasma"),
    "Biotic stress tolerance": ("disease", "biotic stress"),
}


@pytest.mark.parametrize("trait", list(CASES))
def test_biotic_traits_are_recognised(trait):
    key, term = CASES[trait]
    r = traits.expand(trait)
    assert key in r["concepts"], r["concepts"]
    assert term.lower() in [d.lower() for d in r["direct"]]
    assert r["note"] == ""


def test_bacterial_leaf_blight_is_not_a_fungal_leaf_blight():
    assert "leaf_spot_blotch" not in traits.expand("Bacterial leaf blight")["concepts"]


@pytest.mark.parametrize("trait", ["Drought tolerance", "Drought susceptibility index", "Leaf wilting score", "Cysteine content", "Heat tolerance", "Grain yield"])
def test_abiotic_and_agronomic_traits_do_not_pick_up_biotic_entries(trait):
    biotic = {"disease", "insect", "wilt", "nematode", "fusarium", "mildew", "virus", "blast"}
    assert not biotic & set(traits.expand(trait)["concepts"])


def test_concept_keys_are_unique():
    keys = [c["key"] for c in traits.CONCEPTS]
    assert len(keys) == len(set(keys))
