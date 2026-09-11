"""
Tests unitaires pour app/core/stats.py.

Ces tests vérifient chaque fonction indépendamment (longueur, composition,
%GC, poids moléculaire), ainsi que le comportement d'orchestration de
compute_statistics() sur différents types de séquences.
"""

import pytest

from app.core.stats import (
    sequence_length,
    base_composition,
    gc_content,
    molecular_weight_of,
    compute_statistics,
)


def test_sequence_length():
    """La longueur retournée doit correspondre au nombre de caractères."""
    assert sequence_length("ATGC") == 4
    assert sequence_length("") == 0


def test_base_composition_counts_correctly():
    """
    La composition doit compter chaque caractère correctement,
    indépendamment de la casse (normalisation en majuscules).
    """
    result = base_composition("AATTGGCC")
    assert result == {"A": 2, "T": 2, "G": 2, "C": 2}

    result_lower = base_composition("aattggcc")
    assert result_lower == result


def test_gc_content_for_dna():
    """
    Une séquence ADN avec un ratio GC connu doit retourner
    exactement le pourcentage attendu.
    """
    # "ATGC" -> 2 GC sur 4 bases -> 50%
    assert gc_content("ATGC", "DNA") == 50.0


def test_gc_content_for_rna():
    """Le %GC doit aussi être calculable pour une séquence ARN."""
    # "AUGC" -> 2 GC sur 4 bases -> 50%
    assert gc_content("AUGC", "RNA") == 50.0


def test_gc_content_rejects_protein_and_unknown():
    """
    Le %GC n'a pas de sens biologique pour une protéine ou un type
    inconnu : la fonction doit lever une ValueError explicite.
    """
    with pytest.raises(ValueError):
        gc_content("MKTAYIAKQRQISFVK", "protein")

    with pytest.raises(ValueError):
        gc_content("XYZ123", "unknown")


def test_molecular_weight_for_valid_dna():
    """
    Le poids moléculaire d'une séquence ADN valide doit être un nombre
    positif, calculé via Bio.SeqUtils.molecular_weight().
    """
    weight = molecular_weight_of("ATGC", "DNA")
    assert weight is not None
    assert weight > 0


def test_molecular_weight_with_ambiguous_n_raises_error():
    """
    Bio.SeqUtils.molecular_weight() n'accepte que des lettres non
    ambiguës : une séquence contenant 'N' doit lever une ValueError
    lorsqu'on appelle molecular_weight_of() directement.
    """
    with pytest.raises(ValueError):
        molecular_weight_of("ATGCN", "DNA")


def test_compute_statistics_full_dna_sequence():
    """
    Pour une séquence ADN valide, compute_statistics() doit retourner
    un dictionnaire complet avec tous les champs correctement calculés.
    """
    result = compute_statistics("ATGC")

    assert result["length"] == 4
    assert result["type"] == "DNA"
    assert result["composition"] == {"A": 1, "T": 1, "G": 1, "C": 1}
    assert result["gc_content"] == 50.0
    assert result["molecular_weight"] is not None


def test_compute_statistics_handles_ambiguous_n_gracefully():
    """
    compute_statistics() doit rester robuste face à une séquence
    contenant 'N' : le %GC reste calculable (N est ignoré par défaut),
    mais le poids moléculaire doit être None plutôt que de faire
    planter toute l'analyse.
    """
    result = compute_statistics("ATGCN")

    assert result["length"] == 5
    assert result["type"] == "DNA"
    assert result["gc_content"] is not None
    assert result["molecular_weight"] is None