"""
Tests unitaires pour app/core/sequence.py, en particulier
detect_sequence_type() et son traitement du chevauchement ADN/protéine.
"""

from app.core.sequence import detect_sequence_type


def test_detect_sequence_type_short_overlap_is_ambiguous():
    """
    'CAT' n'utilise que des lettres partagées entre ADN et protéine
    (C, A, T) et est trop courte pour trancher : elle doit être
    classée 'ambiguous', plus jamais 'DNA' par défaut.
    """
    result = detect_sequence_type("CAT")
    assert result != "DNA"
    assert result == "ambiguous"


def test_detect_sequence_type_dna_stays_dna():
    """Une séquence ADN suffisamment longue (>= seuil) reste classée DNA."""
    assert detect_sequence_type("ATGC") == "DNA"


def test_detect_sequence_type_rna_stays_rna():
    """La présence de U désambiguïse immédiatement en faveur de l'ARN."""
    assert detect_sequence_type("AUGC") == "RNA"


def test_detect_sequence_type_protein_stays_protein():
    """
    Une séquence contenant au moins une lettre exclusivement protéique
    (ici M, K, Y, I, Q, R, S, F, V) reste classée protein sans ambiguïté.
    """
    assert detect_sequence_type("MKTAYIAKQRQISFVK") == "protein"


def test_detect_sequence_type_empty_stays_unknown():
    """Une séquence vide reste classée unknown (non-régression du bug initial)."""
    assert detect_sequence_type("") == "unknown"