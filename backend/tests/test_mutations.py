"""
Tests unitaires pour app/core/mutations.py.
"""

import pytest

from app.core.alignment import align_sequences
from app.core.mutations import detect_mutations


def test_detect_mutations_identical_sequences():
    """Deux séquences identiques ne doivent produire aucune mutation."""
    alignment = align_sequences("ATGCGT", "ATGCGT", mode="global")
    report = detect_mutations(alignment)

    assert report.total_mutations == 0
    assert report.mutations == []


def test_detect_mutations_single_substitution():
    """
    Une seule différence ponctuelle doit être détectée comme
    substitution, à la bonne position.
    """
    alignment = align_sequences("ATGCGT", "ATGAGT", mode="global")
    report = detect_mutations(alignment)

    assert report.substitutions_count == 1
    assert report.total_mutations == 1

    mutation = report.mutations[0]
    assert mutation.mutation_type == "substitution"
    assert mutation.seq1_char == "C"
    assert mutation.seq2_char == "A"


def test_detect_mutations_insertion_detected():
    """
    Un caractère supplémentaire dans seq2 doit être détecté comme
    insertion, avec seq1_position à None.
    """
    alignment = align_sequences("ATGT", "ATGCT", mode="global")
    report = detect_mutations(alignment)

    insertions = [m for m in report.mutations if m.mutation_type == "insertion"]
    assert len(insertions) == 1
    assert insertions[0].seq1_position is None
    assert insertions[0].seq2_position is not None


def test_detect_mutations_deletion_detected():
    """
    Un caractère manquant dans seq2 (présent dans seq1) doit être
    détecté comme délétion, avec seq2_position à None.
    """
    alignment = align_sequences("ATGCT", "ATGT", mode="global")
    report = detect_mutations(alignment)

    deletions = [m for m in report.mutations if m.mutation_type == "deletion"]
    assert len(deletions) == 1
    assert deletions[0].seq2_position is None
    assert deletions[0].seq1_position is not None


def test_detect_mutations_mixed_case():
    """
    Séquences combinant plusieurs types de mutations : vérifie que
    les compteurs individuels sont cohérents avec total_mutations.
    """
    alignment = align_sequences("ATGCGTAA", "ATGAGTAAT", mode="global")
    report = detect_mutations(alignment)

    assert report.total_mutations == (
        report.substitutions_count + report.insertions_count + report.deletions_count
    )
    assert report.total_mutations > 0


def test_detect_mutations_rejects_local_alignment():
    """Un AlignmentResult en mode local doit être refusé explicitement."""
    alignment = align_sequences("TTTACCGTAAA", "ACCGT", mode="local")

    with pytest.raises(ValueError):
        detect_mutations(alignment)


def test_detect_mutations_works_with_protein():
    """La logique de détection doit fonctionner identiquement pour les protéines."""
    alignment = align_sequences("MKTAYIAK", "MKTAYVAK", mode="global")
    report = detect_mutations(alignment)

    assert report.seq1_type == "protein"
    assert report.seq2_type == "protein"
    assert report.substitutions_count == 1


def test_detect_mutations_position_consistency():
    """
    Vérifie que seq1_position/seq2_position correspondent bien aux
    indices réels dans les séquences d'origine (sans les gaps),
    et non à la position dans l'alignement.
    """
    alignment = align_sequences("ATGT", "ATGCT", mode="global")
    report = detect_mutations(alignment)

    insertion = next(m for m in report.mutations if m.mutation_type == "insertion")
    # L'insertion se trouve après "ATG" (3 caractères) dans seq2 -> index 3.
    assert insertion.seq2_position == 3


def test_detect_mutations_ambiguous_n_counts_as_substitution():
    """
    Documente le comportement actuel (limitation connue) : un 'N'
    face à une base différente est traité comme une substitution
    ordinaire, sans traitement biologique spécial.
    """
    alignment = align_sequences("ATGNGT", "ATGAGT", mode="global")
    report = detect_mutations(alignment)

    substitution = next(
        m for m in report.mutations
        if m.seq1_char == "N" or m.seq2_char == "N"
    )
    assert substitution.mutation_type == "substitution"