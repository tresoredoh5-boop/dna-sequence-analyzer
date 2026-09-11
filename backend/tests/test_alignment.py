"""
Tests unitaires pour app/core/alignment.py.
"""

import pytest

from app.core.alignment import align_sequences


def test_align_identical_dna_sequences():
    """Deux séquences identiques doivent s'aligner sans aucun gap."""
    result = align_sequences("ATGCGT", "ATGCGT", mode="global")

    assert result.aligned_seq1 == result.aligned_seq2 == "ATGCGT"
    assert result.score == 6.0  # 6 correspondances * match_score (1.0)


def test_align_dna_different_lengths():
    """
    Deux séquences DNA de longueurs différentes doivent produire un
    alignement global avec des gaps, de même longueur des deux côtés.
    """
    result = align_sequences("ATGCGTAA", "ATGCGT", mode="global")

    assert len(result.aligned_seq1) == len(result.aligned_seq2)
    assert "-" in result.aligned_seq1 or "-" in result.aligned_seq2


def test_align_local_with_common_subregion():
    """
    En mode local, seule la sous-région commune doit être retenue :
    ici "ACCGT" est un sous-texte exact de "TTTACCGTAAA".
    """
    result = align_sequences("TTTACCGTAAA", "ACCGT", mode="local")

    assert result.aligned_seq1 == "ACCGT"
    assert result.aligned_seq2 == "ACCGT"
    assert result.score == 5.0


def test_align_sequence_with_ambiguous_n():
    """
    Contrairement à molecular_weight_of(), l'alignement ne rejette pas
    les caractères ambigus comme 'N' : il doit réussir normalement.
    """
    result = align_sequences("ATGCN", "ATGCA", mode="global")

    assert len(result.aligned_seq1) == len(result.aligned_seq2)


def test_align_empty_sequence_raises_error():
    """Une séquence vide doit être refusée explicitement."""
    with pytest.raises(ValueError):
        align_sequences("", "ATGC")


def test_align_invalid_mode_raises_error():
    """Un mode autre que 'global'/'local' doit être refusé."""
    with pytest.raises(ValueError):
        align_sequences("ATGC", "ATGC", mode="semi-global")


def test_align_dna_and_protein_rejected():
    """DNA et protein sont des types incompatibles : doit lever ValueError."""
    with pytest.raises(ValueError):
        align_sequences("ATGC", "MKTAYIAKQRQISFVK")


def test_align_dna_and_rna_rejected():
    """
    DNA et RNA sont volontairement incompatibles dans cette version
    (pas de transcription implémentée) : doit lever ValueError.
    """
    with pytest.raises(ValueError):
        align_sequences("ATGC", "AUGC")


def test_align_two_proteins():
    """Deux séquences protéiques doivent s'aligner normalement."""
    result = align_sequences("MKTAYIAKQRQISFVK", "MKTAYIAKQRQISFVA", mode="global")

    assert result.seq1_type == "protein"
    assert result.seq2_type == "protein"
    assert len(result.aligned_seq1) == len(result.aligned_seq2)


def test_aligned_sequences_always_same_length():
    """
    Invariant fondamental : quel que soit le cas, aligned_seq1 et
    aligned_seq2 doivent toujours avoir la même longueur.
    """
    result = align_sequences("ATGCGTAAGT", "ATGCGT", mode="global")
    assert len(result.aligned_seq1) == len(result.aligned_seq2)