"""
Module d'alignement pairwise de séquences biologiques.

Utilise Bio.Align.PairwiseAligner (API moderne de Biopython), pas
Bio.pairwise2, qui est déprécié depuis Biopython 1.80.
"""

from dataclasses import dataclass
from typing import Set

from Bio import Align

from app.core.sequence import detect_sequence_type


VALID_MODES = {"global", "local"}

# Couples de types compatibles pour un alignement biologiquement valide.
# ADN et ARN restent volontairement incompatibles : ce sont deux alphabets
# distincts (T vs U) et aucune transcription ADN<->ARN n'est encore
# implémentée dans le projet.
COMPATIBLE_TYPE_PAIRS: Set[frozenset] = {
    frozenset({"DNA"}),
    frozenset({"RNA"}),
    frozenset({"protein"}),
}


@dataclass(frozen=True)
class AlignmentResult:
    """Résultat figé d'un alignement pairwise entre deux séquences."""
    seq1: str
    seq2: str
    seq1_type: str
    seq2_type: str
    mode: str
    score: float
    aligned_seq1: str
    aligned_seq2: str


def align_sequences(seq1: str, seq2: str, mode: str = "global") -> AlignmentResult:
    """
    Aligne deux séquences biologiques de même nature (DNA-DNA, RNA-RNA,
    ou protein-protein) et retourne le meilleur alignement trouvé.

    Args:
        seq1: première séquence brute.
        seq2: seconde séquence brute.
        mode: "global" (toute la longueur des deux séquences) ou
              "local" (seulement la meilleure sous-région commune).

    Returns:
        AlignmentResult contenant le score et les séquences alignées.

    Raises:
        ValueError: mode invalide, séquence vide, ou types incompatibles.
    """
    if mode not in VALID_MODES:
        raise ValueError(f"Mode d'alignement invalide : '{mode}'. Attendu 'global' ou 'local'.")

    if len(seq1) == 0 or len(seq2) == 0:
        raise ValueError("Impossible d'aligner une séquence vide.")

    seq1_type = detect_sequence_type(seq1)
    seq2_type = detect_sequence_type(seq2)

    if frozenset({seq1_type, seq2_type}) not in COMPATIBLE_TYPE_PAIRS:
        raise ValueError(
            f"Types de séquences incompatibles pour l'alignement : "
            f"'{seq1_type}' et '{seq2_type}'. Seuls DNA-DNA, RNA-RNA et "
            f"protein-protein sont acceptés."
        )

    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]

    # alignment[0] / alignment[1] retournent les séquences alignées,
    # gaps ("-") inclus, toujours de longueur identique entre elles.
    aligned_seq1 = str(best_alignment[0])
    aligned_seq2 = str(best_alignment[1])

    return AlignmentResult(
        seq1=seq1,
        seq2=seq2,
        seq1_type=seq1_type,
        seq2_type=seq2_type,
        mode=mode,
        score=best_alignment.score,
        aligned_seq1=aligned_seq1,
        aligned_seq2=aligned_seq2,
    )