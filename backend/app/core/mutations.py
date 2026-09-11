"""
Module de détection de mutations entre deux séquences alignées.

Ce module ne fait aucun alignement lui-même : il consomme un
AlignmentResult déjà produit par alignment.py, et compare les deux
séquences alignées colonne par colonne pour identifier les mutations.

seq1 est toujours considérée comme la séquence de référence, seq2
comme la séquence comparée. Une insertion/délétion est donc définie
relativement à seq1.
"""

from dataclasses import dataclass
from typing import List, Optional

from app.core.alignment import AlignmentResult


@dataclass(frozen=True)
class Mutation:
    """Une mutation individuelle détectée à une colonne d'alignement."""
    mutation_type: str        # "substitution" | "insertion" | "deletion"
    alignment_position: int   # position dans les séquences alignées (0-based)
    seq1_position: Optional[int]  # position dans seq1 d'origine, None si gap côté seq1
    seq2_position: Optional[int]  # position dans seq2 d'origine, None si gap côté seq2
    seq1_char: str
    seq2_char: str


@dataclass(frozen=True)
class MutationReport:
    """Rapport complet des mutations détectées entre seq1 et seq2."""
    mutations: List[Mutation]
    total_mutations: int
    substitutions_count: int
    insertions_count: int
    deletions_count: int
    seq1_type: str
    seq2_type: str


def detect_mutations(alignment_result: AlignmentResult) -> MutationReport:
    """
    Détecte les mutations (substitutions, insertions, délétions) entre
    les deux séquences d'un AlignmentResult en mode global.

    Args:
        alignment_result: résultat d'alignement produit par align_sequences().

    Returns:
        MutationReport listant toutes les mutations détectées.

    Raises:
        ValueError: si l'alignement fourni n'est pas en mode "global".
    """
    if alignment_result.mode != "global":
        raise ValueError(
            "La détection de mutations nécessite un alignement en mode 'global'."
        )

    aligned_seq1 = alignment_result.aligned_seq1
    aligned_seq2 = alignment_result.aligned_seq2

    mutations: List[Mutation] = []

    # Compteurs de position dans les séquences d'origine (sans les gaps).
    # Ils n'avancent que lorsque le caractère correspondant n'est pas un gap.
    seq1_index = 0
    seq2_index = 0

    for alignment_position, (char1, char2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        is_gap1 = char1 == "-"
        is_gap2 = char2 == "-"

        if is_gap1 and is_gap2:
            # Cas anormal (ne devrait pas se produire avec un aligneur
            # correct) : ignoré, aucun compteur ne bouge.
            continue

        if not is_gap1 and not is_gap2:
            if char1 != char2:
                mutations.append(Mutation(
                    mutation_type="substitution",
                    alignment_position=alignment_position,
                    seq1_position=seq1_index,
                    seq2_position=seq2_index,
                    seq1_char=char1,
                    seq2_char=char2,
                ))
            seq1_index += 1
            seq2_index += 1

        elif is_gap1 and not is_gap2:
            mutations.append(Mutation(
                mutation_type="insertion",
                alignment_position=alignment_position,
                seq1_position=None,
                seq2_position=seq2_index,
                seq1_char=char1,
                seq2_char=char2,
            ))
            seq2_index += 1

        elif not is_gap1 and is_gap2:
            mutations.append(Mutation(
                mutation_type="deletion",
                alignment_position=alignment_position,
                seq1_position=seq1_index,
                seq2_position=None,
                seq1_char=char1,
                seq2_char=char2,
            ))
            seq1_index += 1

    substitutions_count = sum(1 for m in mutations if m.mutation_type == "substitution")
    insertions_count = sum(1 for m in mutations if m.mutation_type == "insertion")
    deletions_count = sum(1 for m in mutations if m.mutation_type == "deletion")

    return MutationReport(
        mutations=mutations,
        total_mutations=len(mutations),
        substitutions_count=substitutions_count,
        insertions_count=insertions_count,
        deletions_count=deletions_count,
        seq1_type=alignment_result.seq1_type,
        seq2_type=alignment_result.seq2_type,
    )