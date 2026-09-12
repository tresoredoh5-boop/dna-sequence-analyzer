"""
Schémas Pydantic pour l'endpoint /api/v1/mutations.
"""

from typing import List, Literal, Optional

from pydantic import BaseModel, Field

from app.core.mutations import MutationReport


class MutationRequest(BaseModel):
    """
    Requête de détection de mutations entre deux séquences.

    L'alignement est effectué en interne par la route (Option A) :
    le client ne fournit jamais d'AlignmentResult directement.
    """
    seq1: str = Field(..., min_length=1)
    seq2: str = Field(..., min_length=1)
    mode: Literal["global", "local"] = "global"


class MutationItem(BaseModel):
    """Une mutation individuelle détectée."""
    mutation_type: str
    alignment_position: int
    seq1_position: Optional[int]
    seq2_position: Optional[int]
    seq1_char: str
    seq2_char: str


class MutationResponse(BaseModel):
    """Rapport complet des mutations détectées entre deux séquences."""
    mutations: List[MutationItem]
    total_mutations: int
    substitutions_count: int
    insertions_count: int
    deletions_count: int
    seq1_type: str
    seq2_type: str


def mutation_report_to_response(report: MutationReport) -> MutationResponse:
    """Convertit un MutationReport (dataclass core) en MutationResponse (API)."""
    return MutationResponse(
        mutations=[
            MutationItem(
                mutation_type=m.mutation_type,
                alignment_position=m.alignment_position,
                seq1_position=m.seq1_position,
                seq2_position=m.seq2_position,
                seq1_char=m.seq1_char,
                seq2_char=m.seq2_char,
            )
            for m in report.mutations
        ],
        total_mutations=report.total_mutations,
        substitutions_count=report.substitutions_count,
        insertions_count=report.insertions_count,
        deletions_count=report.deletions_count,
        seq1_type=report.seq1_type,
        seq2_type=report.seq2_type,
    )