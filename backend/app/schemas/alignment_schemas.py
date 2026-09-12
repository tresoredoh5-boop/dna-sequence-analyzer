"""
Schémas Pydantic pour l'endpoint /api/v1/alignment.
"""

from typing import Literal

from pydantic import BaseModel, Field

from app.core.alignment import AlignmentResult


class AlignmentRequest(BaseModel):
    """Requête d'alignement pairwise entre deux séquences."""
    seq1: str = Field(..., min_length=1)
    seq2: str = Field(..., min_length=1)
    mode: Literal["global", "local"] = "global"


class AlignmentResponse(BaseModel):
    """Réponse contenant le résultat complet d'un alignement."""
    seq1: str
    seq2: str
    seq1_type: str
    seq2_type: str
    mode: str
    score: float
    aligned_seq1: str
    aligned_seq2: str


def alignment_result_to_response(result: AlignmentResult) -> AlignmentResponse:
    """Convertit un AlignmentResult (dataclass core) en AlignmentResponse (API)."""
    return AlignmentResponse(
        seq1=result.seq1,
        seq2=result.seq2,
        seq1_type=result.seq1_type,
        seq2_type=result.seq2_type,
        mode=result.mode,
        score=result.score,
        aligned_seq1=result.aligned_seq1,
        aligned_seq2=result.aligned_seq2,
    )