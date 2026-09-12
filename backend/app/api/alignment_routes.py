"""
Route API pour l'alignement pairwise de deux séquences.
"""

from fastapi import APIRouter

from app.core.alignment import align_sequences
from app.schemas.alignment_schemas import (
    AlignmentRequest,
    AlignmentResponse,
    alignment_result_to_response,
)

router = APIRouter(prefix="/api/v1", tags=["alignment"])


@router.post("/alignment", response_model=AlignmentResponse)
def create_alignment(request: AlignmentRequest) -> AlignmentResponse:
    result = align_sequences(request.seq1, request.seq2, request.mode)

    return alignment_result_to_response(result)