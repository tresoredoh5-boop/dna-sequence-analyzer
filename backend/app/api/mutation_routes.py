"""
Route API pour la détection de mutations entre deux séquences.

Option A (validée) : le client fournit seq1/seq2/mode, jamais un
AlignmentResult. La route effectue l'alignement en interne, puis la
détection de mutations, dans une seule requête HTTP.
"""

from fastapi import APIRouter

from app.core.alignment import align_sequences
from app.core.mutations import detect_mutations
from app.schemas.mutation_schemas import (
    MutationRequest,
    MutationResponse,
    mutation_report_to_response,
)

router = APIRouter(prefix="/api/v1", tags=["mutations"])


@router.post("/mutations", response_model=MutationResponse)
def create_mutation_report(request: MutationRequest) -> MutationResponse:
    alignment_result = align_sequences(request.seq1, request.seq2, request.mode)

    report = detect_mutations(alignment_result)

    return mutation_report_to_response(report)