"""
Route API pour l'analyse d'une séquence unique.

Aucune logique biologique ici : la route construit un SeqRecord minimal
à partir de la chaîne reçue (Option 1, validée), puis délègue tout le
calcul à analyze_sequence() (app/core/analysis.py, non modifié).
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastapi import APIRouter

from app.core.analysis import analyze_sequence
from app.schemas.analysis_schemas import (
    AnalysisRequest,
    AnalysisResponse,
    analysis_result_to_response,
)

router = APIRouter(prefix="/api/v1", tags=["analysis"])


@router.post("/analysis", response_model=AnalysisResponse)
def create_analysis(request: AnalysisRequest) -> AnalysisResponse:
    # SeqRecord minimal : seul le contenu de la séquence importe ici,
    # aucun en-tête FASTA n'est fourni par une requête API brute.
    record = SeqRecord(Seq(request.sequence), id="api_request", description="")

    result = analyze_sequence(record)

    return analysis_result_to_response(result)