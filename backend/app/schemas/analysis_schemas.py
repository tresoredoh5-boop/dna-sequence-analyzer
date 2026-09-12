"""
Schémas Pydantic pour l'endpoint /api/v1/analysis.

Ce fichier ne contient aucune logique biologique : uniquement la forme
des requêtes/réponses JSON et la conversion vers/depuis AnalysisResult.
"""

from typing import Dict, Optional

from pydantic import BaseModel, Field

from app.core.analysis import AnalysisResult


class AnalysisRequest(BaseModel):
    """Requête d'analyse d'une séquence unique."""
    # Pas de min_length=1 ici volontairement : detect_sequence_type("")
    # est un cas déjà géré explicitement par le core (retourne "unknown",
    # ne lève pas d'erreur). Imposer un rejet Pydantic sur ce cas
    # contredirait ce comportement métier déjà validé par 38 tests.
    sequence: str = Field(..., description="Séquence brute (ADN, ARN ou protéine)")


class AnalysisResponse(BaseModel):
    """Réponse contenant les statistiques calculées pour une séquence."""
    sequence: str
    sequence_type: str
    length: int
    composition: Dict[str, int]
    gc_content: Optional[float]
    molecular_weight: Optional[float]


def analysis_result_to_response(result: AnalysisResult) -> AnalysisResponse:
    """Convertit un AnalysisResult (dataclass core) en AnalysisResponse (API)."""
    return AnalysisResponse(
        sequence=result.sequence,
        sequence_type=result.sequence_type,
        length=result.length,
        composition=result.composition,
        gc_content=result.gc_content,
        molecular_weight=result.molecular_weight,
    )