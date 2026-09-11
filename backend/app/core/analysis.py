"""
Module d'orchestration : assemble le chargement, la détection et les
statistiques d'une séquence en un résultat unique et typé.

Ce module ne fait aucun calcul lui-même. Il ne fait que :
1. Convertir un SeqRecord en chaîne brute (une seule fois, ici uniquement).
2. Déléguer les calculs à stats.compute_statistics().
3. Assembler le tout dans un objet AnalysisResult.
"""

from dataclasses import dataclass
from typing import Dict, Optional

from Bio.SeqRecord import SeqRecord

from app.core.stats import compute_statistics


@dataclass(frozen=True)
class AnalysisResult:
    """
    Résultat figé d'une analyse complète de séquence.

    Toute la suite du projet (alignement, mutations, API) doit consommer
    cet objet plutôt que de reconvertir un SeqRecord soi-même.
    """
    record_id: str
    description: str
    sequence: str
    sequence_type: str
    length: int
    composition: Dict[str, int]
    gc_content: Optional[float]
    molecular_weight: Optional[float]


def analyze_sequence(record: SeqRecord) -> AnalysisResult:
    """
    Analyse un SeqRecord et retourne un AnalysisResult complet.

    Args:
        record: séquence chargée par sequence.load_fasta().

    Returns:
        AnalysisResult regroupant métadonnées et statistiques.
    """
    # Seule conversion SeqRecord -> str de tout le pipeline.
    sequence_str = str(record.seq)

    stats = compute_statistics(sequence_str)

    return AnalysisResult(
        record_id=record.id,
        description=record.description,
        sequence=sequence_str,
        sequence_type=stats["type"],
        length=stats["length"],
        composition=stats["composition"],
        gc_content=stats["gc_content"],
        molecular_weight=stats["molecular_weight"],
    )