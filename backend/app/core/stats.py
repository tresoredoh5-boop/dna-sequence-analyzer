"""
Module de calcul de statistiques descriptives sur une séquence biologique.

Ce module suppose que la séquence a déjà été chargée et typée
(voir app/core/sequence.py). Il ne fait aucune détection lui-même,
sauf dans compute_statistics(), qui orchestre l'ensemble.
"""

from collections import Counter
from typing import Dict, Optional

from Bio.SeqUtils import gc_fraction, molecular_weight

from app.core.sequence import detect_sequence_type


def sequence_length(sequence: str) -> int:
    """
    Retourne la longueur de la séquence.

    Args:
        sequence: séquence brute sous forme de chaîne.

    Returns:
        Nombre de caractères dans la séquence (0 si vide).
    """
    return len(sequence)


def base_composition(sequence: str) -> Dict[str, int]:
    """
    Compte le nombre d'occurrences de chaque caractère dans la séquence.

    Fonctionne pour n'importe quel type de séquence (ADN, ARN, protéine) :
    on ne fait aucune hypothèse sur l'alphabet, on compte simplement
    ce qui est présent.

    Args:
        sequence: séquence brute sous forme de chaîne.

    Returns:
        Dictionnaire {caractère: nombre d'occurrences}, ex. {"A": 5, "T": 3}.
    """
    seq_upper = sequence.upper()
    return dict(Counter(seq_upper))


def gc_content(sequence: str, seq_type: str) -> float:
    """
    Calcule le pourcentage de G+C dans une séquence ADN ou ARN.

    Args:
        sequence: séquence brute sous forme de chaîne.
        seq_type: type détecté de la séquence ("DNA", "RNA", "protein", "unknown").

    Returns:
        Pourcentage de G+C, entre 0.0 et 100.0.

    Raises:
        ValueError: si seq_type n'est pas "DNA" ou "RNA" (le %GC n'a pas
                    de sens pour une protéine ou un type inconnu).
    """
    if seq_type not in ("DNA", "RNA"):
        raise ValueError(
            f"Le calcul du %GC n'est pas applicable pour le type '{seq_type}'."
        )

    # gc_fraction retourne une valeur entre 0 et 1 -> on convertit en pourcentage.
    return round(gc_fraction(sequence) * 100, 2)


def molecular_weight_of(sequence: str, seq_type: str) -> Optional[float]:
    """
    Calcule le poids moléculaire approximatif de la séquence, en utilisant
    la fonction native de Biopython (Bio.SeqUtils.molecular_weight).

    Args:
        sequence: séquence brute sous forme de chaîne.
        seq_type: type détecté de la séquence ("DNA", "RNA", "protein", "unknown").

    Returns:
        Poids moléculaire en Daltons, ou None si le type est "unknown"
        (impossible de calculer un poids fiable sans savoir de quoi il s'agit).

    Raises:
        ValueError: si la séquence contient des caractères ambigus (ex. "N"),
                    car Biopython exige un alphabet non ambigu pour ce calcul.
    """
    if seq_type == "unknown":
        return None

    try:
        return round(molecular_weight(sequence, seq_type=seq_type), 2)
    except ValueError as error:
        # On relance une erreur plus explicite pour l'utilisateur du module.
        raise ValueError(
            f"Impossible de calculer le poids moléculaire : {error}"
        ) from error


def compute_statistics(sequence: str) -> Dict:
    """
    Fonction principale : calcule l'ensemble des statistiques disponibles
    pour une séquence donnée, en détectant son type automatiquement.

    C'est cette fonction qui sera appelée depuis l'API plus tard.

    Args:
        sequence: séquence brute sous forme de chaîne.

    Returns:
        Dictionnaire contenant :
            - "length": longueur de la séquence
            - "type": type détecté ("DNA", "RNA", "protein", "unknown")
            - "composition": composition en caractères
            - "gc_content": pourcentage de GC (uniquement si DNA/RNA, sinon None)
            - "molecular_weight": poids moléculaire (None si type "unknown"
              ou si le calcul échoue à cause de caractères ambigus)
    """
    seq_type = detect_sequence_type(sequence)

    result = {
        "length": sequence_length(sequence),
        "type": seq_type,
        "composition": base_composition(sequence),
        "gc_content": None,
        "molecular_weight": None,
    }

    if seq_type in ("DNA", "RNA"):
        result["gc_content"] = gc_content(sequence, seq_type)

    try:
        result["molecular_weight"] = molecular_weight_of(sequence, seq_type)
    except ValueError:
        # Séquence avec caractères ambigus : on garde molecular_weight à None
        # plutôt que de faire planter toute l'analyse.
        result["molecular_weight"] = None

    return result