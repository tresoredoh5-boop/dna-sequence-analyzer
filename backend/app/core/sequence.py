"""
Module de base pour le chargement et l'identification de séquences biologiques.

Ce module ne fait volontairement qu'une seule chose : lire des séquences
et déterminer leur nature (ADN, ARN, protéine). Aucune statistique,
aucun alignement ici — ce sera ajouté dans des modules séparés.
"""

from pathlib import Path
from typing import List, Union, IO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Alphabets de référence utilisés pour la détection du type de séquence.
DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

# Longueur en dessous de laquelle une séquence composée uniquement de lettres
# partagées entre ADN et protéine (A, C, G, T, N) est considérée comme
# réellement ambiguë plutôt que classée par défaut en ADN.
# Heuristique simplificatrice, pas une règle biologique universelle.
AMBIGUOUS_LENGTH_THRESHOLD = 4


def load_fasta(source: Union[str, Path, IO]) -> List[SeqRecord]:
    """
    Charge un fichier (ou flux) au format FASTA et retourne la liste
    des séquences qu'il contient.

    Args:
        source: chemin vers un fichier FASTA, ou objet fichier/flux déjà ouvert.

    Returns:
        Liste d'objets SeqRecord.

    Raises:
        ValueError: si aucune séquence n'a pu être lue depuis la source.
    """
    records = list(SeqIO.parse(source, "fasta"))

    if not records:
        raise ValueError("Aucune séquence valide trouvée dans le fichier FASTA fourni.")

    return records


def detect_sequence_type(sequence: str) -> str:
    """
    Détermine le type biologique d'une séquence : ADN, ARN, protéine,
    ambiguë, ou inconnue.

    Args:
        sequence: la séquence sous forme de chaîne de caractères.

    Returns:
        Une chaîne parmi : "DNA", "RNA", "protein", "ambiguous", "unknown".
    """
    seq_upper = sequence.upper()

    if len(seq_upper) == 0:
        return "unknown"

    unique_chars = set(seq_upper)

    protein_exclusive_letters = PROTEIN_ALPHABET - DNA_ALPHABET - RNA_ALPHABET

    if unique_chars & protein_exclusive_letters:
        return "protein"

    if unique_chars.issubset(RNA_ALPHABET) and "U" in unique_chars:
        return "RNA"

    if unique_chars.issubset(DNA_ALPHABET):
        if len(seq_upper) < AMBIGUOUS_LENGTH_THRESHOLD:
            return "ambiguous"
        return "DNA"

    return "unknown"