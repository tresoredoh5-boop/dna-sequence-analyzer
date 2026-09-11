"""
Tests unitaires pour app/core/analysis.py.

Ces tests vérifient que analyze_sequence() assemble correctement
un SeqRecord en un AnalysisResult cohérent, pour les trois types
de séquences principaux (ADN, ARN, protéine), ainsi que la bonne
transmission des métadonnées et la cohérence interne du résultat.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.core.analysis import analyze_sequence


def test_analyze_sequence_dna():
    """
    Une séquence ADN valide et suffisamment longue (>= 4 caractères,
    voir AMBIGUOUS_LENGTH_THRESHOLD) doit être détectée comme "DNA"
    et posséder un gc_content calculé.
    """
    record = SeqRecord(Seq("ATGCGTACGT"), id="seq_dna", description="Test ADN")

    result = analyze_sequence(record)

    assert result.sequence_type == "DNA"
    assert result.gc_content is not None
    assert result.molecular_weight is not None


def test_analyze_sequence_rna():
    """
    Une séquence ARN (contenant U) doit être détectée comme "RNA"
    et posséder un gc_content calculé, comme pour l'ADN.
    """
    record = SeqRecord(Seq("AUGCGUACGU"), id="seq_rna", description="Test ARN")

    result = analyze_sequence(record)

    assert result.sequence_type == "RNA"
    assert result.gc_content is not None
    assert result.molecular_weight is not None


def test_analyze_sequence_protein():
    """
    Une séquence protéique (contenant des lettres exclusivement
    protéiques) doit être détectée comme "protein". Le gc_content
    n'a pas de sens pour une protéine : il doit rester None.
    """
    record = SeqRecord(Seq("MKTAYIAKQRQISFVK"), id="seq_protein", description="Test protéine")

    result = analyze_sequence(record)

    assert result.sequence_type == "protein"
    assert result.gc_content is None
    assert result.molecular_weight is not None


def test_analyze_sequence_preserves_metadata():
    """
    Les métadonnées du SeqRecord (id et description) doivent être
    transférées telles quelles dans l'AnalysisResult, sans altération.
    """
    record = SeqRecord(
        Seq("ATGCGTACGT"),
        id="my_custom_id",
        description="Séquence de démonstration pour les métadonnées",
    )

    result = analyze_sequence(record)

    assert result.record_id == "my_custom_id"
    assert result.description == "Séquence de démonstration pour les métadonnées"


def test_analyze_sequence_internal_consistency():
    """
    Vérifie la cohérence interne du résultat :
    - la chaîne stockée dans 'sequence' doit correspondre exactement
      à la séquence d'origine (conversion str(record.seq) fidèle) ;
    - 'length' doit correspondre à la longueur réelle de cette chaîne ;
    - 'composition' doit être cohérente avec le contenu de la séquence.
    """
    raw_sequence = "ATGCGTACGT"
    record = SeqRecord(Seq(raw_sequence), id="seq_consistency", description="Test cohérence")

    result = analyze_sequence(record)

    assert result.sequence == raw_sequence
    assert result.length == len(raw_sequence)
    assert sum(result.composition.values()) == result.length