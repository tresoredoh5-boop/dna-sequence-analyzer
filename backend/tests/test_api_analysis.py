"""
Tests d'intégration HTTP pour /api/v1/analysis.

Ces tests vérifient uniquement le branchement API <-> core (codes HTTP,
structure JSON), pas la logique métier déjà couverte par test_analysis.py
et test_stats.py.
"""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_analysis_valid_request_returns_200():
    response = client.post("/api/v1/analysis", json={"sequence": "ATGCGTACGT"})

    assert response.status_code == 200
    data = response.json()
    assert data["sequence_type"] == "DNA"
    assert data["length"] == 10
    assert data["gc_content"] is not None


def test_analysis_response_structure():
    response = client.post("/api/v1/analysis", json={"sequence": "MKTAYIAKQRQISFVK"})

    data = response.json()
    expected_keys = {
        "sequence", "sequence_type", "length",
        "composition", "gc_content", "molecular_weight",
    }
    assert expected_keys.issubset(data.keys())
    assert data["gc_content"] is None  # protéine : pas de %GC


def test_analysis_empty_sequence_returns_200_unknown_type():
    """
    Contrairement à /alignment et /mutations, une séquence vide n'est
    PAS une erreur métier pour l'analyse : detect_sequence_type("")
    retourne explicitement "unknown" (comportement déjà validé par
    38 tests core, cf. correction du bug initial). L'API doit donc
    refléter ce comportement fidèlement : 200, pas 400 ni 422.
    """
    response = client.post("/api/v1/analysis", json={"sequence": ""})

    assert response.status_code == 200
    data = response.json()
    assert data["sequence_type"] == "unknown"
    assert data["length"] == 0
    assert data["gc_content"] is None
    assert data["molecular_weight"] is None


def test_analysis_missing_field_returns_422():
    response = client.post("/api/v1/analysis", json={})

    assert response.status_code == 422


def test_analysis_wrong_type_returns_422():
    response = client.post("/api/v1/analysis", json={"sequence": 12345})

    assert response.status_code == 422