"""
Tests d'intégration HTTP pour /api/v1/mutations.

Ces tests vérifient le pipeline complet côté API :
requête -> align_sequences() -> detect_mutations() -> réponse.
"""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_mutations_valid_request_returns_200():
    response = client.post(
        "/api/v1/mutations",
        json={"seq1": "ATGCGT", "seq2": "ATGAGT", "mode": "global"},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["total_mutations"] == 1
    assert data["substitutions_count"] == 1


def test_mutations_insertion_and_deletion_detected():
    response = client.post(
        "/api/v1/mutations",
        json={"seq1": "ATGCGTAA", "seq2": "ATGAGTAAT", "mode": "global"},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["total_mutations"] == (
        data["substitutions_count"] + data["insertions_count"] + data["deletions_count"]
    )
    assert data["total_mutations"] > 0


def test_mutations_local_mode_returns_400():
    """
    mode="local" est syntaxiquement valide (Literal), donc accepté par
    Pydantic (pas de 422), mais rejeté ensuite par detect_mutations()
    au niveau métier -> 400.
    """
    response = client.post(
        "/api/v1/mutations",
        json={"seq1": "TTTACCGTAAA", "seq2": "ACCGT", "mode": "local"},
    )

    assert response.status_code == 400
    assert "error" in response.json()


def test_mutations_incompatible_types_returns_400():
    response = client.post(
        "/api/v1/mutations",
        json={"seq1": "ATGCGTAA", "seq2": "MKTAYIAKQRQISFVK"},
    )

    assert response.status_code == 400


def test_mutations_missing_field_returns_422():
    response = client.post("/api/v1/mutations", json={"seq1": "ATGC"})

    assert response.status_code == 422