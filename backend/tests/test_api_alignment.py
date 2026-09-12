"""
Tests d'intégration HTTP pour /api/v1/alignment.
"""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_alignment_valid_global_request_returns_200():
    response = client.post(
        "/api/v1/alignment",
        json={"seq1": "ATGCGTAA", "seq2": "ATGCGT", "mode": "global"},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["mode"] == "global"
    assert len(data["aligned_seq1"]) == len(data["aligned_seq2"])


def test_alignment_valid_local_request_returns_200():
    response = client.post(
        "/api/v1/alignment",
        json={"seq1": "TTTACCGTAAA", "seq2": "ACCGT", "mode": "local"},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["mode"] == "local"
    assert data["aligned_seq1"] == "ACCGT"


def test_alignment_invalid_mode_returns_422():
    response = client.post(
        "/api/v1/alignment",
        json={"seq1": "ATGC", "seq2": "ATGC", "mode": "semi-global"},
    )

    assert response.status_code == 422


def test_alignment_incompatible_types_returns_400():
    response = client.post(
        "/api/v1/alignment",
        json={"seq1": "ATGCGTAA", "seq2": "MKTAYIAKQRQISFVK"},
    )

    assert response.status_code == 400
    assert "error" in response.json()


def test_alignment_missing_field_returns_422():
    response = client.post("/api/v1/alignment", json={"seq1": "ATGC"})

    assert response.status_code == 422