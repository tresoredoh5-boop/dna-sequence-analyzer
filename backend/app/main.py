"""
Point d'entrée de l'API BioAI Sequence Analyzer.

Ce fichier ne contient aucune logique métier : il assemble l'application
FastAPI, branche les routers, la configuration CORS, et le gestionnaire
d'erreurs global qui convertit les ValueError métier en HTTP 400.
"""

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from app.api import analysis_routes, alignment_routes, mutation_routes

app = FastAPI(
    title="BioAI Sequence Analyzer API",
    description="API d'analyse de séquences biologiques (ADN, ARN, protéines).",
    version="1.0.0",
)

# CORS : origines explicitement autorisées, pas de "*".
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:5173"],
    allow_credentials=True,
    allow_methods=["POST", "GET"],
    allow_headers=["Content-Type"],
)


@app.exception_handler(ValueError)
def handle_value_error(request: Request, exc: ValueError) -> JSONResponse:
    """
    Convertit toute ValueError levée par le cœur bioinformatique
    (séquence vide, types incompatibles, mode invalide...) en une
    réponse HTTP 400 avec le message d'erreur d'origine.
    """
    return JSONResponse(status_code=400, content={"error": str(exc)})


app.include_router(analysis_routes.router)
app.include_router(alignment_routes.router)
app.include_router(mutation_routes.router)


@app.get("/health", tags=["health"])
def health_check() -> dict:
    """Endpoint simple pour vérifier que l'API répond."""
    return {"status": "ok"}