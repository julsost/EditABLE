# rs3_service.py
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List
from rs3.seq import predict_seq
import math
import re

app = FastAPI(title="RS3 Scoring API", version="1.0.0")

# --- CORS (adjust origins if you prefer to restrict to your site) ---
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # e.g. ["https://<your-gh-username>.github.io"]
    allow_methods=["POST", "GET", "OPTIONS"],
    allow_headers=["*"],
)

# --- Models ---
class ScoreReq(BaseModel):
    sequences: List[str] = Field(..., description="30-mer DNA, A/C/G/T only")
    tracr: str = Field("Hsu2013", description="Hsu2013 or Chen2013")
    as_percent: bool = True  # return 0-100 if True, else 0-1

class ScoreResp(BaseModel):
    scores: List[float]  # aligned with input order

# --- Utils ---
NUC_RE = re.compile(r"^[ACGT]{30}$", re.IGNORECASE)

def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-float(x)))

# --- Routes ---
@app.get("/health")
def health():
    return {"ok": True}

@app.post("/score", response_model=ScoreResp)
def score(req: ScoreReq):
    # validate inputs conservatively
    seqs = []
    for s in req.sequences:
        s = (s or "").upper()
        if not NUC_RE.match(s):
            raise HTTPException(status_code=422, detail=f"Invalid 30-mer: {s}")
        seqs.append(s)

    raw = predict_seq(seqs, sequence_tracr=req.tracr)
    probs = [_sigmoid(x) for x in raw]
    out = [round(p * 100, 2) if req.as_percent else round(p, 6) for p in probs]
    return ScoreResp(scores=out)
