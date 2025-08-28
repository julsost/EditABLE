# rs3_service.py
from fastapi import FastAPI
from pydantic import BaseModel, Field
from typing import List
from rs3.seq import predict_seq
import math

app = FastAPI(title="RS3 Scoring API")

def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-float(x)))

class ScoreReq(BaseModel):
    sequences: List[str] = Field(..., description="30-mer DNA contexts, A/C/G/T, uppercase")
    tracr: str = Field("Hsu2013", description="Hsu2013 or Chen2013")
    as_percent: bool = True  # return 0-100 if True, else 0-1

class ScoreResp(BaseModel):
    scores: List[float]  # aligned to input order

@app.post("/score", response_model=ScoreResp)
def score(req: ScoreReq):
    seqs = [s.upper() for s in req.sequences]
    raw = predict_seq(seqs, sequence_tracr=req.tracr)
    probs = [sigmoid(x) for x in raw]
    out = [round(p * 100, 2) if req.as_percent else round(p, 6) for p in probs]
    return ScoreResp(scores=out)
