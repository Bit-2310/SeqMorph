# SeqMorph_Main.py — Web-ready backend (FastAPI)
from __future__ import annotations

import logging
import random
import time
from collections import Counter
from typing import Dict, Optional, Literal, List

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field, conint, confloat

from library import SequenceValidation
from input_module import SequenceFetcher, is_ncbi_accession, is_uniprot_accession
from sequence_structure import SequenceStructure
from sequence_store import StringStore, ChunkedStore
from mutation import MutationEngine, ContextChooserConfig
from analysis_module import SequenceAnalysisReport
from utils_run import prepare_run_dir, write_manifest, configure_file_logging, RunConfig


# -----------------------------------------------------------------------------
# App + globals
# -----------------------------------------------------------------------------
app = FastAPI(title="SeqMorph API", version="0.3")

def _store_factory(seq: str):
    # Use ChunkedStore for longer inputs; tweak thresholds as you profile
    if len(seq) >= 50_000:
        return ChunkedStore(seq, min_chunk=1024, max_chunk=4096)
    return StringStore(seq)

# Allow duplicate IDs via suffixing so you can fetch/add multiple variants safely
seqs = SequenceStructure(on_duplicate="suffix", store_factory=_store_factory)
fetcher = SequenceFetcher(cache_file="sequence_cache.json", timeout=15)

logger = logging.getLogger("seqmorph")
if not logger.handlers:
    logging.basicConfig(level=logging.INFO)


# -----------------------------------------------------------------------------
# Pydantic models
# -----------------------------------------------------------------------------
class SequenceFetchRequest(BaseModel):
    accession_id: str = Field(..., description="NCBI or UniProt accession")

class SequenceAddRequest(BaseModel):
    accession_id: Optional[str] = Field(
        None, description="Custom ID; if omitted, a sequential ID is used"
    )
    sequence: str = Field(..., description="Raw sequence letters")
    seq_type: Optional[Literal["DNA", "RNA", "Protein"]] = None
    overwrite: bool = False
    description: Optional[str] = None

class MutateAnalyzeRequest(BaseModel):
    accession_id: str
    struct_rate: confloat(ge=0.0, le=100.0) = 0.0
    point_rate: confloat(ge=0.0, le=100.0) = 0.0
    ti_tv_ratio: confloat(ge=0.0) = 2.1
    use_cpg_bias: bool = True
    mean_seg_len: conint(ge=1) = 500
    start: conint(ge=1)
    end: Optional[int] = None
    seed: Optional[int] = None
    save_outputs: bool = False

# Additional, richer response bits
class RunParams(BaseModel):
    window: List[int]
    struct_rate: float
    point_rate: float
    ti_tv_ratio: float
    use_cpg_bias: bool
    mean_seg_len: int
    seed: Optional[int] = None
    store: str

class EventPreview(BaseModel):
    op: Literal["inversion", "duplication", "translocation"]
    start: int
    end: int
    length: int
    insert_pos: Optional[int] = None

class MutateAnalyzeResponse(BaseModel):
    accession_id: str
    seq_type: str
    original_length: int
    mutated_length: int
    delta_length: int
    num_events: int
    num_point_mutations: int
    struct_event_counts: Dict[str, int]
    struct_dup_length_added: int
    params: RunParams
    events_preview: List[EventPreview]
    duration_ms: float
    report: Dict
    output_paths: Optional[Dict] = None  # preserved for compatibility


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


# -----------------------------------------------------------------------------
# Routes
# -----------------------------------------------------------------------------
@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "ok"}


@app.post("/sequence/fetch")
def sequence_fetch(payload: SequenceFetchRequest) -> Dict:
    acc = payload.accession_id.strip()
    if not (is_ncbi_accession(acc) or is_uniprot_accession(acc)):
        raise HTTPException(400, "Accession ID doesn't look like NCBI or UniProt.")
    try:
        seq = fetcher.fetch_sequence(acc)
        seq_type = SequenceValidation.detect_seq(seq)
        used_id = seqs.add_sequence(
            acc, seq, seq_type=seq_type, overwrite=False, metadata={"source": "remote"}
        )
        return {"accession_id": used_id, "seq_type": seq_type, "length": len(seq)}
    except Exception as e:
        raise HTTPException(502, f"Failed to fetch sequence: {e}")


@app.post("/sequence/add")
def sequence_add(payload: SequenceAddRequest) -> Dict:
    seq = (payload.sequence or "").strip().upper()
    if not seq:
        raise HTTPException(400, "Sequence cannot be empty.")
    try:
        stype = payload.seq_type or SequenceValidation.detect_seq(seq)
        if not SequenceValidation.validate(seq, stype):  # type: ignore[arg-type]
            raise HTTPException(400, f"Invalid {stype} sequence.")
        acc = payload.accession_id or f"seq_{len(seqs.list_accessions()) + 1}"
        used_id = seqs.add_sequence(
            acc,
            seq,
            seq_type=stype,
            overwrite=payload.overwrite,
            metadata={"description": payload.description or ""},
        )
        return {"accession_id": used_id, "seq_type": stype, "length": len(seq)}
    except Exception as e:
        raise HTTPException(400, f"Failed to add sequence: {e}")


@app.post("/mutate-and-analyze", response_model=MutateAnalyzeResponse)
def mutate_and_analyze(payload: MutateAnalyzeRequest) -> MutateAnalyzeResponse:
    t0 = time.perf_counter()

    acc = payload.accession_id.strip()
    if acc not in seqs.sequences:
        raise HTTPException(
            404, f"Unknown accession '{acc}'. Use /sequence/fetch or /sequence/add first."
        )

    store = seqs.sequences[acc]
    original = store.to_string()
    seq_meta = seqs.get_metadata(acc)
    seq_type: str = seq_meta.get("type", "DNA")
    n = len(store)

    # 1-based -> 0-based half-open
    start_1 = payload.start
    end_1 = payload.end or n
    if end_1 < start_1:
        raise HTTPException(400, "end must be >= start")
    start0 = max(0, min(n, start_1 - 1))
    end0 = max(0, min(n, end_1))

    # --- Mutations ---
    chooser_cfg = ContextChooserConfig(
        ts_tv_ratio=float(payload.ti_tv_ratio), cpg_enabled=payload.use_cpg_bias
    )
    eng = MutationEngine(seed=payload.seed, chooser_cfg=chooser_cfg)

    # Point mutations
    num_pm = 0
    if payload.point_rate > 0.0 and start0 < end0:
        win_len = end0 - start0
        num_to_mutate = int(win_len * (payload.point_rate / 100.0))
        if num_to_mutate > 0:
            pos_to_mutate = random.sample(
                range(start0, end0), k=min(num_to_mutate, win_len)
            )
            num_pm = eng.mutate_points(store, pos_to_mutate)

    # Structural mutations (store-first)
    events = []
    if payload.struct_rate > 0.0 and start0 < end0:
        events = eng.mutate_structural(
            store,
            win_start=start0,
            win_end=end0,
            rate_pct=float(payload.struct_rate),
            mean_seg_len=int(payload.mean_seg_len),
        )

    mutated = store.to_string()

    # Enriched metrics
    counts = Counter(e.op for e in events)
    dup_added = sum(e.end - e.start for e in events if e.op == "duplication")
    events_preview = [
        EventPreview(
            op=e.op, start=e.start, end=e.end, length=(e.end - e.start), insert_pos=e.insert_pos
        )
        for e in events[: min(10, len(events))]
    ]
    delta_length = len(mutated) - len(original)

    # Optional outputs
    out_paths: Optional[Dict[str, str]] = None
    if payload.save_outputs:
        run_dir = prepare_run_dir(label="mut")
        configure_file_logging(run_dir)
        write_manifest(
            run_dir,
            RunConfig(
                label="mut",
                seed=payload.seed,
                store=store.__class__.__name__,
                struct_rate_pct=float(payload.struct_rate),
                point_rate=float(payload.point_rate),
                ti_tv_ratio=float(payload.ti_tv_ratio),
                use_cpg_bias=payload.use_cpg_bias,
                mean_seg_len=int(payload.mean_seg_len),
            ),
            extra={
                "accession_id": acc,
                "num_events": len(events),
                "num_point_mutations": num_pm,
                "window": [start0, end0],
            },
        )
        # FASTA + events
        (run_dir / f"{acc}_original.fasta").write_text(
            f">{acc} original\n{_wrap_fasta(original)}\n"
        )
        (run_dir / f"{acc}_mutated.fasta").write_text(
            f">{acc} mutated\n{_wrap_fasta(mutated)}\n"
        )
        import json as _json
        evt_dicts = [
            {"op": e.op, "start": e.start, "end": e.end, "insert_pos": e.insert_pos}
            for e in events
        ]
        (run_dir / f"{acc}_events.json").write_text(_json.dumps(evt_dicts, indent=2))
        out_paths = {
            "run_dir": str(run_dir),
            "original_fasta": str(run_dir / f"{acc}_original.fasta"),
            "mutated_fasta": str(run_dir / f"{acc}_mutated.fasta"),
            "events_json": str(run_dir / f"{acc}_events.json"),
            "log": str(run_dir / "seqmorph.log"),
        }

    t1 = time.perf_counter()

    return MutateAnalyzeResponse(
        accession_id=acc,
        seq_type=seq_type,
        original_length=len(original),
        mutated_length=len(mutated),
        delta_length=delta_length,
        num_events=len(events),
        num_point_mutations=num_pm,
        struct_event_counts=dict(counts),
        struct_dup_length_added=dup_added,
        params=RunParams(
            window=[start0, end0],
            struct_rate=float(payload.struct_rate),
            point_rate=float(payload.point_rate),
            ti_tv_ratio=float(payload.ti_tv_ratio),
            use_cpg_bias=payload.use_cpg_bias,
            mean_seg_len=int(payload.mean_seg_len),
            seed=payload.seed,
            store=store.__class__.__name__,
        ),
        events_preview=events_preview,
        duration_ms=(t1 - t0) * 1000.0,
        report=SequenceAnalysisReport(
            original, mutated, sequence_type=seq_type, enable_alignment=False
        ).compare_sequences(),
        output_paths=out_paths,
    )


# -----------------------------------------------------------------------------
# Dev runner
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Run with: python src/SeqMorph_Main.py
    import uvicorn
    uvicorn.run(
        "SeqMorph_Main:app",
        host="127.0.0.1",
        port=8000,
        reload=True,    # handy during active dev; turn off in prod
        log_level="info",
    )
