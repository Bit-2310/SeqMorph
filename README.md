# SeqMorph

SeqMorph is a bioinformatics toolkit for simulating mutations in DNA/RNA/Protein sequences and comparing the results. It focuses on **correctness**, **scalability**, and **reproducibility** using a store-first design (pluggable sequence backends) and a small FastAPI service.

---

## Key features (current)

- 🧬 **Sequence storage, not just strings**
  - `SequenceStructure` registry + `BaseStore` backends:
    - `StringStore` for small inputs
    - `ChunkedStore` for long sequences (O(log k) index mapping)
- 🔁 **Structural mutations** (reproducible with a seed)
  - inversion, duplication, translocation (windowed; non-overlapping sampling)
  - event mix/weights & mean segment length
- 🧪 **Sequence analysis**
  - GC%, k-mers, Shannon entropy, ORFs/translation (RNA handled via U→T normalization)
  - comparison report between original & mutated sequences
- 📦 **Reproducibility & outputs**
  - `manifest.json` (seed, params, backend), per-event logs, FASTA export, `events.json`
- 🌐 **API-first**
  - add/fetch sequences and run mutations via HTTP
  - easy to drive from notebooks, CLI wrappers, or a future GUI

> Note: Point mutations with context policy are designed (engine supports them) and will be exposed via API once we surface the policy knobs.

---

## What’s new

- ✅ Unified `library/` (replaces `gene_library.py` and `protein_library.py`)
- ✅ New `mutation/` package (replaces `mutation_module.py`)
- ✅ FastAPI backend (`SeqMorph_Main.py`) instead of a monolithic CLI
- ✅ Store-first sequence handling (`sequence_store/` backends)
- ✅ Correct CpG “boost” logic (applies only when the **current base is C** and **next is G**)
- 🗑️ Old CLI/GUI temporarily removed; they’ll be rebuilt on the API

---

## Install

Requires **Python 3.11+**.

```bash
# From repo root
pip install -r requirements.txt
```

`requirements.txt` (for reference):

```
biopython>=1.83
requests>=2.31
fastapi>=0.110
pydantic>=2.5
uvicorn[standard]>=0.27
scipy>=1.11
matplotlib>=3.8
```

---

## Quick start (run the API)

```bash
# macOS/Linux
export PYTHONPATH=./src
python src/SeqMorph_Main.py

# Windows PowerShell
$env:PYTHONPATH = "$PWD\src"
python .\src\SeqMorph_Main.py
```

Open **http://127.0.0.1:8000/docs** for interactive docs.

### Example flow (via API)

1) **Add a sequence**
```bash
curl -X POST http://127.0.0.1:8000/sequence/add   -H "Content-Type: application/json"   -d '{"sequence":"ATG" + "ACGT"*50}'
```

2) **Mutate & analyze**
```bash
curl -X POST http://127.0.0.1:8000/mutate-and-analyze   -H "Content-Type: application/json"   -d '{
    "accession_id": "seq_1",
    "struct_rate": 3.0,
    "mean_seg_len": 200,
    "start": 1,
    "seed": 123,
    "save_outputs": true
  }'
```

Response includes: original/mutated lengths, **delta**, per-op counts, total duplicated length added, **full structural events**, a 10-event preview, params (window/rate/seed/store), duration, and an analysis report. When `save_outputs=true`, files are written under `./output/<timestamp>_mut/`:
- `*_original.fasta`, `*_mutated.fasta`
- `*_events.json`
- `manifest.json`
- `seqmorph.log`

---

## DIY “smoke test” (no script in repo)

You can quickly sanity-check the core without any extra files. Pick A (API-only) or B (direct Python).

### A) API-only smoke (works anywhere)

1) Start the server (see **Quick start**).
2) In a new terminal, add & mutate a random DNA:

**macOS/Linux**
```bash
SEQ=$(python - <<'PY'
import random; random.seed(1)
print(''.join(random.choice('ACGT') for _ in range(20000)))
PY
)
curl -s -X POST http://127.0.0.1:8000/sequence/add   -H "Content-Type: application/json"   -d "{"sequence":"$SEQ"}"

curl -s -X POST http://127.0.0.1:8000/mutate-and-analyze   -H "Content-Type: application/json"   -d '{"accession_id":"seq_1","struct_rate":3.0,"mean_seg_len":200,"start":1,"seed":123,"save_outputs":true}'   | python -m json.tool
```

**Windows PowerShell**
```powershell
$seq = python - <<'PY'
import random; random.seed(1)
print(''.join(random.choice('ACGT') for _ in range(20000)))
PY

curl -s -X POST http://127.0.0.1:8000/sequence/add `
  -H "Content-Type: application/json" `
  -d "{""sequence"":""$seq""}"

curl -s -X POST http://127.0.0.1:8000/mutate-and-analyze `
  -H "Content-Type: application/json" `
  -d "{""accession_id"":""seq_1"",""struct_rate"":3.0,""mean_seg_len"":200,""start"":1,""seed"":123,""save_outputs"":true}" `
  | python -m json.tool
```

**What to check**
- `num_events` > 0
- `delta_length == struct_dup_length_added` (inversions/translocations keep length)
- `output_paths.run_dir` exists with FASTA/JSON/log files
- Repeat with the same seed → identical `struct_events`

### B) Direct Python (import modules)

> This checks the storage and mutation engine *without* the HTTP layer.

**macOS/Linux**
```bash
PYTHONPATH=./src python - <<'PY'
from sequence_store import StringStore, ChunkedStore
from sequence_structure import SequenceStructure
from mutation import MutationEngine
import random

rng = random.Random(123)
seq = ''.join(rng.choice('ACGT') for _ in range(20000))

# auto-pick store
def store_factory(s: str):
    return ChunkedStore(s, min_chunk=1024, max_chunk=4096) if len(s) >= 50000 else StringStore(s)

ss = SequenceStructure(store_factory=store_factory)
acc = ss.add_sequence("acc1", seq, seq_type="DNA")
store = ss.sequences[acc]

eng = MutationEngine(seed=123)
events = eng.mutate_structural(store, win_start=0, win_end=len(store), rate_pct=3.0, mean_seg_len=200)
print("events:", len(events), "orig_len:", len(seq), "mut_len:", len(store))
dup_added = sum(e.end - e.start for e in events if e.op == "duplication")
print("delta:", len(store)-len(seq), "dup_added:", dup_added)
PY
```

**Windows PowerShell**
```powershell
$env:PYTHONPATH = "$PWD\src"
python - <<'PY'
from sequence_store import StringStore, ChunkedStore
from sequence_structure import SequenceStructure
from mutation import MutationEngine
import random

rng = random.Random(123)
seq = ''.join(rng.choice('ACGT') for _ in range(20000))

def store_factory(s: str):
    return ChunkedStore(s, min_chunk=1024, max_chunk=4096) if len(s) >= 50000 else StringStore(s)

ss = SequenceStructure(store_factory=store_factory)
acc = ss.add_sequence("acc1", seq, seq_type="DNA")
store = ss.sequences[acc]

eng = MutationEngine(seed=123)
events = eng.mutate_structural(store, win_start=0, win_end=len(store), rate_pct=3.0, mean_seg_len=200)
print("events:", len(events), "orig_len:", len(seq), "mut_len:", len(store))
dup_added = sum(e.end - e.start for e in events if e.op == "duplication")
print("delta:", len(store)-len(seq), "dup_added:", dup_added)
PY
```

---

## Project layout (core)

```
src/
  analysis_module.py          # analysis + comparison report (headless-safe)
  input_module.py             # FASTA I/O + NCBI/UniProt fetching
  library/                    # unified analysis/validation/utils
  mutation/                   # MutationEngine + structural/point primitives
  sequence_store/             # BaseStore + StringStore + ChunkedStore
  sequence_structure.py       # registry + metadata + file import/export
  SeqMorph_Main.py            # FastAPI app (backend “brain”)
```

---

## Reproducibility

- Pass `seed` in `/mutate-and-analyze` to reproduce events deterministically.
- Runs can write:
  - `manifest.json` (seed, backend, rate, mean segment, window, counts)
  - `seqmorph.log` (one line per structural event applied)

---

## Troubleshooting

- **ModuleNotFoundError for local imports**  
  Ensure `PYTHONPATH=./src` (see Quick start). Also make sure `src/__init__.py` **does not exist** (we treat `src` as a source root, not a package).

- **Biopython errors**  
  Install with `pip install biopython`. FASTA parsing depends on it.

- **RNA translation yields X**  
  We normalize U→T before translation and handle reverse-complement with the correct alphabet.

---

## Roadmap

- Expose **point-mutation policies** via API (context/CpG, ti/tv)
- Harden `ChunkedStore` rebalance thresholds with profiling
- New CLI & GUI (call the API)
- Optional packed DNA store (2-bit) for large datasets
- Test suite & benchmarks

---

## Contributing

We welcome issues and PRs. Focus areas:

- `sequence_store/` (performance, packed encodings)
- `mutation/` (policies, additional structural ops)
- `library/analysis.py` (metrics & scalable χ² strategies)
- `SeqMorph_Main.py` (new endpoints, params, CORS for future GUI)

Basic steps:
1. Fork & branch (`feature/xyz`)
2. Keep code PEP8 + typed; add docstrings
3. Prefer backend-shared logic (CLI/GUI must not duplicate logic)
4. Open a PR with a short rationale and example inputs/outputs

---

## License

MIT (see `LICENSE`).