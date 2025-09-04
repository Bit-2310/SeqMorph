from __future__ import annotations
import json, time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Optional
import logging

@dataclass
class RunConfig:
    label: str
    seed: Optional[int] = None
    store: str = "StringStore"
    struct_rate_pct: Optional[float] = None
    point_rate: Optional[float] = None
    ti_tv_ratio: Optional[float] = None
    use_cpg_bias: Optional[bool] = None
    mean_seg_len: Optional[int] = None

def prepare_run_dir(base: str = "output", label: str = "run") -> Path:
    ts = time.strftime("%Y%m%d_%H%M%S")
    p = Path(base) / f"{ts}_{label}"
    p.mkdir(parents=True, exist_ok=True)
    return p

def write_manifest(run_dir: Path, cfg: RunConfig, extra: Optional[Dict[str, Any]] = None) -> None:
    data = asdict(cfg)
    if extra: data.update(extra)
    (run_dir / "manifest.json").write_text(json.dumps(data, indent=2))

def configure_file_logging(run_dir: Path) -> None:
    fh = logging.FileHandler(run_dir / "seqmorph.log")
    fh.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s")
    fh.setFormatter(fmt)
    root = logging.getLogger("seqmorph")
    root.setLevel(logging.INFO)
    root.addHandler(fh)
    logging.getLogger("seqmorph.mutation.structural").addHandler(fh)
