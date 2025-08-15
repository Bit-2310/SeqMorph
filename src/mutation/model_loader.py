from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Mapping, Any
import json

K3 = Tuple[str,str,str]
K5 = Tuple[str,str,str,str,str]

@dataclass(slots=True)
class LoadedModel:
    name: str
    version: str
    k: int
    baseline: Dict[str, Dict[str,float]]
    ctx3: Dict[K3, Dict[str,float]]
    ctx5: Dict[K5, Dict[str,float]]

def load_model(path: str | Path) -> LoadedModel:
    p = Path(path)
    obj: Dict[str, Any] = json.loads(p.read_text(encoding="utf-8"))
    name = obj.get("name","model"); version = obj.get("version","0")
    k = int(obj.get("k",3))
    baseline = {k.upper(): {kk.upper(): float(vv) for kk,vv in v.items()}
                for k,v in obj["baseline"].items()}
    ctx3: Dict[K3, Dict[str,float]] = {}
    for k3, dist in obj.get("contexts_3mer", {}).items():
        L,B,R = k3.split()
        ctx3[(L.upper(),B.upper(),R.upper())] = {b.upper(): float(w) for b,w in dist.items()}
    ctx5: Dict[K5, Dict[str,float]] = {}
    for k5, dist in obj.get("contexts_5mer", {}).items():
        a,b,c,d,e = k5.split()
        ctx5[(a.upper(),b.upper(),c.upper(),d.upper(),e.upper())] = {b2.upper(): float(w) for b2,w in dist.items()}
    return LoadedModel(name, version, k, baseline, ctx3, ctx5)
