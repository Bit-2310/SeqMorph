from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, Dict, Tuple, Sequence, Mapping
import random, bisect

DNA: set[str] = {"A", "C", "G", "T", "N"}

class ChooseBase(Protocol):
    def choose(self, left: str, base: str, right: str, rng: random.Random) -> str: ...

def _sanitize(b: str) -> str:
    if not b: return "N"
    b = b.upper()
    if b == "U": b = "T"
    return b if b in DNA else "N"

def _rc(b: str) -> str:
    return {"A":"T","T":"A","C":"G","G":"C","N":"N"}[b]

def _cum(items: Sequence[tuple[str,float]]) -> tuple[list[float], list[str]]:
    total = 0.0; cum: list[float] = []; labels: list[str] = []
    for lab,w in items:
        if w > 0:
            total += w; cum.append(total); labels.append(lab)
    return cum, labels

@dataclass(slots=True)
class ContextChooserConfig:
    # k-mer cascade: try 5-mer -> 3-mer -> baseline(1-mer)
    rc_symmetry: bool = True                  # enforce reverse-complement symmetry
    ts_tv_ratio: float = 1.0                  # >1 boosts transitions (A<->G, C<->T)
    # optional CpG bump layered on top of table (keep for convenience)
    cpg_enabled: bool = True
    cpg_strength: float = 3.0

class ContextAwareChooser(ChooseBase):
    """
    Data-driven k-mer chooser:
      - Looks up weight tables in order: 5-mer -> 3-mer -> baseline(1-mer)
      - Optional reverse-complement symmetry at load time
      - Optional transition/transversion re-weighting
      - Optional CpG bump (B=='C' and R=='G') layered after table weights
    """
    __slots__ = ("cfg","mat5","mat3","baseline")

    def __init__(
        self,
        cfg: ContextChooserConfig,
        matrix_5mer: Mapping[Tuple[str,str,str,str,str], Mapping[str,float]] | None,
        matrix_3mer: Mapping[Tuple[str,str,str], Mapping[str,float]] | None,
        baseline_1mer: Mapping[str, Mapping[str,float]],
    ):
        self.cfg = cfg
        self.mat5 = dict(matrix_5mer or {})
        self.mat3 = dict(matrix_3mer or {})
        self.baseline = {k: dict(v) for k, v in baseline_1mer.items()}

        if self.cfg.rc_symmetry:
            # add reverse-complement entries if missing
            # 5-mer
            m5 = {}
            for key, dist in self.mat5.items():
                rc_key = tuple(_rc(x) for x in key[::-1])
                if rc_key not in self.mat5:
                    m5[rc_key] = { _rc(b): w for b,w in dist.items() }
            self.mat5.update(m5)
            # 3-mer
            m3 = {}
            for key, dist in self.mat3.items():
                rc_key = tuple(_rc(x) for x in key[::-1])
                if rc_key not in self.mat3:
                    m3[rc_key] = { _rc(b): w for b,w in dist.items() }
            self.mat3.update(m3)

    @staticmethod
    def _is_transition(orig: str, mut: str) -> bool:
        return (orig, mut) in {("A","G"),("G","A"),("C","T"),("T","C")}

    def _reweight_tstv(self, base: str, weights: Dict[str,float]) -> None:
        r = self.cfg.ts_tv_ratio
        if r == 1.0: return
        for b in list(weights.keys()):
            if self._is_transition(base, b):
                weights[b] *= r

    def choose(self, left: str, base: str, right: str, rng: random.Random) -> str:
        L = _sanitize(left); B = _sanitize(base); R = _sanitize(right)

        # Lookup cascade: 5-mer -> 3-mer -> baseline(1-mer)
        # For 5-mer we synthesize N-L-B-R-N to avoid extra indexing
        probs = self.mat5.get(("N", L, B, R, "N"))
        if probs is None:
            probs = self.mat3.get((L, B, R))
        if probs is None:
            probs = self.baseline.get(B, {})

        # start from table weights (copy so we can tweak)
        cand = [x for x in ("A","C","G","T") if x != B]
        weights: Dict[str,float] = {b: float(probs.get(b, 0.0)) for b in cand}

        # optional TT/TV scaling
        self._reweight_tstv(B, weights)

        # optional CpG bump layered after
        if self.cfg.cpg_enabled and B == "C" and R == "G":
            weights["T"] = weights.get("T", 0.0) * max(self.cfg.cpg_strength, 0.0001)

        # sample
        items = [(b,w) for b,w in weights.items() if w > 0]
        if not items:
            return B
        cum, labels = _cum(items)
        r = rng.random() * cum[-1]
        i = bisect.bisect_left(cum, r)
        if i >= len(labels): i = len(labels) - 1
        return labels[i]
