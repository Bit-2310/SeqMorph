from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Dict, List, Optional

from .types import StructuralEvent, OpName
from .structural import sample_structural_events, apply_structural_events_on_store
from .point import apply_point_mutations_on_store, ChooseBase

from sequence_store import BaseStore



@dataclass
class MutationReport:
    """Minimal mutation report; extend as needed."""
    point_mutations: int
    structural_events: List[StructuralEvent]


class MutationEngine:
    """
    Store-first mutation engine.
    - Uses a shared RNG (seed once for reproducibility).
    - Structural ops operate on BaseStore backends (ChunkedStore recommended).
    - Point mutations are injected via a callback policy.
    """

    def __init__(self, rng: Optional[random.Random] = None, seed: Optional[int] = None):
        self.rng = rng or random.Random()
        if seed is not None:
            self.rng.seed(seed)

    # ---- point mutations ----
    def mutate_points(
        self,
        store: BaseStore,
        positions: List[int],
        choose_base: ChooseBase,
    ) -> int:
        return apply_point_mutations_on_store(self.rng, store, positions, choose_base)

    # ---- structural mutations ----
    def mutate_structural(
        self,
        store: BaseStore,
        win_start: int,
        win_end: int,
        rate_pct: float,
        mean_seg_len: int = 500,
        mix: Optional[Dict[OpName, float]] = None,
    ) -> List[StructuralEvent]:
        events = sample_structural_events(
            self.rng, win_start, win_end, rate_pct, mean_seg_len, mix  # type: ignore[arg-type]
        )
        apply_structural_events_on_store(store, events)
        return events

    # ---- combined convenience ----
    def mutate(
        self,
        store: BaseStore,
        *,
        point_positions: Optional[List[int]] = None,
        choose_base: Optional[callable] = None,
        struct_window: Optional[tuple[int, int]] = None,
        rate_pct: float = 0.0,
        mean_seg_len: int = 500,
        mix: Optional[Dict[OpName, float]] = None,
    ) -> MutationReport:
        pm = 0
        if point_positions and choose_base:
            pm = self.mutate_points(store, point_positions, choose_base)

        se: List[StructuralEvent] = []
        if struct_window and rate_pct > 0.0:
            s, e = struct_window
            se = self.mutate_structural(store, s, e, rate_pct, mean_seg_len, mix)
        return MutationReport(point_mutations=pm, structural_events=se)
