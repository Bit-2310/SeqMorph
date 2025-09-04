from __future__ import annotations

import logging
import random
from typing import Dict, List, Tuple, Optional
from mutation.types import OpName, StructuralEvent
from sequence_store import BaseStore


log = logging.getLogger("seqmorph.mutation.structural")


def _sample_non_overlapping_segments(
    rng: random.Random,
    win_start: int,
    win_end: int,
    target_total_len: int,
    mean_seg_len: int,
) -> List[Tuple[int, int]]:
    """
    Sample non-overlapping [s, e) segments inside [win_start, win_end) whose
    combined length is ~ target_total_len. Lengths are drawn around mean_seg_len.
    """
    segments: List[Tuple[int, int]] = []
    usable = max(0, win_end - win_start)
    remaining = max(0, min(target_total_len, usable))
    if usable <= 1 or remaining == 0:
        return segments

    # Avoid long loops on tiny windows
    max_trials = 10 * max(1, usable // max(1, mean_seg_len))
    trials = 0
    taken: List[Tuple[int, int]] = []

    while remaining > 0 and trials < max_trials:
        trials += 1
        seg_len = max(1, int(rng.gauss(mu=mean_seg_len, sigma=max(1, mean_seg_len // 3))))
        seg_len = min(seg_len, remaining, usable)
        s = rng.randrange(win_start, max(win_start + 1, win_end - seg_len + 1))
        e = s + seg_len

        # reject if overlaps any taken segment
        if any(not (e <= ts or te <= s) for ts, te in taken):
            continue

        taken.append((s, e))
        segments.append((s, e))
        remaining -= seg_len

    segments.sort()
    return segments


def sample_structural_events(
    rng: random.Random,
    win_start: int,
    win_end: int,
    rate_pct: float,
    mean_seg_len: int = 500,
    mix: Optional[Dict[OpName, float]] = None,
) -> List[StructuralEvent]:
    """
    Mixed-model sampler for structural edits inside [win_start, win_end).

    - rate_pct:     expected fraction (0..100) of the window to be affected
    - mean_seg_len: typical event segment length (nt)
    - mix:          op weights among {"inversion","duplication","translocation"}

    Returns a non-overlapping list of StructuralEvent with absolute coordinates.
    """
    if mix is None:
        mix = {"inversion": 0.34, "duplication": 0.33, "translocation": 0.33}

    win_len = max(0, win_end - win_start)
    if win_len <= 1 or rate_pct <= 0.0:
        return []

    # expected edited length
    target_total = int((rate_pct / 100.0) * win_len)
    if target_total <= 0:
        # probabilistic single event
        if rng.random() > (rate_pct / 100.0):
            return []
        target_total = min(mean_seg_len, win_len)

    segments = _sample_non_overlapping_segments(rng, win_start, win_end, target_total, mean_seg_len)
    if not segments:
        return []

    ops, weights = zip(*mix.items())
    events: List[StructuralEvent] = []
    for s, e in segments:
        op: OpName = rng.choices(ops, weights=weights, k=1)[0]  # type: ignore[assignment]
        if op == "inversion":
            events.append(StructuralEvent("inversion", s, e, None))
        elif op == "duplication":
            # insert right after the duplicated segment
            events.append(StructuralEvent("duplication", s, e, e))
        else:  # "translocation"
            # move to end of the window (simple and index-stable)
            events.append(StructuralEvent("translocation", s, e, win_end))
    return events


def apply_structural_events_on_store(store: BaseStore, events: List[StructuralEvent]) -> None:
    """
    Apply events to the store in-place.

    - We apply from right to left (descending by start) so earlier coordinates
      remain valid while we mutate the right side first.
    - We also defensively clamp events to store bounds and drop invalid ones.
    """
    if not events:
        return

    # Defensive bounds & shape check
    n = len(store)
    cleaned: List[StructuralEvent] = []
    for ev in events:
        s = max(0, min(ev.start, n))
        e = max(0, min(ev.end, n))
        if s >= e:
            continue
        ins = ev.insert_pos
        if ev.op in ("duplication", "translocation"):
            if ins is None or not (0 <= ins <= n):
                continue
        cleaned.append(StructuralEvent(ev.op, s, e, ins))
    if not cleaned:
        return

    for ev in sorted(cleaned, key=lambda x: (x.start, x.end), reverse=True):
        if ev.op == "inversion":
            store.invert(ev.start, ev.end)
        elif ev.op == "duplication":
            store.duplicate(ev.start, ev.end, ev.insert_pos)  # type: ignore[arg-type]
        elif ev.op == "translocation":
            store.translocate(ev.start, ev.end, ev.insert_pos)  # type: ignore[arg-type]
        else:  # pragma: no cover
            raise ValueError(f"Unknown structural op: {ev.op}")
        log.info(
            "struct-op=%s start=%d end=%d insert=%s new_len=%d",
            ev.op, ev.start, ev.end, ev.insert_pos, len(store)
        )
