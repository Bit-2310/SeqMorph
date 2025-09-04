from __future__ import annotations

import random
from typing import Iterable, Optional

from sequence_store import BaseStore
from mutation.types import ChooseBase


def apply_point_mutations_on_store(
    rng: random.Random,
    store: BaseStore,
    positions: Iterable[int],
    choose_base: ChooseBase,
) -> int:
    """
    Apply point mutations at specified positions using a provided policy.
    Returns the number of positions changed.
    """
    changed = 0
    for pos in positions:
        cur = store.get(pos)
        prev = store.get(pos - 1) if pos - 1 >= 0 else None
        nxt = store.get(pos + 1) if pos + 1 < len(store) else None
        new = choose_base(pos, cur, prev, nxt, rng)
        if new and new != cur:
            store.set(pos, new)
            changed += 1
    return changed
