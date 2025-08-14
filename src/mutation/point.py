from __future__ import annotations

import random
from typing import Iterable, Optional, Protocol

from sequence_store import BaseStore



class ChooseBase(Protocol):
    """Policy callable for point mutations."""
    def __call__(self, pos: int, current: str,
                 prev: Optional[str], nxt: Optional[str]) -> str: ...


def apply_point_mutations_on_store(
    rng: random.Random,
    store: BaseStore,
    positions: Iterable[int],
    choose_base: ChooseBase,
) -> int:
    """
    Apply point mutations at specified positions using a provided policy.
    Returns the number of positions changed.

    Note: 'rng' is supplied for policies that capture it via closure; this
    function itself does not use it directly.
    """
    changed = 0
    for pos in positions:
        cur = store.get(pos)
        prev = store.get(pos - 1) if pos - 1 >= 0 else None
        nxt = store.get(pos + 1) if pos + 1 < len(store) else None
        new = choose_base(pos, cur, prev, nxt)
        if new and new != cur:
            store.set(pos, new)
            changed += 1
    return changed
