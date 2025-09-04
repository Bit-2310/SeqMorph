from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

# Allowed structural operation names
OpName = Literal["inversion", "duplication", "translocation"]


@dataclass(frozen=True)
class StructuralEvent:
    """
    One structural edit on a half-open interval [start, end).

    Fields
    ------
    op : OpName
        Operation name.
    start, end : int
        0-based coordinates in the current sequence, end-exclusive.
    insert_pos : Optional[int]
        Required for duplication/translocation; absolute 0-based index.
    """
    op: OpName
    start: int
    end: int
    insert_pos: Optional[int] = None


import random
from typing import Protocol

class ChooseBase(Protocol):
    """Policy callable for point mutations."""
    def __call__(self, pos: int, current: str,
                 prev: Optional[str], nxt: Optional[str], rng: random.Random) -> str: ...
