from __future__ import annotations
import logging
from typing import List
from sequence_store.base import BaseStore

logger = logging.getLogger("seqmorph.store")

class StringStore(BaseStore):
    """
    Baseline storage using a Python list of single-character strings.
    Complexity is O(n) for structural edits, but the API is stable.
    We’ll swap this for a chunked/Fenwick-backed store next.
    """

    __slots__ = ("_data",)

    def __init__(self, sequence: str = "") -> None:
        self._data: List[str] = list(sequence or "")

    def __len__(self) -> int:
        return len(self._data)

    # ----- views -----
    def to_string(self) -> str:
        return "".join(self._data)

    def slice(self, start: int, end: int) -> str:
        return "".join(self._data[start:end])

    def get(self, pos: int) -> str:
        return self._data[pos]

    # ----- point edits -----
    def set(self, pos: int, base: str) -> None:
        b = (base or "")[:1]
        if not b:
            raise ValueError("Base must be a non-empty character.")
        old = self._data[pos]
        self._data[pos] = b
        logger.debug("set(pos=%d, %s→%s)", pos, old, b)

    def replace_range(self, start: int, end: int, substring: str) -> None:
        self._data[start:end] = list(substring or "")
        logger.debug("replace_range([%d:%d], len=%d)", start, end, len(substring or ""))

    def insert(self, pos: int, substring: str) -> None:
        if not substring:
            return
        self._data[pos:pos] = list(substring)
        logger.debug("insert(pos=%d, len=%d)", pos, len(substring))

    def delete(self, start: int, end: int) -> None:
        if end <= start:
            return
        del self._data[start:end]
        logger.debug("delete([%d:%d])", start, end)

    # ----- structural edits -----
    def invert(self, start: int, end: int) -> None:
        if end - start <= 1:
            return
        seg = self._data[start:end]
        self._data[start:end] = seg[::-1]
        logger.info("invert([%d:%d]) len=%d", start, end, end - start)

    def duplicate(self, start: int, end: int, insert_pos: int) -> None:
        if end <= start:
            return
        seg = self._data[start:end]
        self._data[insert_pos:insert_pos] = seg
        logger.info("duplicate([%d:%d]→%d) len=%d", start, end, insert_pos, len(seg))

    def translocate(self, start: int, end: int, insert_pos: int) -> None:
        if end <= start:
            return
        seg = self._data[start:end]
        # Remove original
        del self._data[start:end]
        # Adjust insertion index if we removed before it
        if insert_pos > start:
            insert_pos -= (end - start)
        self._data[insert_pos:insert_pos] = seg
        logger.info("translocate([%d:%d]→%d) len=%d", start, end, insert_pos, len(seg))
