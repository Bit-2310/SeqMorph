from __future__ import annotations

import logging
from typing import List

from sequence_store.base import BaseStore

logger = logging.getLogger("seqmorph.store")


class ChunkedStore(BaseStore):
    """
    Chunked, mutable sequence storage with O(log k) index mapping,
    where k = number of chunks (k << sequence length).
    Edits touch only the affected chunks, then we rebalance locally.

    Representation
    --------------
    - self._chunks: List[List[str]]  # each chunk is a list of single-char strings
    - self._bit: Fenwick tree (1-based) over chunk lengths for fast prefix sums
    - self._min_chunk, self._max_chunk: soft bounds for (re)balancing

    Complexity (amortized)
    ----------------------
    - Index → (chunk, offset): O(log k)
    - Point set/get:          O(log k)
    - Insert/Delete:          O(log k + m) where m = length inserted/deleted
    - Invert/Duplicate/Txn:   O(log k + m) due to copy of the affected region

    Notes
    -----
    - This is a solid baseline. We rebuild the Fenwick tree whenever the
      number of chunks changes (cheap relative to sequence length).
    - Later, we can add 2-bit packing by swapping chunk payloads.
    """

    __slots__ = ("_chunks", "_bit", "_min_chunk", "_max_chunk")

    def __init__(self, sequence: str = "", *, min_chunk: int = 1024, max_chunk: int = 4096) -> None:
        if min_chunk <= 0 or max_chunk <= 0 or min_chunk > max_chunk:
            raise ValueError("min_chunk and max_chunk must be positive and min_chunk <= max_chunk")
        self._min_chunk = int(min_chunk)
        self._max_chunk = int(max_chunk)
        self._chunks: List[List[str]] = []
        if sequence:
            self._chunks = [list(sequence)]
        self._bit: List[int] = []
        self._rebuild_bit()

    # ------------------------------------------------------------------ #
    # Protocol: basic views
    # ------------------------------------------------------------------ #
    def __len__(self) -> int:
        return self._prefix_sum(len(self._chunks))

    def to_string(self) -> str:
        # faster than joining chars individually
        return "".join("".join(chk) for chk in self._chunks)

    def slice(self, start: int, end: int) -> str:
        self._bounds_check_closed(start, end)
        if end <= start:
            return ""
        out: List[str] = []
        ci_s, off_s = self._locate(start)
        ci_e, off_e = self._locate(end)
        if ci_s == ci_e:
            out.append("".join(self._chunks[ci_s][off_s:off_e]))
        else:
            out.append("".join(self._chunks[ci_s][off_s:]))
            for ci in range(ci_s + 1, ci_e):
                out.append("".join(self._chunks[ci]))
            out.append("".join(self._chunks[ci_e][:off_e]))
        return "".join(out)

    def get(self, pos: int) -> str:
        ci, off = self._locate_strict(pos)
        return self._chunks[ci][off]

    # ------------------------------------------------------------------ #
    # Protocol: point & range edits
    # ------------------------------------------------------------------ #
    def set(self, pos: int, base: str) -> None:
        if not base:
            raise ValueError("Base must be a non-empty character")
        ci, off = self._locate_strict(pos)
        old = self._chunks[ci][off]
        self._chunks[ci][off] = base[0]
        logger.debug("set(pos=%d, %s→%s)", pos, old, base[0])

    def replace_range(self, start: int, end: int, substring: str) -> None:
        # delete then insert
        self.delete(start, end)
        if substring:
            self.insert(start, substring)

    def insert(self, pos: int, substring: str) -> None:
        if not substring:
            return
        # Fast path: empty sequence
        if not self._chunks:
            self._chunks = [list(substring)]
            self._rebuild_bit()
            self._rebalance_around(0)
            logger.info("insert(pos=%d, len=%d) [init]", pos, len(substring))
            return

        if pos == len(self):
            # append at end
            self._chunks.append(list(substring))
            self._rebuild_bit()
            self._rebalance_around(len(self._chunks) - 1)
            logger.info("insert(pos=end, len=%d)", len(substring))
            return

        ci, off = self._locate(pos)
        # split the chunk at insertion point
        left = self._chunks[ci][:off]
        right = self._chunks[ci][off:]
        new_chunk = list(substring)
        # replace current chunk with left, then insert new, then right
        self._chunks[ci:ci + 1] = [left, new_chunk, right]
        self._rebuild_bit()
        # rebalance around the 3 affected positions
        self._rebalance_around(max(0, ci - 1))
        self._rebalance_around(ci)
        self._rebalance_around(min(len(self._chunks) - 1, ci + 1))
        logger.info("insert(pos=%d, len=%d)", pos, len(substring))

    def delete(self, start: int, end: int) -> None:
        self._bounds_check_closed(start, end)
        if end <= start:
            return
        ci_s, off_s = self._locate(start)
        ci_e, off_e = self._locate(end)

        if ci_s == ci_e:
            # delete within a single chunk
            chunk = self._chunks[ci_s]
            del chunk[off_s:off_e]
            self._rebuild_bit()
            self._rebalance_around(ci_s)
        else:
            # trim tail of start chunk
            del self._chunks[ci_s][off_s:]
            # trim head of end chunk
            del self._chunks[ci_e][:off_e]
            # drop full chunks in between
            if ci_e - ci_s > 1:
                del self._chunks[ci_s + 1:ci_e]
            self._rebuild_bit()
            # rebalance around the joined region
            self._rebalance_around(max(0, ci_s - 1))
            self._rebalance_around(min(len(self._chunks) - 1, ci_s))
        logger.info("delete([%d:%d])", start, end)

    # ----- structural ops implemented via slice/replace ----- #
    def invert(self, start: int, end: int) -> None:
        seg = self.slice(start, end)
        self.replace_range(start, end, seg[::-1])
        logger.info("invert([%d:%d]) len=%d", start, end, end - start)

    def duplicate(self, start: int, end: int, insert_pos: int) -> None:
        seg = self.slice(start, end)
        self.insert(insert_pos, seg)
        logger.info("duplicate([%d:%d]→%d) len=%d", start, end, insert_pos, end - start)

    def translocate(self, start: int, end: int, insert_pos: int) -> None:
        seg_len = max(0, end - start)
        if seg_len == 0:
            return
        seg = self.slice(start, end)
        self.delete(start, end)
        if insert_pos > start:
            insert_pos -= seg_len
        self.insert(insert_pos, seg)
        logger.info("translocate([%d:%d]→%d) len=%d", start, end, insert_pos, seg_len)

    # ------------------------------------------------------------------ #
    # Internals: locate, bounds, Fenwick, rebalancing
    # ------------------------------------------------------------------ #
    def _bounds_check_closed(self, start: int, end: int) -> None:
        n = len(self)
        if start < 0 or end < 0 or start > n or end > n:
            raise IndexError(f"slice bounds out of range: [0,{n}] got [{start},{end}]")

    def _locate_strict(self, pos: int) -> tuple[int, int]:
        n = len(self)
        if pos < 0 or pos >= n:
            raise IndexError(f"index out of range: 0 <= pos < {n}, got {pos}")
        return self._locate(pos)

    def _locate(self, pos: int) -> tuple[int, int]:
        """
        Map absolute index 0..len to (chunk_index, offset) where offset may equal
        chunk length when pos == len and we are at the end.
        """
        n = len(self)
        if pos < 0 or pos > n:
            raise IndexError(f"index out of range: 0 <= pos <= {n}, got {pos}")
        if not self._chunks:
            return (0, 0)

        if pos == n:
            ci = len(self._chunks) - 1
            return (ci, len(self._chunks[ci]))

        # Fenwick lower_bound: first index with prefix_sum > pos
        ci = self._bit_lower_bound(pos)
        prefix_before = self._prefix_sum(ci)
        return (ci, pos - prefix_before)

    # ----- Fenwick (BIT) helpers over chunk lengths ----- #
    def _rebuild_bit(self) -> None:
        size = len(self._chunks)
        self._bit = [0] * (size + 1)  # 1-based
        for i, chunk in enumerate(self._chunks, start=1):
            val = len(chunk)
            j = i
            while j <= size:
                self._bit[j] += val
                j += j & -j

    def _prefix_sum(self, i: int) -> int:
        """Sum of first i chunks (0..i-1)."""
        s = 0
        while i > 0:
            s += self._bit[i]
            i -= i & -i
        return s

    def _bit_lower_bound(self, pos: int) -> int:
        """
        Return the 0-based chunk index ci such that:
            prefix_sum(ci) <= pos < prefix_sum(ci + 1)
        Assumes 0 <= pos < total_length (pos == len handled by caller).
        """
        bit = self._bit
        n = len(bit) - 1  # number of chunks, Fenwick is 1-based
        if n <= 0:
            return 0

        idx = 0    # 1-based Fenwick index of the last prefix <= pos
        acc = 0    # accumulated sum at 'idx'
        step = 1 << (n.bit_length() - 1)

        # Binary lifting over the Fenwick tree
        while step:
            nxt = idx + step
            if nxt <= n and acc + bit[nxt] <= pos:
                idx = nxt
                acc += bit[nxt]
            step >>= 1

        # Convert 1-based Fenwick index to 0-based chunk index
        # idx is the count of chunks whose total length <= pos
        ci = min(idx, n - 1)
        return ci


    # ----- rebalancing ----- #
    def _rebalance_around(self, ci: int) -> None:
        """Split overly large chunks and merge tiny neighbors near ci."""
        if not self._chunks:
            return
        ci = max(0, min(ci, len(self._chunks) - 1))

        # Split while too large
        while len(self._chunks[ci]) > self._max_chunk:
            chunk = self._chunks[ci]
            mid = len(chunk) // 2
            left, right = chunk[:mid], chunk[mid:]
            self._chunks[ci:ci + 1] = [left, right]
            self._rebuild_bit()
            # continue rebalancing on the left piece
            if ci + 1 < len(self._chunks) and len(self._chunks[ci + 1]) > self._max_chunk:
                ci = ci + 1

        # Merge with previous/next if both are small
        def small(i: int) -> bool:
            return 0 <= i < len(self._chunks) and len(self._chunks[i]) < self._min_chunk

        merged = True
        while merged and len(self._chunks) > 1:
            merged = False
            if small(ci) and small(ci - 1):
                self._chunks[ci - 1].extend(self._chunks[ci])
                del self._chunks[ci]
                ci -= 1
                merged = True
            elif small(ci) and small(ci + 1):
                self._chunks[ci].extend(self._chunks[ci + 1])
                del self._chunks[ci + 1]
                merged = True
        self._rebuild_bit()
