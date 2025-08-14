# src/sequence_structure.py
from __future__ import annotations

import hashlib
import os
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from library import SequenceValidation
from sequence_store import BaseStore, StringStore


class SequenceStructure:
    """
    Registry for biological sequences with pluggable storage backends.

    - Stores sequences keyed by accession ID.
    - Persists per-sequence metadata (type, length, sha1, free-form fields).
    - Provides simple I/O helpers (FASTA in/out; plain text in).
    - Delegates sequence storage/editing to a BaseStore implementation.

    Coordinates: 0-based, half-open [start, end).
    """

    def __init__(
        self,
        on_duplicate: str = "error",
        store_factory: Optional[Callable[[str], BaseStore]] = None,
    ):
        """
        Parameters
        ----------
        on_duplicate : {"error","overwrite","suffix"}
            Behavior when adding an accession that already exists:
              - "error"      -> raise ValueError (default)
              - "overwrite"  -> replace existing sequence + metadata
              - "suffix"     -> auto-suffix _1, _2, ... until unique
        store_factory : Callable[[str], BaseStore], optional
            Factory to wrap raw strings into a storage backend. Defaults to StringStore.
        """
        if on_duplicate not in {"error", "overwrite", "suffix"}:
            raise ValueError('on_duplicate must be "error", "overwrite", or "suffix"')

        self.sequences: Dict[str, BaseStore] = {}
        self.metadata: Dict[str, Dict[str, Any]] = {}
        self._on_duplicate = on_duplicate
        self._store_factory = store_factory or (lambda s: StringStore(s))

    # ---------------------------------------------------------------------
    # Core add/update/remove
    # ---------------------------------------------------------------------
    def add_sequence(
        self,
        accession_id: str,
        sequence: str,
        seq_type: Optional[str] = None,
        *,
        overwrite: bool = False,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Add a sequence to the registry.

        Returns
        -------
        str
            The accession ID used (may be suffixed if on_duplicate == "suffix").
        """
        if not accession_id:
            raise ValueError("accession_id must be a non-empty string")

        seq = (sequence or "").strip()
        if not seq:
            raise ValueError("sequence is empty")

        # Auto-detect type if not provided; validate if provided.
        if seq_type is None:
            seq_type = SequenceValidation.detect_seq(seq)
        else:
            if not SequenceValidation.validate(seq, seq_type=seq_type):  # type: ignore[arg-type]
                raise ValueError(f"Invalid {seq_type} sequence")

        acc = self._resolve_duplicate_id(accession_id, allow_overwrite=overwrite)
        self.sequences[acc] = self._store_factory(seq)
        self.metadata[acc] = {
            "type": seq_type,
            "length": len(seq),
            "sha1": hashlib.sha1(seq.encode("utf-8")).hexdigest(),
            **(metadata or {}),
        }
        return acc

    def update_sequence(
        self,
        accession_id: str,
        sequence: str,
        seq_type: Optional[str] = None,
        *,
        metadata_updates: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Replace an existing sequence payload and refresh metadata."""
        self._require_exists(accession_id)

        seq = (sequence or "").strip()
        if not seq:
            raise ValueError("sequence is empty")

        if seq_type is None:
            seq_type = SequenceValidation.detect_seq(seq)
        else:
            if not SequenceValidation.validate(seq, seq_type=seq_type):  # type: ignore[arg-type]
                raise ValueError(f"Invalid {seq_type} sequence")

        self.sequences[accession_id] = self._store_factory(seq)
        self.metadata[accession_id].update(
            {
                "type": seq_type,
                "length": len(seq),
                "sha1": hashlib.sha1(seq.encode("utf-8")).hexdigest(),
            }
        )
        if metadata_updates:
            self.metadata[accession_id].update(metadata_updates)

    def remove_sequence(self, accession_id: str) -> None:
        """Remove a sequence and its metadata."""
        self._require_exists(accession_id)
        del self.sequences[accession_id]
        del self.metadata[accession_id]

    # ---------------------------------------------------------------------
    # Accessors
    # ---------------------------------------------------------------------
    def get_sequence(self, accession_id: str) -> str:
        """Return the stored sequence as a string."""
        self._require_exists(accession_id)
        return self.sequences[accession_id].to_string()

    def get_metadata(self, accession_id: str) -> Dict[str, Any]:
        """Get metadata dict for an accession (do not mutate in place)."""
        self._require_exists(accession_id)
        return dict(self.metadata[accession_id])

    def set_metadata(self, accession_id: str, **fields: Any) -> None:
        """Update metadata fields for an accession."""
        self._require_exists(accession_id)
        self.metadata[accession_id].update(fields)

    def list_accessions(self) -> List[str]:
        """List all accession IDs (in insertion order)."""
        return list(self.sequences.keys())

    # ---------------------------------------------------------------------
    # I/O
    # ---------------------------------------------------------------------
    def load_from_file(
        self,
        filepath: str,
        seq_type: Optional[str] = None,
        id_prefix: Optional[str] = None,
    ) -> int:
        """
        Load sequences from a file.

        - .fasta / .fa : parse all records via Biopython (lazy import).
        - otherwise    : read entire file as a single sequence (plain text).

        Returns
        -------
        int
            Number of sequences loaded.
        """
        filepath = os.fspath(filepath)
        ext = os.path.splitext(filepath)[1].lower()

        if ext in (".fasta", ".fa"):
            # Lazy import: Biopython is only required when FASTA I/O is used
            from Bio import SeqIO  # type: ignore

            count = 0
            for record in SeqIO.parse(filepath, "fasta"):
                acc = f"{id_prefix or ''}{record.id}"
                self.add_sequence(acc, str(record.seq), seq_type, overwrite=False)
                count += 1
            if count == 0:
                raise ValueError("No sequences found in the FASTA file.")
            return count

        # Plain text: entire file is a single sequence
        with open(filepath, "r", encoding="utf-8") as f:
            seq = f.read().strip()
        if not seq:
            raise ValueError("Plain text file contained an empty sequence")
        acc = f"{id_prefix or ''}{os.path.basename(filepath)}"
        self.add_sequence(acc, seq, seq_type, overwrite=False)
        return 1

    def export_to_file(
        self,
        filepath: str,
        accession_ids: Optional[List[str]] = None,
        wrap: int = 60,
    ) -> int:
        """
        Export sequences to a FASTA-like file.

        Parameters
        ----------
        filepath : str
            Output path.
        accession_ids : list[str] | None
            If provided, export only these accessions (in order). Otherwise export all.
        wrap : int
            Line wrap width; set 0 or negative to disable wrapping.

        Returns
        -------
        int
            Number of sequences written.
        """
        ids = accession_ids or self.list_accessions()
        if not ids:
            return 0

        with open(filepath, "w", encoding="utf-8") as fh:
            for acc in ids:
                self._require_exists(acc)
                seq = self.sequences[acc].to_string()
                meta = self.metadata[acc]
                desc = meta.get("description", "")
                fh.write(f">{acc} {desc}\n")
                if wrap and wrap > 0:
                    for i in range(0, len(seq), wrap):
                        fh.write(seq[i : i + wrap] + "\n")
                else:
                    fh.write(seq + "\n")
        return len(ids)

    # ---------------------------------------------------------------------
    # Serialization / iteration / stats
    # ---------------------------------------------------------------------
    def to_json(self) -> Dict[str, Dict]:
        """Return a serializable snapshot of all sequences + metadata."""
        return {
            "sequences": {acc: store.to_string() for acc, store in self.sequences.items()},
            "metadata": self.metadata,
        }

    def __iter__(self) -> Iterable[Tuple[str, str, Dict[str, Any]]]:
        """Iterate (accession_id, sequence_string, metadata_dict)."""
        for acc_id, store in self.sequences.items():
            yield acc_id, store.to_string(), self.metadata[acc_id]

    def get_full_sequence_length(self, accession_id: Optional[str] = None) -> int:
        """
        If accession_id is provided, return its length.
        Otherwise, return the sum of lengths of all sequences.
        """
        if accession_id is not None:
            self._require_exists(accession_id)
            return len(self.sequences[accession_id])  # BaseStore.__len__
        return sum(len(store) for store in self.sequences.values())

    # ---------------------------------------------------------------------
    # Internals
    # ---------------------------------------------------------------------
    def _require_exists(self, accession_id: str) -> None:
        if accession_id not in self.sequences:
            raise KeyError(f"Unknown accession_id: {accession_id}")

    def _resolve_duplicate_id(self, accession_id: str, *, allow_overwrite: bool) -> str:
        """
        Apply the on_duplicate policy and return a usable accession ID.
        """
        exists = accession_id in self.sequences
        if not exists:
            return accession_id

        if allow_overwrite or self._on_duplicate == "overwrite":
            # Caller wants to overwrite or policy allows overwrite.
            return accession_id

        if self._on_duplicate == "error":
            raise ValueError(f"accession_id already exists: {accession_id}")

        # Suffix mode: accession, accession_1, accession_2, ...
        base = accession_id
        i = 1
        new_id = f"{base}_{i}"
        while new_id in self.sequences:
            i += 1
            new_id = f"{base}_{i}"
        return new_id
