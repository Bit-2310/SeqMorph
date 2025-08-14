from __future__ import annotations
import re
from typing import ClassVar, Literal

class SequenceValidation:
    """Validation and type detection for DNA, RNA and Protein sequences."""
    DNA_REGEX: ClassVar[re.Pattern[str]] = re.compile(r"^[ATGCRYSWKMBDHVN]+$", re.IGNORECASE)
    RNA_REGEX: ClassVar[re.Pattern[str]] = re.compile(r"^[AUGCRYSWKMBDHVN]+$", re.IGNORECASE)
    PROTEIN_REGEX: ClassVar[re.Pattern[str]] = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$", re.IGNORECASE)

    DNA_BASES: ClassVar[set[str]] = {"A", "T", "G", "C"}
    RNA_BASES: ClassVar[set[str]] = {"A", "U", "G", "C"}
    PROTEIN_BASES: ClassVar[set[str]] = {
        "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"
    }

    @staticmethod
    def validate(sequence: str, seq_type: Literal["DNA","RNA","Protein"]="DNA") -> bool:
        s = (sequence or "").strip().replace("\n", "").upper()
        if seq_type == "DNA":     return bool(SequenceValidation.DNA_REGEX.fullmatch(s))
        if seq_type == "RNA":     return bool(SequenceValidation.RNA_REGEX.fullmatch(s))
        if seq_type == "Protein": return bool(SequenceValidation.PROTEIN_REGEX.fullmatch(s))
        raise ValueError(f"Unknown sequence type: {seq_type}")

    @staticmethod
    def detect_seq(sequence: str) -> Literal["DNA","RNA","Protein"]:
        s = (sequence or "").strip().replace("\n", "").upper()
        if not s: raise ValueError("Empty sequence")
        if SequenceValidation.DNA_REGEX.fullmatch(s):     return "DNA"
        if SequenceValidation.RNA_REGEX.fullmatch(s):     return "RNA"
        if SequenceValidation.PROTEIN_REGEX.fullmatch(s): return "Protein"
        raise ValueError("Cannot determine sequence type")
