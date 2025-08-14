from __future__ import annotations
import random
from typing import Dict, Literal

class GeneUtils:
    """Low-level helpers used across modules."""

    @staticmethod
    def reverse_complement(sequence: str, seq_type: Literal["DNA","RNA"]="DNA") -> str:
        s = (sequence or "").upper()
        if seq_type == "DNA": return s.translate(str.maketrans("ATGC","TACG"))[::-1]
        if seq_type == "RNA": return s.translate(str.maketrans("AUGC","UACG"))[::-1]
        raise ValueError("Invalid seq_type for reverse_complement")

    @staticmethod
    def get_standard_codon_table() -> Dict[str,str]:
        return {
            'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
        }

    @staticmethod
    def transition_base(base: str) -> str:
        return {"A":"G","G":"A","C":"T","T":"C","U":"C"}.get(base.upper(), base)

    @staticmethod
    def transversion_base(base: str) -> str:
        opts = {"A":["C","T"],"G":["C","T"],"C":["A","G"],"T":["A","G"],"U":["A","G"]}.get(base.upper(), [base])
        return random.choice(opts)
