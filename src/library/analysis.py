from __future__ import annotations
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Literal, Tuple
from math import log2
from .validation import SequenceValidation
from .utils import GeneUtils

class SequenceAnalysis:
    """Nucleotide and AA level analyses used by CLI/GUI/API."""

    @staticmethod
    def gc_content(sequence: str) -> float:
        s = (sequence or "").upper()
        return 0.0 if not s else 100.0 * (s.count("G")+s.count("C"))/len(s)

    @staticmethod
    def nucleotide_frequency(sequence: str) -> Dict[str,int]:
        c = Counter((sequence or "").upper())
        return {k:int(c.get(k,0)) for k in ("A","C","G","T","U")}

    @staticmethod
    def kmer_count(sequence: str, k: int = 3, *, alphabet: Optional[Iterable[str]]=None) -> Dict[str,int]:
        if k <= 0: return {}
        s = (sequence or "").upper()
        if len(s) < k: raise ValueError("Sequence shorter than k-mer size")
        out: Dict[str,int] = defaultdict(int)
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            if alphabet and any(ch not in alphabet for ch in kmer): continue
            out[kmer] += 1
        return dict(out)

    @staticmethod
    def shannon_entropy(sequence: str, *, seq_type: Optional[Literal["DNA","RNA","Protein"]]=None) -> float:
        s = (sequence or "").upper()
        if not s: return 0.0
        stype = seq_type or SequenceValidation.detect_seq(s)
        if stype == "Protein":
            raise ValueError("Shannon entropy expects a nucleotide sequence.")
        counts = Counter(ch for ch in s if ch in "ATGCU")
        total = sum(counts.values())
        if total == 0: return 0.0
        return -sum((c/total) * log2(c/total) for c in counts.values())

    @staticmethod
    def translate(sequence: str, frame: int = 0, *, to_stop: bool=False) -> str:
        s = (sequence or "").upper()
        table = GeneUtils.get_standard_codon_table()
        if frame not in (0,1,2): raise ValueError("frame must be 0, 1 or 2")
        i, aa = frame, []
        while i + 2 < len(s):
            a = table.get(s[i:i+3], "X")
            if to_stop and a == "*": break
            aa.append(a); i += 3
        return "".join(aa)

    @staticmethod
    def orf_detection(sequence: str, seq_type: Literal["DNA","RNA"]="DNA",
                      min_length: int = 30, both_strands: bool=True) -> List[Dict]:
        s = (sequence or "").upper()
        if seq_type == "RNA":
            start_set, stop_set = {"AUG"}, {"UAA","UAG","UGA"}
            fwd, rev = s, GeneUtils.reverse_complement(s, "RNA")
        else:
            start_set, stop_set = {"ATG"}, {"TAA","TAG","TGA"}
            fwd, rev = s, GeneUtils.reverse_complement(s, "DNA")
        strands: List[Tuple[str,str]] = [("+", fwd)] + ([("-", rev)] if both_strands else [])
        out: List[Dict] = []
        for label, seq in strands:
            for frame in (0,1,2):
                i = frame
                while i + 2 < len(seq):
                    if seq[i:i+3] in start_set:
                        j = i + 3
                        while j + 2 < len(seq):
                            if seq[j:j+3] in stop_set:
                                nt_len = j + 3 - i
                                if nt_len >= min_length:
                                    if label == "+": start, end = i, j+3
                                    else: start, end = len(s)-(j+3), len(s)-i
                                    out.append({
                                        "strand": label, "frame": frame,
                                        "start": start, "end": end, "nt_length": nt_len,
                                        "protein": SequenceAnalysis.translate(
                                            s[start:end] if label=="+" else GeneUtils.reverse_complement(s[start:end])
                                        )
                                    })
                                break
                            j += 3
                    i += 3
        return out
