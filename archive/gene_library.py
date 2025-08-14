import re
import random
from collections import defaultdict, Counter
from typing import Dict, List, Literal, Optional


# -----------------------------
# Validation & detection helpers
# -----------------------------
class SequenceValidation:
    # Yes, we allow IUPAC ambiguity for nucleotides; proteins are the 20 canonical AAs.
    DNA_REGEX = re.compile(r"^[ATGCRYSWKMBDHVN]+$", re.IGNORECASE)
    RNA_REGEX = re.compile(r"^[AUGCRYSWKMBDHVN]+$", re.IGNORECASE)
    PROTEIN_REGEX = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$", re.IGNORECASE)

    DNA_BASES = {"A", "T", "G", "C"}
    RNA_BASES = {"A", "U", "G", "C"}
    PROTEIN_BASES = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
        "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
    }

    @staticmethod
    def validate(sequence: str, seq_type: str = "DNA") -> bool:
        """Quick sanity check. We normalize to uppercase and strip whitespace because... humans."""
        sequence = sequence.strip().replace("\n", "").upper()
        if seq_type == "DNA":
            return bool(SequenceValidation.DNA_REGEX.fullmatch(sequence))
        elif seq_type == "RNA":
            return bool(SequenceValidation.RNA_REGEX.fullmatch(sequence))
        elif seq_type == "Protein":
            return bool(SequenceValidation.PROTEIN_REGEX.fullmatch(sequence))
        raise ValueError(f"Unknown sequence type: {seq_type}")

    @staticmethod
    def detect_seq(sequence: str) -> str:
        """Try to guess the type. This is best-effort: ambiguous stuff can still fool us."""
        sequence = sequence.strip().replace("\n", "").upper()
        if SequenceValidation.DNA_REGEX.fullmatch(sequence):
            return "DNA"
        elif SequenceValidation.RNA_REGEX.fullmatch(sequence):
            return "RNA"
        elif SequenceValidation.PROTEIN_REGEX.fullmatch(sequence):
            return "Protein"
        raise ValueError("Cannot determine sequence type")


# -----------------------------
# Analysis utilities
# -----------------------------
class SequenceAnalysis:
    @staticmethod
    def gc_content(sequence: str) -> float:
        """GC% for nucleotides. If you pass protein here, the number is meaningless (don’t)."""
        seq = sequence.upper()
        if not seq:
            return 0.0
        gc = seq.count("G") + seq.count("C")
        return (gc / len(seq)) * 100

    @staticmethod
    def nucleotide_frequency(sequence: str) -> Dict[str, int]:
        """
        Count bases A/T/G/C/U. For RNA, T will be ~0; for DNA, U will be ~0.
        Yes, we keep both to make UI tables simple.
        """
        seq = sequence.upper()
        return {base: seq.count(base) for base in "ATGCU"}

    @staticmethod
    def kmer_count(sequence: str, k: int = 3) -> Dict[str, int]:
        """
        Count k-mers using Counter (because speed + readability).
        We throw if the sequence is shorter than k — caller can decide how to handle that.
        """
        seq = sequence.upper()
        if len(seq) < k:
            raise ValueError("Sequence shorter than k-mer size")
        return dict(Counter(seq[i:i + k] for i in range(len(seq) - k + 1)))

    @staticmethod
    def amino_acid_frequency(protein_seq: str) -> Dict[str, int]:
        """Counts per amino acid. Simple, fast, friendly."""
        seq = protein_seq.upper()
        return {aa: seq.count(aa) for aa in SequenceValidation.PROTEIN_BASES}

    @staticmethod
    def codon_usage(sequence: str) -> Dict[str, int]:
        """
        Count 3-mer codons (frame 0). We don’t assume ORFs here — it’s raw frame-0 usage.
        UI can still use this as a rough fingerprint.
        """
        seq = sequence.upper()
        codon_count: Dict[str, int] = defaultdict(int)
        # Stop before trailing leftovers; only full codons count
        for i in range(0, (len(seq) // 3) * 3, 3):
            codon = seq[i:i + 3]
            codon_count[codon] += 1
        return dict(codon_count)

    @staticmethod
    def shannon_entropy(sequence: str, *, seq_type: Optional[str] = None) -> float:
        """
        Shannon entropy for nucleotide sequences. If you hand us a protein sequence,
        we raise — not because we’re mean, but because mixing apples and ribosomes makes charts lie.
        """
        from math import log2

        seq = sequence.upper()
        if not seq:
            return 0.0

        # If caller didn’t say, we guess. If it looks like protein, we bail.
        stype = (seq_type or SequenceValidation.detect_seq(seq))
        if stype == "Protein":
            raise ValueError("Shannon entropy here expects a nucleotide sequence, not protein.")

        freq = SequenceAnalysis.nucleotide_frequency(seq)
        total = sum(freq.values())
        if total == 0:
            return 0.0

        # Classic Shannon entropy (bits). We skip zeros because log2(0) is a drama queen.
        return -sum((c / total) * log2(c / total) for c in freq.values() if c)

    @staticmethod
    def orf_detection(sequence: str, seq_type: str = "DNA") -> List[Dict[str, object]]:
        """
        Find ORFs on both strands, all three frames.
        We return a list of dicts so the GUI can show frame/strand/coords nicely.

        NOTE: start/end are indexes relative to the scanned strand’s frame slice,
        not absolute genome coordinates. If you need absolute positions, shout and we’ll map them.
        """
        stype = seq_type.upper()
        if stype == "RNA":
            start_codons = {"AUG"}
            stop_codons = {"UAA", "UAG", "UGA"}
        else:  # DNA
            start_codons = {"ATG"}
            stop_codons = {"TAA", "TAG", "TGA"}

        fwd = sequence.upper()
        rev = GeneUtils.reverse_complement(fwd, stype)
        results: List[Dict[str, object]] = []

        def scan(seq: str, strand: Literal["+", "-"]) -> None:
            for frame in (0, 1, 2):
                frame_seq = seq[frame:]
                i = 0
                while i <= len(frame_seq) - 3:
                    codon = frame_seq[i:i + 3]
                    if codon in start_codons:
                        # Walk until we hit a stop codon or run out of road
                        j = i + 3
                        while j <= len(frame_seq) - 3:
                            stop = frame_seq[j:j + 3]
                            if stop in stop_codons:
                                orf_seq = frame_seq[i:j + 3]
                                results.append({
                                    "orf_seq": orf_seq,
                                    "frame": frame,
                                    "strand": strand,
                                    "start_idx": i,   # 0-based in the frame slice
                                    "end_idx": j + 3  # exclusive
                                })
                                i = j  # hop ahead to avoid redundant starts inside ORFs
                                break
                            j += 3
                    i += 3

        scan(fwd, "+")
        scan(rev, "-")
        return results

    @staticmethod
    def translate_sequence(sequence: str, seq_type: str = "DNA") -> str:
        """
        Translate nucleotide sequence to protein using the standard table.
        RNA is converted to DNA alphabet first because the lookup is DNA-based.
        """
        codon_table = GeneUtils.get_standard_codon_table()
        seq = sequence.upper()
        if seq_type.upper() == "RNA":
            seq = seq.replace("U", "T")

        protein = []
        # Only full codons; trailing leftovers become nothing, not ‘X’
        for i in range(0, (len(seq) // 3) * 3, 3):
            codon = seq[i:i + 3]
            protein.append(codon_table.get(codon, "X"))
        return "".join(protein)


# -----------------------------
# Low-level gene helpers
# -----------------------------
class GeneUtils:
    @staticmethod
    def transition(codon: str) -> str:
        """Flip a random position by transition (purine<->purine or pyrimidine<->pyrimidine)."""
        codon = codon.upper()
        if len(codon) != 3:
            raise ValueError("Codon must be 3 bases")
        pos = random.randint(0, 2)
        base = codon[pos]
        transition_map = {"A": "G", "G": "A", "C": "T", "T": "C", "U": "C"}  # U -> C for RNA
        return codon[:pos] + transition_map.get(base, base) + codon[pos + 1:]

    @staticmethod
    def transversion(codon: str) -> str:
        """Flip a random position by transversion (purine<->pyrimidine)."""
        codon = codon.upper()
        if len(codon) != 3:
            raise ValueError("Codon must be 3 bases")
        pos = random.randint(0, 2)
        base = codon[pos]
        transversion_map = {
            "A": ["C", "T"], "G": ["C", "T"],
            "C": ["A", "G"], "T": ["A", "G"], "U": ["A", "G"]
        }
        new_base = random.choice(transversion_map.get(base, [base]))
        return codon[:pos] + new_base + codon[pos + 1:]

    @staticmethod
    def is_valid_codon(codon: str, seq_type: str = "DNA") -> bool:
        """Three letters. Right alphabet. That’s the tweet."""
        cod = codon.strip().upper()
        if len(cod) != 3:
            return False
        valid_bases = SequenceValidation.DNA_BASES if seq_type.upper() == "DNA" else SequenceValidation.RNA_BASES
        return all(b in valid_bases for b in cod)

    @staticmethod
    def reverse_complement(sequence: str, seq_type: str = "DNA") -> str:
        """Reverse-complement for DNA or RNA. Ambiguity codes aren’t handled (yet)."""
        seq = sequence.upper()
        if seq_type.upper() == "DNA":
            return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
        elif seq_type.upper() == "RNA":
            return seq.translate(str.maketrans("AUGC", "UACG"))[::-1]
        raise ValueError("Invalid sequence type")

    @staticmethod
    def get_standard_codon_table() -> Dict[str, str]:
        """Standard translation table (DNA alphabet). Stars are stops."""
        return {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
