import re
import random
from collections import defaultdict

class SequenceValidation:
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
        if seq_type == "DNA":
            return bool(SequenceValidation.DNA_REGEX.fullmatch(sequence))
        elif seq_type == "RNA":
            return bool(SequenceValidation.RNA_REGEX.fullmatch(sequence))
        elif seq_type == "Protein":
            return bool(SequenceValidation.PROTEIN_REGEX.fullmatch(sequence))
        raise ValueError(f"Unknown sequence type: {seq_type}")

    @staticmethod
    def detect_seq(sequence: str) -> str:
        if SequenceValidation.DNA_REGEX.fullmatch(sequence):
            return "DNA"
        elif SequenceValidation.RNA_REGEX.fullmatch(sequence):
            return "RNA"
        elif SequenceValidation.PROTEIN_REGEX.fullmatch(sequence):
            return "Protein"
        raise ValueError("Cannot determine sequence type")
class SequenceAnalysis:
    @staticmethod
    def gc_content(sequence: str) -> float:
        gc = sequence.upper().count("G") + sequence.upper().count("C")
        return (gc / len(sequence)) * 100 if sequence else 0.0

    @staticmethod
    def nucleotide_frequency(sequence: str) -> dict:
        return {base: sequence.upper().count(base) for base in "ATGCU"}

    @staticmethod
    def kmer_count(sequence: str, k: int = 3) -> dict:
        if len(sequence) < k:
            raise ValueError("Sequence shorter than k-mer size")
        kmers = defaultdict(int)
        for i in range(len(sequence) - k + 1):
            kmers[sequence[i:i+k].upper()] += 1
        return kmers

    @staticmethod
    def orf_detection(sequence: str, seq_type: str = "DNA") -> list:
        start_codons = {"ATG", "AUG"} if seq_type == "RNA" else {"ATG"}
        stop_codons = {"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"}

        rev_seq = GeneUtils.reverse_complement(sequence, seq_type)
        orfs = []

        for seq in [sequence, rev_seq]:
            for frame in [0, 1, 2]:
                frame_seq = seq[frame:]
                for i in range(0, len(frame_seq)-2, 3):
                    codon = frame_seq[i:i+3].upper()
                    if codon in start_codons:
                        for j in range(i+3, len(frame_seq)-2, 3):
                            stop = frame_seq[j:j+3].upper()
                            if stop in stop_codons:
                                orfs.append(frame_seq[i:j+3])
                                break
        return orfs

class GeneUtils:
    @staticmethod
    def transition(codon: str) -> str:
        if len(codon) != 3:
            raise ValueError("Codon must be 3 bases")
        pos = random.randint(0, 2)
        base = codon[pos]
        transition_map = {"A": "G", "G": "A", "C": "T", "T": "C", "U": "C"}
        return codon[:pos] + transition_map.get(base, base) + codon[pos+1:]

    @staticmethod
    def transversion(codon: str) -> str:
        if len(codon) != 3:
            raise ValueError("Codon must be 3 bases")
        pos = random.randint(0, 2)
        base = codon[pos]
        transversion_map = {
            "A": ["C", "T"], "G": ["C", "T"],
            "C": ["A", "G"], "T": ["A", "G"], "U": ["A", "G"]
        }
        new_base = random.choice(transversion_map.get(base, [base]))
        return codon[:pos] + new_base + codon[pos+1:]

    @staticmethod
    def is_valid_codon(codon: str, seq_type: str = "DNA") -> bool:
        valid_bases = SequenceValidation.DNA_BASES if seq_type == "DNA" else SequenceValidation.RNA_BASES
        return all(base.upper() in valid_bases for base in codon)

    @staticmethod
    def reverse_complement(sequence: str, seq_type: str = "DNA") -> str:
        if seq_type == "DNA":
            return sequence.translate(str.maketrans("ATGC", "TACG"))[::-1]
        elif seq_type == "RNA":
            return sequence.translate(str.maketrans("AUGC", "UACG"))[::-1]
        raise ValueError("Invalid sequence type")