import re
import random
from collections import defaultdict
from typing import Dict, List

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
        sequence = sequence.strip().replace("\n", "").upper()
        if SequenceValidation.DNA_REGEX.fullmatch(sequence):
            return "DNA"
        elif SequenceValidation.RNA_REGEX.fullmatch(sequence):
            return "RNA"
        elif SequenceValidation.PROTEIN_REGEX.fullmatch(sequence):
            return "Protein"
        elif "U" in sequence and "T" not in sequence:
            return "RNA"
        elif "T" in sequence and "U" not in sequence:
            return "DNA"
        elif set(sequence).issubset(SequenceValidation.PROTEIN_BASES):
            return "Protein"
        raise ValueError("Cannot determine sequence type")


class SequenceAnalysis:
    @staticmethod
    def gc_content(sequence: str, verbose: bool = False):
        gc = sequence.upper().count("G") + sequence.upper().count("C")
        percent = (gc / len(sequence)) * 100 if sequence else 0.0
        return (percent, f"{gc} G/C bases out of {len(sequence)}") if verbose else percent

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
    def amino_acid_frequency(protein_seq: str) -> dict:
        return {aa: protein_seq.upper().count(aa) for aa in SequenceValidation.PROTEIN_BASES}

    @staticmethod
    def amino_acid_composition(protein_seq: str) -> dict:
        freq = SequenceAnalysis.amino_acid_frequency(protein_seq)
        total = sum(freq.values())
        return {aa: (count / total) * 100 for aa, count in freq.items()} if total else {}

    @staticmethod
    def codon_usage(sequence: str) -> dict:
        codon_count = defaultdict(int)
        sequence = sequence.upper()
        usable_length = len(sequence) - (len(sequence) % 3)
        sequence = sequence[:usable_length]
        for i in range(0, usable_length, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                codon_count[codon] += 1
        return dict(codon_count)

    @staticmethod
    def shannon_entropy(sequence: str) -> float:
        from math import log2
        freq = SequenceAnalysis.nucleotide_frequency(sequence)
        total = sum(freq.values())
        if total == 0:
            return 0.0
        entropy = -sum((count/total) * log2(count/total) for count in freq.values() if count)
        return entropy

    @staticmethod
    def orf_detection(sequence: str, seq_type: str = "DNA", min_length: int = 30) -> list:
        start_codons = {"AUG"} if seq_type == "RNA" else {"ATG"}
        stop_codons = {"UAA", "UAG", "UGA", "TAA", "TAG", "TGA"}
        rev_seq = GeneUtils.reverse_complement(sequence, seq_type)
        orfs = []

        for strand, seq in [("forward", sequence), ("reverse", rev_seq)]:
            for frame in [0, 1, 2]:
                for i in range(frame, len(seq)-2, 3):
                    codon = seq[i:i+3].upper()
                    if codon in start_codons:
                        for j in range(i+3, len(seq)-2, 3):
                            stop = seq[j:j+3].upper()
                            if stop in stop_codons:
                                orf_seq = seq[i:j+3]
                                if len(orf_seq) >= min_length:
                                    orfs.append({
                                        "frame": frame,
                                        "strand": strand,
                                        "start": i,
                                        "end": j+3,
                                        "sequence": orf_seq
                                    })
                                break
        return orfs

    @staticmethod
    def translate_sequence(sequence: str, seq_type: str = "DNA") -> str:
        codon_table = GeneUtils.get_standard_codon_table()
        sequence = sequence.upper().replace("U", "T") if seq_type == "RNA" else sequence.upper()
        usable_length = len(sequence) - (len(sequence) % 3)
        sequence = sequence[:usable_length]
        protein = ""
        for i in range(0, usable_length, 3):
            codon = sequence[i:i+3]
            protein += codon_table.get(codon, "X")
        return protein

    @staticmethod
    def reverse_translate(protein_seq: str) -> str:
        AMINO_TO_CODONS = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'N': ['AAT', 'AAC'],
            'D': ['GAT', 'GAC'],
            'C': ['TGT', 'TGC'],
            'Q': ['CAA', 'CAG'],
            'E': ['GAA', 'GAG'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'K': ['AAA', 'AAG'],
            'M': ['ATG'],
            'F': ['TTT', 'TTC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            '*': ['TAA', 'TAG', 'TGA']
        }
        dna = ""
        for aa in protein_seq.upper():
            dna += random.choice(AMINO_TO_CODONS.get(aa, ['NNN']))
        return dna


class GeneUtils:
    @staticmethod
    def transition(codon: str) -> str:
        codon = codon.upper()
        if len(codon) != 3:
            raise ValueError("Codon must be 3 bases")
        pos = random.randint(0, 2)
        base = codon[pos]
        transition_map = {"A": "G", "G": "A", "C": "T", "T": "C", "U": "C"}
        return codon[:pos] + transition_map.get(base, base) + codon[pos+1:]

    @staticmethod
    def transversion(codon: str) -> str:
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
        return codon[:pos] + new_base + codon[pos+1:]

    @staticmethod
    def is_valid_codon(codon: str, seq_type: str = "DNA") -> bool:
        codon = codon.strip().upper()
        if len(codon) != 3:
            return False
        valid_bases = SequenceValidation.DNA_BASES if seq_type == "DNA" else SequenceValidation.RNA_BASES
        return all(base in valid_bases for base in codon)

    @staticmethod
    def reverse_complement(sequence: str, seq_type: str = "DNA") -> str:
        sequence = sequence.upper()
        if seq_type == "DNA":
            return sequence.translate(str.maketrans("ATGC", "TACG"))[::-1]
        elif seq_type == "RNA":
            return sequence.translate(str.maketrans("AUGC", "UACG"))[::-1]
        raise ValueError("Invalid sequence type")

    @staticmethod
    def get_standard_codon_table() -> Dict[str, str]:
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
