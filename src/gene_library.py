import random

class SequenceValidation:
    """
    Handles validation and type detection for DNA, RNA, and ambiguous sequences.
    """
    DNA_BASES = {"A", "T", "G", "C"}
    RNA_BASES = {"A", "U", "G", "C"}
    PROTEIN_BASES = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
        "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
    }
    AMBIGUITY_CODES = {
        "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"],
        "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
        "B": ["C", "G", "T"], "D": ["A", "G", "T"],
        "H": ["A", "C", "T"], "V": ["A", "C", "G"], "N": ["A", "C", "G", "T"]
    }

    @staticmethod
    def validate(sequence, seq_type="DNA"):
        """
        Validate a sequence based on the type (DNA, RNA, or Protein).
        :param sequence: The sequence to validate.
        :param seq_type: The type of the sequence ('DNA', 'RNA', or 'Protein').
        :return: True if the sequence is valid.
        """
        sequence = sequence.upper()
        if seq_type == "DNA":
            return set(sequence).issubset(SequenceValidation.DNA_BASES.union(SequenceValidation.AMBIGUITY_CODES))
        elif seq_type == "RNA":
            return set(sequence).issubset(SequenceValidation.RNA_BASES.union(SequenceValidation.AMBIGUITY_CODES))
        elif seq_type == "Protein":
            return set(sequence).issubset(SequenceValidation.PROTEIN_BASES)
        else:
            raise ValueError(f"Unknown sequence type: {seq_type}")

    @staticmethod
    def detect_seq(sequence):
        """
        Detect the type of a sequence (DNA, RNA, Protein).
        :param sequence: The sequence to detect.
        :return: The type of sequence ('DNA', 'RNA', or 'Protein').
        """
        sequence = sequence.upper()
        if set(sequence).issubset(SequenceValidation.DNA_BASES.union(SequenceValidation.AMBIGUITY_CODES)):
            return "DNA"
        elif set(sequence).issubset(SequenceValidation.RNA_BASES.union(SequenceValidation.AMBIGUITY_CODES)):
            return "RNA"
        elif set(sequence).issubset(SequenceValidation.PROTEIN_BASES):
            return "Protein"
        else:
            raise ValueError("Invalid sequence: Cannot determine if DNA, RNA, or Protein.")


class SequenceAnalysis:
    """
    Provides analytical methods for DNA, RNA, and Protein sequences.
    """
    @staticmethod
    def gc_content(sequence):
        gc_count = sum(1 for base in sequence.upper() if base in "GC")
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    @staticmethod
    def nucleotide_frequency(sequence):
        sequence = sequence.upper()
        valid_bases = {"A", "T", "G", "C", "U"}
        return {base: sequence.count(base) for base in valid_bases if base in sequence}
    
    @staticmethod
    def kmer_count(sequence, k=3):
        if len(sequence) < k:
            raise ValueError("Sequence length is shorter than k-mer size.")
        return {sequence[i:i + k]: sequence.count(sequence[i:i + k]) for i in range(len(sequence) - k + 1)}

    @staticmethod
    def orf_detection(sequence):
        """
        Detects Open Reading Frames (ORFs) in a DNA or RNA sequence.
        ORFs start with 'ATG' (start codon) and end with any stop codon.
        """
        start_codon = "ATG"
        stop_codons = {"TAA", "TAG", "TGA"}
        sequence = sequence.upper()

        orfs = []
        for i in range(len(sequence) - 2):
            if sequence[i:i + 3] == start_codon:
                for j in range(i + 3, len(sequence) - 2, 3):
                    if sequence[j:j + 3] in stop_codons:
                        orfs.append(sequence[i:j + 3])
                        break
        return orfs


class GeneUtils:
    """
    Manages mutation-specific logic for DNA and RNA sequences.
    """
    @staticmethod
    def transition(codon):
        """
        Performs a transition mutation (A ↔ G, C ↔ T/U).
        """
        transition_map = {"A": "G", "G": "A", "C": "T", "T": "C", "U": "C"}
        return "".join(transition_map.get(base, base) for base in codon)

    @staticmethod
    def transversion(codon):
        """
        Performs a transversion mutation (A ↔ C/T, G ↔ C/T).
        """
        transversion_map = {"A": ["C", "T"], "G": ["C", "T"], "C": ["A", "G"], "T": ["A", "G"], "U": ["A", "G"]}
        return "".join(random.choice(transversion_map.get(base, [base])) for base in codon)
    
    @staticmethod
    def is_valid_codon(codon, seq_type="DNA"):
        """
        Checks if the codon is valid based on the sequence type (DNA or RNA).
        """
        codon = codon.upper()
        if seq_type == "DNA":
            return codon in SequenceValidation.DNA_BASES
        elif seq_type == "RNA":
            return codon in SequenceValidation.RNA_BASES
        return False