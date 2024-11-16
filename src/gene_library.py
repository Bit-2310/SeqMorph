# gene_library.py : Library for DNA and RNA-specific operations
class GeneUtils:
    """
    A unified library for handling DNA and RNA sequences.
    """
    CODON_TABLE_DNA = {
                                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                                'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
                                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                                'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
                                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                        }
    STOP_CODONS_DNA = {'TAA', 'TAG', 'TGA'}

    START_CODON_RNA = 'AUG'
    CODON_TABLE_RNA = { codon.replace('T', 'U'): aa for codon, aa in CODON_TABLE_DNA.items() }
    STOP_CODONS_RNA = { codon.replace('T', 'U') for codon in STOP_CODONS_DNA }

    # Ambiguity codes and their possible matches
    AMBIGUITY_CODES = {
        "R": ["A", "G"],  # Purine
        "Y": ["C", "T"],  # Pyrimidine
        "S": ["G", "C"],  # Strong interaction
        "W": ["A", "T"],  # Weak interaction
        "K": ["G", "T"],  # Keto
        "M": ["A", "C"],  # Amino
        "B": ["C", "G", "T"],  # Not A
        "D": ["A", "G", "T"],  # Not C
        "H": ["A", "C", "T"],  # Not G
        "V": ["A", "C", "G"],  # Not T
        "N": ["A", "C", "G", "T"],  # Any base
    }

    TRANSITION_MAP = {"A": "G", "G": "A", "C": "T", "T": "C", "U": "C"}
    TRANSVERSION_MAP = {
        "A": ["C", "T"],
        "G": ["C", "T"],
        "C": ["A", "G"],
        "T": ["A", "G"],
        "U": ["A", "G"]
    }

    @staticmethod
    def detect_sequence_type(sequence):
        """
        Detect whether a sequence is DNA or RNA.
        Args:
            sequence (str): Input nucleotide sequence.
        Returns:
            str: "DNA" or "RNA".
        """
        dna_bases = set("ATGC")
        rna_bases = set("AUGC")
        sequence_set = set(sequence.upper())

        if sequence_set.issubset(dna_bases):
            return "DNA"
        elif sequence_set.issubset(rna_bases):
            return "RNA"
        else:
            raise ValueError("Invalid sequence: Cannot determine if DNA or RNA.")

    @staticmethod
    def complement(sequence):
        """
        Generate the complement of a DNA or RNA sequence.
        Args:
            sequence (str): Input sequence.
        Returns:
            str: Complemented sequence.
        """
        dna_complement_map = str.maketrans("ATGC", "TACG")
        rna_complement_map = str.maketrans("AUGC", "UACG")

        if "T" in sequence.upper():
            return sequence.upper().translate(dna_complement_map)
        elif "U" in sequence.upper():
            return sequence.upper().translate(rna_complement_map)
        else:
            raise ValueError("Invalid sequence: Cannot determine if DNA or RNA.")

    @staticmethod
    def reverse_complement(sequence):
        return GeneUtils.complement(sequence)[::-1]

    @staticmethod
    def gc_content(sequence):
        """
        Calculate the GC content of a DNA or RNA sequence.
        Args:
            sequence (str): Input sequence.
        Returns:
            float: GC content as a percentage.
        """
        gc_count = sum(1 for base in sequence.upper() if base in "GC")
        return (gc_count / len(sequence)) * 100 if sequence else 0

    @staticmethod
    def nucleotide_frequency(sequence):
        """
        Calculate the nucleotide frequency of a DNA or RNA sequence.
        Args:
            sequence (str): Input sequence.
        Returns:
            dict: Dictionary of nucleotide frequencies.
        """
        sequence = sequence.upper()
        valid_bases = {"A", "T", "G", "C", "U"}
        if not set(sequence).issubset(valid_bases):
            raise ValueError("Invalid sequence: Contains non-nucleotide characters.")
        return {base: sequence.count(base) for base in valid_bases if base in sequence}

    @staticmethod
    def transcribe(dna_sequence):
        return dna_sequence.upper().replace("T", "U")

    @staticmethod
    def reverse_transcribe(rna_sequence):
        return rna_sequence.upper().replace("U", "T")

    @staticmethod
    def codon_splitter(sequence):
        """
        Split a DNA or RNA sequence into codons (groups of three bases).
        Args:
            sequence (str): Input sequence.
        Returns:
            list[str]: List of codons.
        """
        if len(sequence) % 3 != 0:
            raise ValueError("Sequence length is not a multiple of 3, incomplete codons found.")
        return [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

    @staticmethod
    def transition_counts(sequence):
        """
        Count the number of transitions (A ↔ G, C ↔ T/U) in a DNA or RNA sequence.
        Args:
            sequence (str): Input sequence.
        Returns:
            int: Number of transitions.
        """
        transition_count = 0 
        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C'), ('C', 'U'), ('U', 'C')}
        for i in range(len(sequence) - 1):
            pair = (sequence[i].upper(), sequence[i + 1].upper())
            if pair in transitions:
                transition_count += 1

        return transition_count

    @staticmethod
    def transversion_counts(sequence):
        """
        Count the number of transversions (A ↔ C/T/U, G ↔ C/T/U) in a DNA or RNA sequence.
        Args:
            sequence (str): Input sequence.
        Returns:
            int: Number of transversions.
        """
        transversion_counts = 0
        purines = {'A', 'G'}
        pyrimidines = {'C', 'T', 'U'}
        for i in range(len(sequence) - 1):
            first = sequence[i].upper()
            second = sequence[i + 1].upper()

            if (first in purines and second in pyrimidines) or (first in pyrimidines and second in purines):
                transversion_count += 1

        return transversion_count

    @staticmethod
    def validate(sequence):
        valid_bases = {"A", "T", "G", "C", "U"}.union(GeneUtils.AMBIGUITY_CODES.keys())
        return set(sequence.upper()).issubset(valid_bases)
    
    def resolve_ambiguity(base):
        return GeneUtils.AMBIGUITY_CODES.get(base.upper(), [base.upper()])
    
    @staticmethod
    def count_ambiguities(sequence):
        return sum(1 for base in sequence if base.upper() in GeneUtils.AMBIGUITY_CODES)
    
    @staticmethod
    def is_transition(codon1, codon2):
        """
        Check if a mutation from codon1 to codon2 is a transition.
        :param codon1: Original codon.
        :param codon2: Mutated codon.
        :return: True if transition, False otherwise.
        """
        purines = {'A', 'G'}
        pyrimidines = {'C', 'T'}
        transitions = {
            ('A', 'G'), ('G', 'A'),
            ('C', 'T'), ('T', 'C')
        }

        for b1, b2 in zip(codon1, codon2):
            if (b1, b2) in transitions:
                return True
        return False

    @staticmethod
    def is_transversion(codon1, codon2):
        """
        Check if a mutation from codon1 to codon2 is a transversion.
        :param codon1: Original codon.
        :param codon2: Mutated codon.
        :return: True if transversion, False otherwise.
        """
        purines = {'A', 'G'}
        pyrimidines = {'C', 'T'}
        for b1, b2 in zip(codon1, codon2):
            if (b1 in purines and b2 in pyrimidines) or (b1 in pyrimidines and b2 in purines):
                return True
        return False
