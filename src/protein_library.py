# protein_library.py: Library for Protein-specific operations
import random

class ProteinUtils:

    @staticmethod
    def hydrophobicity(protein_sequence):
        hydrophobicity_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }

        score = sum(hydrophobicity_scale.get(aa, 0) for aa in protein_sequence)
        return score / len(protein_sequence) if protein_sequence else 0

    @staticmethod
    def validate(sequence):
        valid_bases = set("ACDEFGHIKLMNPQRSTVWY")
        return set(sequence.upper()).issubset(valid_bases)