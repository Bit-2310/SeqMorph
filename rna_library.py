# rna_library.py: Library for RNA-specific operations

class RNAUtils:

    @staticmethod
    def transcribe(dna_sequence):
        """
        Transcribe a DNA sequence into an RNA sequence.
        """
        return dna_sequence.upper().replace("T", "U")

    @staticmethod
    def reverse_transcribe(rna_sequence):
        """
        Reverse transcribe an RNA sequence into a DNA sequence.
        """
        return rna_sequence.upper().replace("U", "T")

    @staticmethod
    def gc_content(sequence):
        """
        Calculate the GC content of a DNA or RNA sequence.
        """
        gc_count = sum(1 for base in sequence.upper() if base in "GC")
        return (gc_count / len(sequence)) * 100 if sequence else 0

    @staticmethod
    def is_sense_strand(sequence, reference_mrna):
        """
        Determine if the given RNA sequence matches the sense strand of the reference mRNA.
        """
        return sequence == reference_mrna

    @staticmethod
    def is_antisense_strand(sequence, reference_mrna):
        """
        Determine if the given RNA sequence is the antisense strand of the reference mRNA.
        """
        complement_map = str.maketrans("AUGC", "UACG")
        antisense = reference_mrna.translate(complement_map)[::-1]
        return sequence == antisense

    @staticmethod
    def codon_splitter(rna_sequence):
        """
        Split an RNA sequence into codons (groups of three bases).
        """
        if len(rna_sequence) % 3 != 0:
            raise ValueError("RNA sequence length is not a multiple of 3, incomplete codons found.")
        return [rna_sequence[i:i + 3] for i in range(0, len(rna_sequence), 3)]

    @staticmethod
    def translate(rna_sequence):
        """
        Translate an RNA sequence into a protein sequence.
        """
        codon_table = {
            'AUG': 'M',
            'UUU': 'F', 'UUC': 'F',
            'UUA': 'L', 'UUG': 'L',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'UAU': 'Y', 'UAC': 'Y',
            'CAU': 'H', 'CAC': 'H',
            'CAA': 'Q', 'CAG': 'Q',
            'AAU': 'N', 'AAC': 'N',
            'AAA': 'K', 'AAG': 'K',
            'GAU': 'D', 'GAC': 'D',
            'GAA': 'E', 'GAG': 'E',
            'UGU': 'C', 'UGC': 'C',
            'UGG': 'W',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'
        }

        protein = []
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i + 3]
            amino_acid = codon_table.get(codon, 'X')
            if amino_acid == 'STOP':
                break
            protein.append(amino_acid)
        return ''.join(protein)

if __name__ == "__main__":
    dna_sequence = "ATGCGTACGTTAGC"
    print("Original DNA Sequence:", dna_sequence)
    print("Transcribed RNA Sequence:", RNAUtils.transcribe(dna_sequence))

    rna_sequence = "AUGCGUACAGUUCAGCACGUAGCUAGAGCUAGCUAGCUAGCGUUAGC"
    print("Original RNA Sequence:", rna_sequence)
    print("Reverse Transcribed DNA Sequence:", RNAUtils.reverse_transcribe(rna_sequence))
    print("GC Content (%):", RNAUtils.gc_content(rna_sequence))
   # print("Codons:", RNAUtils.codon_splitter(rna_sequence))
    print("Translated Protein:", RNAUtils.translate(rna_sequence))
