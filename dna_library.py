# dna_library.py: Library for DNA-specific operations

class DNAUtils:
    @staticmethod
    def complement(sequence):
        complement_map = str.maketrans("ATGC", "TACG")
        return sequence.upper().translate(complement_map)

    @staticmethod
    def reverse_complement(sequence):
        return DNAUtils.complement(sequence)[::-1]

    @staticmethod
    def gc_content(sequence):
        gc_count = sum(1 for base in sequence.upper() if base in "GC")
        return (gc_count / len(sequence)) * 100

    @staticmethod
    def nucleotide_frequency(sequence):
        sequence = sequence.upper()
        if not set(sequence).issubset({"A", "T", "G", "C"}):
            raise ValueError("Invalid DNA sequence: Contains non-DNA characters.")
        return {base: sequence.count(base) for base in "ATGC"}

class DNAtoRNA:
    @staticmethod
    def rev_transcribe(sequence):
        return sequence.upper().replace("T", "U")

if __name__ == "__main__":
    sequence = "ATGCGTACGTTAGC"
    print("Original Sequence:", sequence)
    print("Complement:", DNAUtils.complement(sequence))
    print("Reverse Complement:", DNAUtils.reverse_complement(sequence))
    print("GC Content (%):", DNAUtils.gc_content(sequence))
    print("Nucleotide Frequency:", DNAUtils.nucleotide_frequency(sequence))
