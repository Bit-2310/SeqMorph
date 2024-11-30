from collections import Counter
from gene_library import SequenceAnalysis as SeqA

class SequenceAnalysisReport:
    """
    Class to calculate and compare sequence metrics for both normal and mutated sequences,
    including protein-specific analysis.
    """
    def __init__(self, normal_sequence, mutated_sequence, sequence_type="DNA"):
        self.normal_sequence = normal_sequence
        self.mutated_sequence = mutated_sequence
        self.sequence_type = sequence_type

    def gc_content(self, sequence):
        """Calculate GC content of the sequence."""
        return SeqA.gc_content(sequence)

    def nucleotide_frequency(self, sequence):
        """Calculate nucleotide frequency."""
        return SeqA.nucleotide_frequency(sequence)

    def kmer_frequency(self, sequence, k=3):
        """Calculate k-mer frequency."""
        return SeqA.kmer_count(sequence, k)

    def orf_detection(self, sequence):
        """Detect Open Reading Frames (ORFs)."""
        return SeqA.orf_detection(sequence)

    def compare_sequences(self):
        """Compare the metrics of normal vs mutated sequences and return a report."""
        normal_gc = self.gc_content(self.normal_sequence)
        mutated_gc = self.gc_content(self.mutated_sequence)

        normal_nucleotide_freq = self.nucleotide_frequency(self.normal_sequence)
        mutated_nucleotide_freq = self.nucleotide_frequency(self.mutated_sequence)

        normal_kmers = self.kmer_frequency(self.normal_sequence, k=3)
        mutated_kmers = self.kmer_frequency(self.mutated_sequence, k=3)

        normal_orfs = self.orf_detection(self.normal_sequence)
        mutated_orfs = self.orf_detection(self.mutated_sequence)

        # Generate report comparing each metric
        report = {
            "GC Content": {"Normal": normal_gc, "Mutated": mutated_gc},
            "Nucleotide Frequency": {
                "Normal": normal_nucleotide_freq,
                "Mutated": mutated_nucleotide_freq
            },
            "K-mer Frequency": {"Normal": normal_kmers, "Mutated": mutated_kmers},
            "Open Reading Frames": {"Normal": normal_orfs, "Mutated": mutated_orfs}
        }

        return report

    def plot_comparison(self, report):
        """Generate visual comparison of normal vs mutated sequence metrics."""
        import matplotlib.pyplot as plt

        # GC Content Comparison
        plt.bar(['Normal', 'Mutated'], [report["GC Content"]["Normal"], report["GC Content"]["Mutated"]])
        plt.title('GC Content Comparison')
        plt.ylabel('% GC Content')
        plt.show()

        # Nucleotide Frequency Comparison (bar chart for A, T, G, C)
        normal_freq = report["Nucleotide Frequency"]["Normal"]
        mutated_freq = report["Nucleotide Frequency"]["Mutated"]
        labels = ['A', 'T', 'G', 'C']
        normal_values = [normal_freq.get(base, 0) for base in labels]
        mutated_values = [mutated_freq.get(base, 0) for base in labels]

        x = range(len(labels))
        width = 0.35
        fig, ax = plt.subplots()
        ax.bar(x, normal_values, width, label='Normal')
        ax.bar([p + width for p in x], mutated_values, width, label='Mutated')

        ax.set_ylabel('Frequency')
        ax.set_title('Nucleotide Frequency Comparison')
        ax.set_xticks([p + width / 2 for p in x])
        ax.set_xticklabels(labels)
        ax.legend()

        plt.show()

        # K-mer Frequency Comparison
        normal_kmers = report["K-mer Frequency"]["Normal"]
        mutated_kmers = report["K-mer Frequency"]["Mutated"]

        kmers = list(set(normal_kmers).union(mutated_kmers))  # Combine unique k-mers
        normal_kmer_counts = [normal_kmers.get(kmer, 0) for kmer in kmers]
        mutated_kmer_counts = [mutated_kmers.get(kmer, 0) for kmer in kmers]

        plt.bar(kmers, normal_kmer_counts, width=0.4, label='Normal', align='center')
        plt.bar(kmers, mutated_kmer_counts, width=0.4, label='Mutated', align='edge')
        plt.title('K-mer Frequency Comparison')
        plt.ylabel('Frequency')
        plt.xticks(rotation=90)
        plt.legend()

        plt.show()
