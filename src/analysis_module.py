import os
from collections import Counter
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, ttest_ind
from gene_library import SequenceAnalysis as SeqA
from Bio import Align, pairwise2
from scipy.stats import ttest_1samp
class SequenceAnalysisReport:
    """
    Analyzes and compares normal vs mutated sequences, generating reports and visualizations.
    Supports DNA, RNA, and protein sequences. Includes statistical tests for significance.
    """
    def __init__(self, normal_sequence: str, mutated_sequence: str, sequence_type: str = "DNA"):
        """
        Initialize with normal and mutated sequences.
        """
        if not normal_sequence or not mutated_sequence:
            raise ValueError("Sequences cannot be empty.")
        if sequence_type not in ["DNA", "RNA", "Protein"]:
            raise ValueError("Invalid sequence type. Must be DNA, RNA, or Protein.")

        self.normal_sequence = normal_sequence
        self.mutated_sequence = mutated_sequence
        self.sequence_type = sequence_type
        self._cached_results = {}  # Cache results for reuse

    def align_sequences(self) -> dict:
        """Performs global sequence alignment using Bio.Align.PairwiseAligner."""
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -1

        alignment = aligner.align(self.normal_sequence, self.mutated_sequence)[0]
        return {
            "Original": alignment.target,
            "Mutated": alignment.query,
            "Score": alignment.score
    }

    def compare_sequences(self) -> dict:
        """
        Compare normal vs mutated sequences and return a structured report.
        """
        if "report" in self._cached_results:
            return self._cached_results["report"]

        normal_orfs = SeqA.orf_detection(self.normal_sequence)
        mutated_orfs = SeqA.orf_detection(self.mutated_sequence)

        report = {
            "GC Content": {
                "Normal": SeqA.gc_content(self.normal_sequence),
                "Mutated": SeqA.gc_content(self.mutated_sequence)
            },
            "Nucleotide Frequency": {
                "Normal": SeqA.nucleotide_frequency(self.normal_sequence),
                "Mutated": SeqA.nucleotide_frequency(self.mutated_sequence)
            },
            "K-mer Frequency": {
                "Normal": SeqA.kmer_count(self.normal_sequence),
                "Mutated": SeqA.kmer_count(self.mutated_sequence)
            },
            "Alignment": self.align_sequences(),
            "Open Reading Frames": {
                "Count": {
                    "Normal": len(normal_orfs),
                    "Mutated": len(mutated_orfs)
                },
                "Sequences": {
                    "Normal": normal_orfs,
                    "Mutated": mutated_orfs
                }
            }
        }

        self._cached_results["report"] = report
        report["Statistical Tests"] = {
            "GC Content": self._compare_gc_content(report),
            "Nucleotide Frequency": self._compare_nucleotide_freq(report),
            "K-mer Frequency": self._compare_kmer_freq(report)
        }
        return report
    
    def _compare_gc_content(self, report: Dict) -> Dict:
        """Perform t-test for GC content differences."""
        normal_gc = report["GC Content"]["Normal"]
        mutated_gc = report["GC Content"]["Mutated"]
        if normal_gc == 0 or mutated_gc == 0:
            return {"t_stat": None, "p_value": 1.0}  # Prevent division by zero
        t_stat, p_value = ttest_1samp([mutated_gc], normal_gc)
        return {"t_stat": t_stat, "p_value": p_value}

    def _compare_nucleotide_freq(self, report: Dict) -> Dict:
        """Perform chi-square test for nucleotide frequency differences."""
        normal_freq = list(report["Nucleotide Frequency"]["Normal"].values())
        mutated_freq = list(report["Nucleotide Frequency"]["Mutated"].values())
        if any(val == 0 for val in normal_freq) or any(val == 0 for val in mutated_freq):
            return {"chi2_stat": None, "p_value": 1.0}  # Avoid division errors
        chi2, p, _, _ = chi2_contingency([normal_freq, mutated_freq])
        return {"chi2_stat": chi2, "p_value": p}

    def _compare_kmer_freq(self, report: Dict) -> Dict:
        """Perform chi-square test for k-mer frequency differences."""
        normal_kmers = report["K-mer Frequency"]["Normal"]
        mutated_kmers = report["K-mer Frequency"]["Mutated"]
        kmers = list(set(normal_kmers).union(mutated_kmers))
        normal_counts = [normal_kmers.get(kmer, 0) for kmer in kmers]
        mutated_counts = [mutated_kmers.get(kmer, 0) for kmer in kmers]
        chi2, p, _, _ = chi2_contingency([normal_counts, mutated_counts])
        return {"chi2_stat": chi2, "p_value": p}

    def plot_comparison(self, report: Optional[Dict] = None, save_path: Optional[str] = None):
        """
        Generate visual comparison of normal vs mutated sequence metrics.
        """
        if report is None:
            report = self.compare_sequences()

        if save_path:
            os.makedirs(os.path.dirname(save_path), exist_ok=True)

        self._plot_gc_content(report, save_path)
        self._plot_nucleotide_frequency(report, save_path)
        self._plot_kmer_frequency(report, save_path)

    def _plot_gc_content(self, report: Dict, save_path: Optional[str] = None):
        """Plot and save GC content comparison."""
        labels = ['Normal', 'Mutated']
        values = [report["GC Content"]["Normal"], report["GC Content"]["Mutated"]]

        plt.bar(labels, values)
        plt.title('GC Content Comparison')
        plt.ylabel('% GC Content')
        if save_path:
            plt.savefig(f"{save_path}_gc_content.png")
        plt.show()

    def _plot_nucleotide_frequency(self, report: Dict, save_path: Optional[str] = None):
        """Plot and save nucleotide frequency comparison."""
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

        if save_path:
            plt.savefig(f"{save_path}_nucleotide_frequency.png")
        plt.show()

    def _plot_kmer_frequency(self, report: Dict, save_path: Optional[str] = None):
        """Plot and save k-mer frequency comparison."""
        normal_kmers = report["K-mer Frequency"]["Normal"]
        mutated_kmers = report["K-mer Frequency"]["Mutated"]

        kmers = list(set(normal_kmers).union(mutated_kmers))
        normal_kmer_counts = [normal_kmers.get(kmer, 0) for kmer in kmers]
        mutated_kmer_counts = [mutated_kmers.get(kmer, 0) for kmer in kmers]

        plt.bar(kmers, normal_kmer_counts, width=0.4, label='Normal', align='center')
        plt.bar(kmers, mutated_kmer_counts, width=0.4, label='Mutated', align='edge')
        plt.title('K-mer Frequency Comparison')
        plt.ylabel('Frequency')
        plt.xticks(rotation=90)
        plt.legend()

        if save_path:
            plt.savefig(f"{save_path}_kmer_frequency.png")
        plt.show()
