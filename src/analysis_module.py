import os
from library import SequenceAnalysis as SeqA
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
try:
    from Bio import Align as _BioAlign
except Exception:
    _BioAlign = None


# =======================
# Context-based mutation model
# =======================
class ContextModel:
    """
    Tiny helper to load a trinucleotide context → mutation probability table.
    Think COSMIC-style signatures but simpler. If you don’t give us a table,
    we politely do nothing (fallbacks still work).
    """
    def __init__(self, table_path: Optional[str] = None, boost_cpg: bool = False):
        self.context_probs: Dict[tuple, float] = {}
        self.boost_cpg = boost_cpg
        if table_path:
            self.load_table(table_path)

    def load_table(self, path: str) -> None:
        """
        CSV format:
        context,from,to,probability
        e.g. ACA,A,C,0.0025
        """
        import csv
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                context = row["context"].upper()
                base_from = row["from"].upper()
                base_to = row["to"].upper()
                prob = float(row["probability"])
                self.context_probs[(context, base_from, base_to)] = prob

    def get_probs(self, context: str, base_from: str) -> Dict[str, float]:
        """
        Return a {to_base: prob} dict for a context and from-base.
        If you turned on CpG boosting, we nudge those a bit.
        """
        probs: Dict[str, float] = {}
        for base_to in ["A", "C", "G", "T"]:
            if base_to == base_from:
                continue
            p = self.context_probs.get((context, base_from, base_to), 0.0)
            # CpG hotspot = current base is C and next base is G (prev + curr + next)
            if self.boost_cpg and base_from == "C" and len(context) >= 3 and context[1:3] == "CG":
                p *= 1.5
            probs[base_to] = p
        return probs

    def sample_mutation(self, context: str, base_from: str) -> str:
        """
        Draw a base according to the table. If we don’t have data, pick a random valid base.
        """
        import random
        probs = self.get_probs(context, base_from)
        total = sum(probs.values())
        if total <= 0:
            return random.choice([b for b in "ACGT" if b != base_from])
        bases, weights = zip(*[(b, p / total) for b, p in probs.items()])
        return random.choices(bases, weights=weights, k=1)[0]


class SequenceAnalysisReport:
    """
    Compare a “normal” vs “mutated” sequence and produce a structured report.
    Supports DNA / RNA / Protein. For proteins, we skip ORFs (because, well, they’re for nucleotides).
    """
    def __init__(
        self,
        normal_sequence: str,
        mutated_sequence: str,
        sequence_type: str = "DNA",
        enable_alignment: bool = True
    ):
        if not normal_sequence or not mutated_sequence:
            raise ValueError("Sequences cannot be empty.")
        if sequence_type not in ["DNA", "RNA", "Protein"]:
            raise ValueError("Invalid sequence type. Must be DNA, RNA, or Protein.")

        self.normal_sequence = normal_sequence
        self.mutated_sequence = mutated_sequence
        self.sequence_type = sequence_type
        # Only enable alignment if caller asked for it AND Biopython Align is available
        self.enable_alignment = bool(enable_alignment and _BioAlign is not None)
        self._cached_results: Dict[str, Dict] = {}  # because recomputing is boring

    # -------- Alignment (desktop/server only) --------
    def align_sequences(self) -> dict:
        """
        Global alignment using Biopython’s PairwiseAligner.
        If not available or disabled, return a small note dict.
        """
        if not self.enable_alignment:
            reason = "disabled by flag" if _BioAlign is not None else "biopython not installed"
            return {"note": f"alignment skipped ({reason})"}

        aligner = _BioAlign.PairwiseAligner()
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

    # -------- Main comparison --------
    def compare_sequences(self) -> dict:
        """
        Crunch the key metrics and return a dict the UI can love.
        """
        if "report" in self._cached_results:
            return self._cached_results["report"]

        # Use correct ORF logic for DNA vs RNA; skip for proteins
        if self.sequence_type == "Protein":
            normal_orfs: List[str] = []
            mutated_orfs: List[str] = []
        else:
            normal_orfs = SeqA.orf_detection(self.normal_sequence, seq_type=self.sequence_type)
            mutated_orfs = SeqA.orf_detection(self.mutated_sequence, seq_type=self.sequence_type)

        # k-mers: if sequences are too short, just hand back {}
        normal_kmers = SeqA.kmer_count(self.normal_sequence) if len(self.normal_sequence) >= 3 else {}
        mutated_kmers = SeqA.kmer_count(self.mutated_sequence) if len(self.mutated_sequence) >= 3 else {}

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
                "Normal": normal_kmers,
                "Mutated": mutated_kmers
            },
            # Alignment can be expensive and not always available; guard it.
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

        # Cache before we add the “stats” block so callers can reuse the heavy bits
        self._cached_results["report"] = report

        # Keep the stats honest and interpretable.
        # - GC Content: show simple delta (no bogus t-test with n=1)
        # - Chi-square for nt/k-mers only if counts are sane
        report["Statistical Tests"] = {
            "GC Content": {
                "delta_pct": report["GC Content"]["Mutated"] - report["GC Content"]["Normal"]
            },
            "Nucleotide Frequency": self._compare_nucleotide_freq(report),
            "K-mer Frequency": self._compare_kmer_freq(report),
        }
        return report

    # -------- Stats helpers --------
    def _compare_nucleotide_freq(self, report: Dict) -> Dict:
        """
        Chi-square on raw nucleotide counts. Drop categories that are zero in BOTH rows
        (e.g., 'U' for DNA) to avoid zero-expected cells. If it still fails, return a note.
        """
        normal = report["Nucleotide Frequency"]["Normal"]
        mutated = report["Nucleotide Frequency"]["Mutated"]

        labels_all = ["A", "T", "G", "C", "U"]
        # keep any base that appears in at least one row
        labels = [b for b in labels_all if (normal.get(b, 0) + mutated.get(b, 0)) > 0]

        # need at least 2 categories with nonzero totals to run chi-square meaningfully
        if len(labels) < 2:
            return {"note": "not enough categories for chi-square", "labels": labels}

        normal_freq = [normal.get(b, 0) for b in labels]
        mutated_freq = [mutated.get(b, 0) for b in labels]

        # sanity: both rows must have some counts
        if sum(normal_freq) == 0 or sum(mutated_freq) == 0:
            return {"note": "insufficient counts for chi-square", "labels": labels}

        try:
            chi2, p, _, _ = chi2_contingency([normal_freq, mutated_freq])
            return {"chi2_stat": float(chi2), "p_value": float(p), "labels": labels}
        except Exception as e:
            return {"note": f"chi-square failed: {type(e).__name__}", "labels": labels}


    def _compare_kmer_freq(self, report: Dict) -> Dict:
        """
        Chi-square over the union of observed k-mers. If there are too many categories,
        this gets unwieldy; UIs should show “top-N” anyway.
        """
        normal_kmers: Dict[str, int] = report["K-mer Frequency"]["Normal"] or {}
        mutated_kmers: Dict[str, int] = report["K-mer Frequency"]["Mutated"] or {}
        if not normal_kmers and not mutated_kmers:
            return {"note": "no k-mers (sequence too short)"}

        kmers = list(set(normal_kmers).union(mutated_kmers))
        normal_counts = [int(normal_kmers.get(k, 0)) for k in kmers]
        mutated_counts = [int(mutated_kmers.get(k, 0)) for k in kmers]

        if sum(normal_counts) == 0 or sum(mutated_counts) == 0:
            return {"note": "insufficient k-mer counts for chi-square"}

        chi2, p, _, _ = chi2_contingency([normal_counts, mutated_counts])
        return {"chi2_stat": float(chi2), "p_value": float(p), "num_kmers": len(kmers)}

    # -------- Plotting (desktop/server only) --------
    def plot_comparison(self, report: Optional[Dict] = None, save_path: Optional[str] = None):
        """
        Create simple plots. If save_path is given, we save files and close the figures
        (no blocking popups); otherwise we call plt.show().
        """
        if report is None:
            report = self.compare_sequences()

        if save_path:
            out_dir = os.path.dirname(save_path) or "."
            os.makedirs(out_dir, exist_ok=True)

        self._plot_gc_content(report, save_path)
        self._plot_nucleotide_frequency(report, save_path)
        self._plot_kmer_frequency(report, save_path)

    def _plot_gc_content(self, report: Dict, save_path: Optional[str] = None):
        labels = ['Normal', 'Mutated']
        values = [report["GC Content"]["Normal"], report["GC Content"]["Mutated"]]

        plt.figure()
        plt.bar(labels, values)
        plt.title('GC Content Comparison')
        plt.ylabel('% GC Content')

        if save_path:
            plt.savefig(f"{save_path}_gc_content.png", bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def _plot_nucleotide_frequency(self, report: Dict, save_path: Optional[str] = None):
        normal_freq = report["Nucleotide Frequency"]["Normal"]
        mutated_freq = report["Nucleotide Frequency"]["Mutated"]
        labels = ['A', 'T', 'G', 'C', 'U']  # We include U for RNA; it will be ~0 for DNA.

        normal_values = [normal_freq.get(base, 0) for base in labels]
        mutated_values = [mutated_freq.get(base, 0) for base in labels]

        x = range(len(labels))
        width = 0.35
        fig, ax = plt.subplots()
        ax.bar(x, normal_values, width, label='Normal')
        ax.bar([p + width for p in x], mutated_values, width, label='Mutated')
        ax.set_ylabel('Count')
        ax.set_title('Nucleotide Frequency Comparison')
        ax.set_xticks([p + width / 2 for p in x])
        ax.set_xticklabels(labels)
        ax.legend()

        if save_path:
            fig.savefig(f"{save_path}_nucleotide_frequency.png", bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()

    def _plot_kmer_frequency(self, report: Dict, save_path: Optional[str] = None):
        normal_kmers = report["K-mer Frequency"]["Normal"] or {}
        mutated_kmers = report["K-mer Frequency"]["Mutated"] or {}

        kmers = list(set(normal_kmers).union(mutated_kmers))
        # If there are too many k-mers this gets ugly; callers should slice top-N before plotting.
        normal_kmer_counts = [normal_kmers.get(kmer, 0) for kmer in kmers]
        mutated_kmer_counts = [mutated_kmers.get(kmer, 0) for kmer in kmers]

        plt.figure(figsize=(max(6, len(kmers) * 0.2), 4))
        plt.bar(kmers, normal_kmer_counts, width=0.4, label='Normal', align='center')
        plt.bar(kmers, mutated_kmer_counts, width=0.4, label='Mutated', align='edge')
        plt.title('K-mer Frequency Comparison')
        plt.ylabel('Count')
        plt.xticks(rotation=90)
        plt.legend()

        if save_path:
            plt.savefig(f"{save_path}_kmer_frequency.png", bbox_inches="tight")
            plt.close()
        else:
            plt.show()
