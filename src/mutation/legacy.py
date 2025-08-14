# src/mutation/legacy.py
import warnings
warnings.warn(
    "mutation.legacy is deprecated and will be removed in the next release. "
    "Use mutation.engine.MutationEngine instead.",
    DeprecationWarning,
    stacklevel=2,
)
import random
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
from library import SequenceValidation, GeneUtils
from analysis_module import ContextModel

# DNA codon-to-amino-acid lookup table
# Used for synonymous/non-synonymous analysis
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}


class MarkovModelHandler:
    """Tiny helper that stores a 4x4 base transition matrix (A/C/G/T)."""
    BASE_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    def __init__(self):
        self.markov_matrix: Optional[np.ndarray] = None

    def set_markov_matrix(self, matrix: np.ndarray):
        if matrix.shape != (4, 4):
            raise ValueError("Markov matrix must be 4x4")
        # Normalize rows so they sum to 1.0 (like a good Markov citizen).
        self.markov_matrix = matrix / matrix.sum(axis=1, keepdims=True)

    def get_context_rate(self, prev_base: str, current_base: str) -> float:
        # No NumPy truthiness traps; check for None explicitly.
        if self.markov_matrix is None:
            return 1.0  # neutral multiplier if no model loaded
        p = self.BASE_INDEX.get(prev_base)
        c = self.BASE_INDEX.get(current_base)
        if p is None or c is None:
            return 1.0
        return float(self.markov_matrix[p, c])


class MutationEngine:
    """
    Runs different mutation types over a sequence slice.
    We log what we do because future debugging-you will thank us.
    """
    def __init__(self, sequence_structure, context_model: ContextModel = None):
        self.seqs = sequence_structure
        self.mutation_log: List[Dict] = []
        self.markov = MarkovModelHandler()
        self.context_model = context_model
        self._setup_defaults()

    def _setup_defaults(self):
        # Real genomes often have ~2:1 transitions:transversions
        self.set_ti_tv_ratio(2.1)
        # Simple hand-tuned context bumps (CpG hotspots, etc.)
        self.context_boost = {'CG': 10.0, 'TA': 0.5}
        # Motifs we’ll occasionally drop in for insertions
        self.motifs = ["AT", "CG", "TATA", "GCGC", "AATT"]

    def set_ti_tv_ratio(self, ratio: float):
        self.ti_prob = ratio / (ratio + 1.0)
        self.tv_prob = 1.0 / (ratio + 1.0)

    def _safe_sample(self, population_range, k):
        """Pick k unique positions safely without raising if not enough sites."""
        population = list(population_range)
        if k <= 0 or not population:
            return []
        k = min(k, len(population))
        return random.sample(population, k)

    def apply_mutations(self, accession_id: str, mutation_rate: float,
                        mutation_type: str, start: int, end: int) -> str:
        dispatch = {
            "point": self._point_mutations,
            "insertion": self._insertion_mutation,
            "deletion": self._deletion_mutation,
            "substitution": self._substitution_mutation,
            "inversion": self._inversion_mutation,
            "duplication": self._duplication_mutation,
            "translocation": self._translocation_mutation
        }
        sequence = self.seqs.get_sequence(accession_id)

        if mutation_type not in dispatch:
            raise ValueError(f"Unknown mutation type: {mutation_type}")

        # Be kind to callers: clamp and sanity-check the window.
        start = max(1, int(start))
        end = min(len(sequence), int(end))
        if end < start:
            raise ValueError(f"Invalid window: start={start}, end={end}")

        return dispatch[mutation_type](sequence, float(mutation_rate), start, end)

    # -------------------------
    # Point-like mutations
    # -------------------------
    def _point_mutations(self, sequence: str, rate: float,
                         start: int, end: int) -> str:
        """
        Apply per-base mutations with context boosts and optional Markov weighting.
        We also support a trinucleotide context model if provided.
        """
        seq = np.array(list(sequence), dtype='U1')
        positions = np.arange(start - 1, end)

        # Base mutation probability per site (as fraction)
        rates = np.full(len(positions), rate / 100.0, dtype=float)

        # Boost by simple dinucleotide context and optional Markov prev->curr transition prob
        for i, pos in enumerate(positions):
            if pos > 0:
                context2 = (seq[pos - 1] + seq[pos]).upper()
                rates[i] *= self.context_boost.get(context2, 1.0)
                rates[i] *= self.markov.get_context_rate(seq[pos - 1].upper(), seq[pos].upper())

        # Keep probabilities sane
        rates = np.clip(rates, 0.0, 1.0)

        # Sample which positions mutate
        mutate_mask = np.random.rand(len(positions)) < rates

        for idx in positions[mutate_mask]:
            old_base = str(seq[idx]).upper()

            # If we have a context model and a full 3-mer context, use it
            if self.context_model and 0 < idx < len(seq) - 1:
                tri_context = f"{seq[idx - 1]}{seq[idx]}{seq[idx + 1]}".upper()
                new_base = self.context_model.sample_mutation(tri_context, old_base)
            else:
                new_base = self._mutate_base(old_base)

            seq[idx] = new_base
            kind = self._classify_kind(old_base, str(new_base))

            self.mutation_log.append({
                "type": "point",
                "kind": kind,  # "transition" or "transversion" (or "none" in weird cases)
                "position": int(idx),
                "before": old_base,
                "after": str(new_base),
                **({"context": tri_context} if self.context_model and 0 < idx < len(seq) - 1 else {})
            })

        return ''.join(seq)

    def _mutate_base(self, base: str) -> str:
        """Flip a single base respecting the configured Ti/Tv ratio."""
        if random.random() < self.ti_prob:
            return self._transition(base)
        return self._transversion(base)

    def _transition(self, base: str) -> str:
        return {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C', 'U': 'C'}.get(base, base)

    def _transversion(self, base: str) -> str:
        options = {'A': ['C', 'T'], 'G': ['C', 'T'], 'C': ['A', 'G'], 'T': ['A', 'G'], 'U': ['A', 'G']}
        return random.choice(options.get(base, [base]))

    def _classify_kind(self, a: str, b: str) -> str:
        """Label a→b as transition/transversion for reporting."""
        a = a.upper(); b = b.upper()
        if a == b:
            return "none"
        ti_pairs = {frozenset(("A", "G")), frozenset(("C", "T"))}
        # Be nice to RNA sequences: treat U like T for classification
        if "U" in (a, b):
            a = 'T' if a == 'U' else a
            b = 'T' if b == 'U' else b
        return "transition" if frozenset((a, b)) in ti_pairs else "transversion"

    # -------------------------
    # Structural-ish mutations
    # -------------------------
    def _insertion_mutation(self, sequence: str, rate: float,
                            start: int, end: int) -> str:
        """Insert short motifs at selected positions."""
        mutated = list(sequence)
        window_len = max(0, end - (start - 1))
        num_insertions = int((end - start + 1) * (rate / 100.0))

        # Try to bias toward AT-rich windows; fall back to random sites
        rel_sites = self._find_at_rich_regions(sequence[start - 1:end])
        if not rel_sites:
            rel_sites = self._safe_sample(range(0, window_len), num_insertions)

        abs_sites = [start - 1 + s for s in rel_sites]
        for site in sorted(abs_sites, reverse=True):
            motif = random.choice(self.motifs)
            # Insert as characters, not a single string element
            mutated[site:site] = list(motif)
            self.mutation_log.append({
                "type": "insertion",
                "position": int(site),
                "value": motif
            })
        return ''.join(mutated)

    def _deletion_mutation(self, sequence: str, rate: float,
                           start: int, end: int) -> str:
        """Delete bases, favoring repetitive regions (slippage model-ish)."""
        mutated = list(sequence)
        window_len = max(0, end - (start - 1))
        num_deletions = int((end - start + 1) * (rate / 100.0))

        rel_sites = self._find_repetitive_regions(sequence[start - 1:end])
        if not rel_sites:
            rel_sites = self._safe_sample(range(0, window_len), num_deletions)

        abs_sites = [start - 1 + s for s in rel_sites]
        for site in sorted(set(abs_sites), reverse=True):
            if 0 <= site < len(mutated):
                removed = mutated[site]
                del mutated[site]
                self.mutation_log.append({
                    "type": "deletion",
                    "position": int(site),
                    "value": removed
                })
        return ''.join(mutated)

    def _substitution_mutation(self, sequence: str, rate: float,
                               start: int, end: int) -> str:
        """Alias to point mutations for caller convenience."""
        return self._point_mutations(sequence, rate, start, end)

    def _inversion_mutation(self, sequence: str, rate: float,
                             start: int, end: int) -> str:
        seg = sequence[start - 1:end]
        self.mutation_log.append({
            "type": "inversion",
            "start": int(start),
            "end": int(end),
            "length": len(seg)
        })
        return sequence[:start - 1] + seg[::-1] + sequence[end:]

    def _duplication_mutation(self, sequence: str, rate: float,
                              start: int, end: int) -> str:
        segment = sequence[start - 1:end]
        self.mutation_log.append({
            "type": "duplication",
            "start": int(start),
            "end": int(end),
            "length": len(segment)
        })
        return sequence[:end] + segment + sequence[end:]

    def _translocation_mutation(self, sequence: str, rate: float,
                                start: int, end: int) -> str:
        segment = sequence[start - 1:end]
        remaining = sequence[:start - 1] + sequence[end:]
        insert_pos = random.randint(0, len(remaining))
        self.mutation_log.append({
            "type": "translocation",
            "start": int(start),
            "end": int(end),
            "insert_pos": int(insert_pos),
            "length": len(segment)
        })
        return remaining[:insert_pos] + segment + remaining[insert_pos:]

    # -------------------------
    # Small analyses / helpers
    # -------------------------
    def analyze_synonymous(self, original: str, mutated: str) -> Dict:
        """Count how many codon changes preserve amino acid vs change it."""
        syn_count = 0
        non_syn_count = 0
        for i in range(0, len(original), 3):
            orig_codon = original[i:i + 3].upper()
            mut_codon = mutated[i:i + 3].upper()
            if len(orig_codon) == 3 and len(mut_codon) == 3:
                orig_aa = CODON_TABLE.get(orig_codon)
                mut_aa = CODON_TABLE.get(mut_codon)
                if orig_aa and mut_aa:
                    if orig_aa == mut_aa:
                        syn_count += 1
                    else:
                        non_syn_count += 1
        return {"synonymous": syn_count, "non_synonymous": non_syn_count}

    def _find_at_rich_regions(self, segment: str, threshold=0.6) -> List[int]:
        """Return indices (relative to the slice) that live inside AT-heavy windows."""
        sites: List[int] = []
        for i in range(len(segment) - 3):
            window = segment[i:i + 4]
            at_content = sum(1 for b in window if b in "AT") / 4.0
            if at_content >= threshold:
                sites.extend(range(i, i + 4))
        return list(set(sites))

    def _find_repetitive_regions(self, segment: str, min_repeats=3) -> List[int]:
        """Heuristic: collect indices participating in short repeats (2..5-mer × repeats)."""
        sites = set()
        for i in range(len(segment) - 1):
            for rlen in range(2, 6):
                motif = segment[i:i + rlen]
                if motif and motif * min_repeats in segment:
                    sites.update(range(i, i + rlen))
        return list(sites)

    def get_mutation_log(self) -> List[Dict]:
        return self.mutation_log

    def generate_report(self) -> Dict:
        """Tiny summary for dashboards."""
        ti = sum(1 for m in self.mutation_log if m.get('kind') == 'transition')
        tv = sum(1 for m in self.mutation_log if m.get('kind') == 'transversion')
        return {
            "total_mutations": len(self.mutation_log),
            "ti_tv_ratio": (ti / tv) if tv else float('inf') if ti > 0 else 0.0,
            "hotspots": self._find_hotspots(),
            "sequence_length": self.seqs.get_full_sequence_length()
        }

    def _find_hotspots(self, window=10) -> List[int]:
        """
        Look for clusters of point mutations within +/- window.
        It’s a quick heuristic; feel free to replace with a nicer kernel later.
        """
        positions = [m.get("position") for m in self.mutation_log if "position" in m]
        positions = [p for p in positions if p is not None]
        hotspots = []
        for i, p in enumerate(positions):
            neighborhood = positions[max(0, i - 2): i + 3]
            if any(abs(p - o) <= window for o in neighborhood if o is not None):
                hotspots.append(p)
        return hotspots
