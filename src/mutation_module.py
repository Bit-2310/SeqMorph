import random
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple
from gene_library import SequenceValidation, GeneUtils

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
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

class MarkovModelHandler:
    BASE_INDEX = {'A':0, 'C':1, 'G':2, 'T':3}
    
    def __init__(self):
        self.markov_matrix = None

    def set_markov_matrix(self, matrix: np.ndarray):
        if matrix.shape != (4, 4):
            raise ValueError("Markov matrix must be 4x4")
        self.markov_matrix = matrix / matrix.sum(axis=1, keepdims=True)

    def get_context_rate(self, prev_base: str, current_base: str) -> float:
        if not self.markov_matrix or prev_base not in self.BASE_INDEX:
            return 0.0
        return self.markov_matrix[self.BASE_INDEX[prev_base], self.BASE_INDEX[current_base]]

class MutationEngine:
    def __init__(self, sequence_structure):
        self.seqs = sequence_structure
        self.mutation_log = []
        self.markov = MarkovModelHandler()
        self._setup_defaults()

    def _setup_defaults(self):
        self.set_ti_tv_ratio(2.1)
        self.context_boost = {'CG': 10.0, 'TA': 0.5}
        self.motifs = ["AT", "CG", "TATA", "GCGC", "AATT"]

    def set_ti_tv_ratio(self, ratio: float):
        self.ti_prob = ratio / (ratio + 1)
        self.tv_prob = 1 / (ratio + 1)

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
        return dispatch[mutation_type](sequence, mutation_rate, start, end)

    def _point_mutations(self, sequence: str, rate: float,
                        start: int, end: int) -> str:
        seq = np.array(list(sequence), dtype='U1')
        positions = np.arange(start-1, end)
        
        # Context-aware mutation rates
        rates = np.full(len(positions), rate/100)
        for i, pos in enumerate(positions):
            if pos > 0:
                context = seq[pos-1] + seq[pos]
                rates[i] *= self.context_boost.get(context, 1.0)
            if self.markov.markov_matrix and pos > 0:
                rates[i] *= self.markov.get_context_rate(seq[pos-1], seq[pos])
        
        # Vectorized mutation application
        mutate_mask = np.random.rand(len(positions)) < rates
        for idx in positions[mutate_mask]:
            seq[idx] = self._mutate_base(seq[idx])
        return ''.join(seq)

    def _mutate_base(self, base: str) -> str:
        if random.random() < self.ti_prob:
            return self._transition(base)
        return self._transversion(base)

    def _transition(self, base: str) -> str:
        return {'A':'G', 'G':'A', 'C':'T', 'T':'C', 'U':'C'}.get(base, base)

    def _transversion(self, base: str) -> str:
        options = {'A':['C','T'], 'G':['C','T'], 
                 'C':['A','G'], 'T':['A','G'], 'U':['A','G']}
        return random.choice(options.get(base, [base]))

    # Original mutation types preserved with optimizations
    def _insertion_mutation(self, sequence: str, rate: float,
                           start: int, end: int) -> str:
        mutated = list(sequence)
        num_insertions = int((end - start + 1) * (rate / 100))
        sites = self._find_at_rich_regions(sequence[start-1:end]) or \
                random.sample(range(start-1, end), num_insertions)
        
        for site in sorted(sites, reverse=True):
            mutated.insert(site, random.choice(self.motifs))
        return ''.join(mutated)

    def _deletion_mutation(self, sequence: str, rate: float,
                          start: int, end: int) -> str:
        mutated = list(sequence)
        num_deletions = int((end - start + 1) * (rate / 100))
        sites = self._find_repetitive_regions(sequence[start-1:end]) or \
                random.sample(range(start-1, end), num_deletions)
        
        for site in sorted(sites, reverse=True):
            del mutated[site]
        return ''.join(mutated)

    def _substitution_mutation(self, sequence: str, rate: float,
                              start: int, end: int) -> str:
        return self._point_mutations(sequence, rate, start, end)

    def _inversion_mutation(self, sequence: str, rate: float,
                           start: int, end: int) -> str:
        return sequence[:start-1] + \
               sequence[start-1:end][::-1] + \
               sequence[end:]

    def _duplication_mutation(self, sequence: str, rate: float,
                             start: int, end: int) -> str:
        segment = sequence[start-1:end]
        return sequence[:end] + segment + sequence[end:]

    def _translocation_mutation(self, sequence: str, rate: float,
                               start: int, end: int) -> str:
        segment = sequence[start-1:end]
        remaining = sequence[:start-1] + sequence[end:]
        insert_pos = random.randint(0, len(remaining))
        return remaining[:insert_pos] + segment + remaining[insert_pos:]

    # Original analysis methods
    def analyze_synonymous(self, original: str, mutated: str) -> Dict:
        syn_count = 0
        non_syn_count = 0
        for i in range(0, len(original), 3):
            orig_codon = original[i:i+3]
            mut_codon = mutated[i:i+3]
            if len(orig_codon) == 3 and len(mut_codon) == 3:
                orig_aa = CODON_TABLE.get(orig_codon)
                mut_aa = CODON_TABLE.get(mut_codon)
                if orig_aa and mut_aa:
                    if orig_aa == mut_aa: syn_count +=1
                    else: non_syn_count +=1
        return {"synonymous": syn_count, "non_synonymous": non_syn_count}

    # Helper methods from original code
    def _find_at_rich_regions(self, segment: str, threshold=0.6) -> List[int]:
        sites = []
        for i in range(len(segment) - 3):
            window = segment[i:i+4]
            at_content = sum(1 for b in window if b in "AT") / 4
            if at_content >= threshold:
                sites.extend(range(i, i+4))
        return list(set(sites))

    def _find_repetitive_regions(self, segment: str, min_repeats=3) -> List[int]:
        sites = set()
        for i in range(len(segment) - 1):
            for rlen in range(2, 6):
                if segment.count(segment[i:i+rlen] * min_repeats) > 0:
                    sites.update(range(i, i+rlen))
        return list(sites)

    def get_mutation_log(self) -> List[Dict]:
        return self.mutation_log

    def generate_report(self) -> Dict:
        ti = sum(1 for m in self.mutation_log if m['type'] == 'transition')
        tv = len(self.mutation_log) - ti
        return {
            "total_mutations": len(self.mutation_log),
            "ti_tv_ratio": ti / tv if tv else 0,
            "hotspots": self._find_hotspots(),
            "sequence_length": self.seqs.get_full_sequence_length()
        }

    def _find_hotspots(self, window=10) -> List[int]:
        positions = [m["position"] for m in self.mutation_log]
        return [p for i, p in enumerate(positions)
                if any(abs(p - o) <= window 
                      for o in positions[max(0,i-2):i+3])]