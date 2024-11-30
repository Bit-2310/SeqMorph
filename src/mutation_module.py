import random
import numpy as np
from gene_library import SequenceValidation  # Import the specific function used
import sequence_structure as SeqStruct  # Rename to something more concise if widely used
#------------------------------------------------------------
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
#------------------------------------------------------------
class MarkovModelHandler:
    def __init__(self):
        # Initialize the Markov matrix to None by default
        self.markov_matrix = None

    def set_markov_matrix(self, markov_matrix):
        """
        Set the Markov matrix for calculating mutation rates dynamically.
        """
        if not isinstance(markov_matrix, np.ndarray):
            raise TypeError("Markov matrix must be a NumPy ndarray.")
        if markov_matrix.shape != (4, 4):
            raise ValueError("Markov matrix must be a 4x4 matrix.")
        
        # Additional check to ensure that all values are in a valid probability range [0, 1]
        if not ((0 <= markov_matrix).all() and (markov_matrix <= 1).all()):
            raise ValueError("All elements of the Markov matrix must be between 0 and 1.")
        
        self.markov_matrix = markov_matrix

    def calculate_markov_mutation_rate(self, base, previous_base):
        """
        Calculate the mutation rate for the given base using the Markov model.
        """
        base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        # Check if the matrix and valid bases are provided
        if self.markov_matrix is None:
            return 0.0  # Default mutation rate if no Markov model is provided
        if previous_base not in base_to_index or base not in base_to_index:
            return 0.0  # Return default if bases are not valid
        
        return self.markov_matrix[base_to_index[previous_base], base_to_index[base]]
#------------------------------------------------------------
class SequenceAnalyzer:
    def __init__(self, codon_table):
        self.codon_table = codon_table

    def analyze_synonymous(self, original_sequence, mutated_sequence):
        """
        Analyze mutations for synonymous vs. non-synonymous changes.
        """
        synonymous_count = 0
        non_synonymous_count = 0
        # Traverse the sequences by codon (3 bases at a time)
        for i in range(0, len(original_sequence) - 2, 3):
            original_codon = original_sequence[i:i+3]
            mutated_codon = mutated_sequence[i:i+3]
            if len(original_codon) == 3 and len(mutated_codon) == 3:
                original_aa = self.codon_table.get(original_codon, None)
                mutated_aa = self.codon_table.get(mutated_codon, None)
                if original_aa and mutated_aa:
                    if original_aa == mutated_aa:
                        synonymous_count += 1
                    else:
                        non_synonymous_count += 1
        return synonymous_count, non_synonymous_count

    def generate_mutation_report(self, synonymous_count, non_synonymous_count):
        """
        Generate a detailed mutation report.
        """
        return {
            "synonymous_mutations": synonymous_count,
            "non_synonymous_mutations": non_synonymous_count
        }
#------------------------------------------------------------
class MutationUtils:
    def __init__(self, sequence_structure):
        """
        Initializes the MutationUtils with a Sequence_Structure object.
        :param sequence_structure: The Sequence_Structure object that contains the sequence data.
        """
        self.seqs = sequence_structure
        self.mutation_log = []  # Log of all mutations applied to this sequence
        self.synonymous_count = 0
        self.non_synonymous_count = 0
        self.markov_handler = MarkovModelHandler()  # For managing Markov Matrix related operations
        self.transition_prob = 0.5  # Default probability for transitions
        self.transversion_prob = 0.5  # Default probability for transversions

    def set_markov_matrix(self, markov_matrix):
        """
        Set the Markov matrix to be used in mutations.
        """
        self.markov_handler.set_markov_matrix(markov_matrix)
    
    def validate_sequence(self, accession_id):
        """
        Validate the sequence using utility libraries.
        :param accession_id: The ID of the sequence to validate.
        :return: True if valid, else raises ValueError.
        """
        sequence = self.seqs.get_sequence(accession_id)
        if not sequence:
            raise ValueError(f"Sequence {accession_id} not found.")
        
        seq_type = self.seqs.MetaData.get(accession_id, {}).get('Type')
        if not seq_type:
            raise ValueError(f"Invalid sequence type for {accession_id}.")
        
        # Validate by sequence type (DNA, RNA, Protein)
        return GeneUtils.SequenceValidation.validate(sequence, seq_type)

    def mutation_inputs(self, seqs, accession_id):
        """
        Gather mutation inputs (mutation rate, sequence range, and mutation type) from the user.
        :param seqs: The Sequence_Structure object containing sequence data.
        :return: mutation_rate, range (start_pos, end_pos), mutation_type
        """
        # Step 1: Rate of Mutation (default to 5% if no input is provided)
        mutation_rate = input("Enter the mutation rate(%)[default is 5%]: ")
        mutation_rate = float(mutation_rate) if mutation_rate else 5.0  # Default is 5%

        # Step 2: Range of Sequence (start and end positions)
        seq_length = seqs.get_full_sequence_length(accession_id)
        print(f"Full sequence length: {seq_length} bp")
        start_pos = int(input(f"Enter the start position (1 to {seq_length}): "))
        end_pos = int(input(f"Enter the end position (start to {seq_length}): "))
        
        if start_pos < 1 or end_pos > seq_length or start_pos > end_pos:
            print("Invalid positions entered. Exiting.")
            return None, None, None
        
        # Step 3: Type of Mutation
        print("Select mutation type:\n1. Point Mutation\n2. Insertion\n3. Deletion\n4. Substitution\n5. Inversion\n6. Duplication\n7. Translocation")
        mutation_type = input("Enter mutation type (1, 2, 3, 4, 5, 6, or 7): ")
        
        if mutation_type not in ["1", "2", "3", "4", "5", "6", "7"]:
            print("Invalid mutation type selected. Exiting.")
            return None, None, None
        
        mutation_type_map = {
            "1": "point",
            "2": "insertion",
            "3": "deletion",
            "4": "substitution",
            "5": "inversion",
            "6": "duplication",
            "7": "translocation"
        }
        mutation_type = mutation_type_map[mutation_type]
        
        return mutation_rate, (start_pos, end_pos), mutation_type

    def apply_mutation(self, accession_id, mutation_rate, start_pos, end_pos, mutation_type):
        """
        Apply mutations to a specified range of a sequence.
        """
        sequence = self.seqs.get_sequence(accession_id)

        if sequence is None:
            raise ValueError(f"Sequence {accession_id} not found.")
        
        if start_pos < 1 or end_pos > self.seqs.get_full_sequence_length(accession_id):
            raise ValueError(f"Invalid range: {start_pos} - {end_pos} for sequence {accession_id}.")

        mutation_methods = {
            "point": self.apply_point_mutation,
            "insertion": self.apply_insertion_mutation,
            "deletion": self.apply_deletion_mutation,
            "substitution": self.apply_substitution_mutation,
            "inversion": self.apply_inversion_mutation,
            "duplication": self.apply_duplication_mutation,
            "translocation": self.apply_translocation_mutation
        }

        if mutation_type in mutation_methods:
            mutated_sequence = mutation_methods[mutation_type](sequence, mutation_rate, start_pos, end_pos)
        else:
            raise ValueError("Invalid mutation type.")

        # Log mutation details
        self._log_mutation(accession_id, mutated_sequence, mutation_rate, start_pos, end_pos, mutation_type)

        # Update metadata
        self.seqs.update_metadata(accession_id, mutated_sequence)

        return mutated_sequence

    def apply_point_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply a point mutation to the sequence within the specified range, considering transition and transversion probabilities.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The start position of mutation.
        :param end_pos: The end position of mutation.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        
        for i in range(start_pos - 1, end_pos):  # Adjust for 0-indexed position
            # Calculate dynamic mutation rate if Markov matrix is provided
            previous_base = sequence[i - 1] if i > 0 else None
            dynamic_mutation_rate = (
                self.markov_handler.calculate_markov_mutation_rate(mutated_sequence[i], previous_base) 
                if self.markov_handler.markov_matrix is not None 
                else mutation_rate / 100.0
            )
            
            # Apply mutation based on the calculated rate
            if random.random() < dynamic_mutation_rate:
                if random.random() < self.transition_prob:
                    # Apply transition mutation
                    mutated_sequence[i] = self.get_transition_base(mutated_sequence[i])
                else:
                    # Apply transversion mutation
                    mutated_sequence[i] = self.get_transversion_base(mutated_sequence[i])
                    
        return "".join(mutated_sequence)

    def get_transition_base(self, base):
        """
        Get a valid transition base for the given nucleotide.
        :param base: The original nucleotide base (A, T, G, C).
        :return: The transitioned base.
        """
        transition_map = {
            'A': 'G',  # A transitions to G
            'G': 'A',  # G transitions to A
            'C': 'T',  # C transitions to T
            'T': 'C'   # T transitions to C
        }
        return transition_map.get(base, base)

    def get_transversion_base(self, base):
        """
        Get a valid transversion base for the given nucleotide.
        :param base: The original nucleotide base (A, T, G, C).
        :return: The transversioned base.
        """
        transversion_map = {
            'A': ['C', 'T'],  # A transverts to C or T
            'G': ['C', 'T'],  # G transverts to C or T
            'C': ['A', 'G'],  # C transverts to A or G
            'T': ['A', 'G']   # T transverts to A or G
        }
        if base in transversion_map:
            return random.choice(transversion_map[base])
        return base

    def apply_inversion_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply inversion mutation to a specified range of the sequence.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage). (This parameter is not used for inversion)
        :param start_pos: The starting position in the sequence to begin inversion.
        :param end_pos: The ending position of the sequence for inversion.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        mutated_sequence[start_pos - 1:end_pos] = reversed(mutated_sequence[start_pos - 1:end_pos])
        return "".join(mutated_sequence)

    def apply_duplication_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply duplication mutation by duplicating a segment of the sequence and appending it to the end.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The starting position in the sequence to begin duplication.
        :param end_pos: The ending position of the sequence for duplication.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        duplicated_segment = mutated_sequence[start_pos - 1:end_pos]
        mutated_sequence = mutated_sequence[:end_pos] + duplicated_segment + mutated_sequence[end_pos:]
        return "".join(mutated_sequence)

    def apply_translocation_mutation(self, sequence, mutation_rate ,start_pos, end_pos):
        """
        Apply translocation mutation by moving a segment from one part of the sequence to another.
        :param sequence: The sequence to mutate.
        :param start_pos: The starting position in the sequence to begin translocation.
        :param end_pos: The ending position of the sequence for translocation.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        segment = mutated_sequence[start_pos - 1:end_pos]
        del mutated_sequence[start_pos - 1:end_pos]
        insertion_pos = random.randint(0, len(mutated_sequence))
        mutated_sequence = mutated_sequence[:insertion_pos] + segment + mutated_sequence[insertion_pos:]
        return "".join(mutated_sequence)

    def apply_point_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply a point mutation to the sequence within the specified range, considering transition and transversion probabilities.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The start position of mutation.
        :param end_pos: The end position of mutation.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        for i in range(start_pos - 1, end_pos):  # Adjust for 0-indexed position
            previous_base = sequence[i - 1] if i > 0 else None
            dynamic_mutation_rate = (
                self.markov_handler.calculate_markov_mutation_rate(mutated_sequence[i], previous_base)
                if self.markov_handler.markov_matrix is not None
                else mutation_rate / 100.0
            )
            if random.random() < dynamic_mutation_rate:
                if random.random() < self.transition_prob:
                    mutated_sequence[i] = self.get_transition_base(mutated_sequence[i])
                else:
                    mutated_sequence[i] = self.get_transversion_base(mutated_sequence[i])
        return "".join(mutated_sequence)

    def apply_insertion_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply insertion mutations to the sequence within the specified range, focusing on AT-rich regions or repetitive motifs.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The start position of mutation.
        :param end_pos: The end position of mutation.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        possible_range_size = end_pos - start_pos + 1
        num_insertions = int(possible_range_size * (mutation_rate / 100.0))

        # Ensure we don't try to insert more than possible
        num_insertions = min(num_insertions, possible_range_size)

        if num_insertions > 0:
            # Step 1: Identify potential AT-rich regions or repetitive motifs
            insertion_sites = self.find_at_rich_regions(sequence[start_pos-1:end_pos])
            
            # Step 2: Select insertion sites, with a higher probability for AT-rich regions
            if insertion_sites:
                selected_sites = random.sample(insertion_sites, min(num_insertions, len(insertion_sites)))
            else:
                # If no AT-rich regions are found, insert randomly within the range
                selected_sites = random.sample(range(start_pos - 1, end_pos), num_insertions)

            # Step 3: Insert randomly chosen nucleotides or motifs at the selected sites
            for site in selected_sites:
                insert_base = self.get_biologically_relevant_insert()
                mutated_sequence.insert(site, insert_base)

        return "".join(mutated_sequence)

    def find_at_rich_regions(self, sequence_segment, threshold=0.6):
        """
        Identify regions with high AT content in the given sequence segment.
        :param sequence_segment: The segment of the sequence to analyze.
        :param threshold: Minimum AT content required to consider the region AT-rich.
        :return: List of indices where AT-rich regions are found.
        """
        at_rich_sites = []
        sequence_length = len(sequence_segment)

        # Scan sequence in windows of four bases
        for i in range(sequence_length - 3):
            window = sequence_segment[i:i + 4]
            at_content = sum(1 for base in window if base in "AT") / 4.0
            if at_content >= threshold:
                at_rich_sites.extend(range(i, i + 4))

        # Remove duplicates and return unique sites
        return list(set(at_rich_sites))
    
    def get_biologically_relevant_insert(self):
        """
        Get a biologically relevant motif or nucleotide for insertion.
        :return: A random nucleotide or a short motif for insertion.
        """
        motifs = ["AT", "CG", "TATA", "GCGC", "AATT"]  # Common motifs for insertion
        if random.random() < 0.5:
            # Insert a single nucleotide
            return random.choice("ATGC")
        else:
            # Insert a motif
            return random.choice(motifs)
 
    def apply_deletion_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply deletion mutations to the sequence within the specified range, focusing on repetitive regions.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The start position of mutation.
        :param end_pos: The end position of mutation.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        possible_range_size = end_pos - start_pos + 1
        num_deletions = int(possible_range_size * (mutation_rate / 100.0))

        # Ensure we don't try to delete more than possible
        num_deletions = min(num_deletions, possible_range_size)

        if num_deletions > 0:
            # Step 1: Identify potential repetitive regions to focus deletions on
            repetitive_sites = self.find_repetitive_regions(sequence[start_pos-1:end_pos])
            
            # Step 2: Select deletion sites with a higher probability in repetitive regions
            if repetitive_sites:
                deletion_sites = random.sample(repetitive_sites, min(num_deletions, len(repetitive_sites)))
            else:
                # If no repetitive regions found, delete randomly within the range
                deletion_sites = random.sample(range(start_pos - 1, end_pos), num_deletions)
            
            # Step 3: Delete the selected sites
            for site in sorted(deletion_sites, reverse=True):  # Deleting from the end to avoid index issues
                del mutated_sequence[site]

        return "".join(mutated_sequence)

    def find_repetitive_regions(self, sequence_segment, min_repeats=3):
        """
        Identify regions with repetitive motifs in the given sequence segment.
        :param sequence_segment: The segment of the sequence to analyze.
        :param min_repeats: Minimum number of repeat units to be considered repetitive.
        :return: List of indices where repetitive regions are found.
        """
        repetitive_sites = set()
        sequence_length = len(sequence_segment)

        # Scan for repeated patterns of length 2-5 bases
        for i in range(sequence_length - 1):
            for repeat_length in range(2, 6):
                if i + repeat_length <= sequence_length:
                    repeat_unit = sequence_segment[i:i + repeat_length]
                    # Check if repeat_unit appears consecutively min_repeats times
                    if sequence_segment.count(repeat_unit * min_repeats) > 0:
                        repetitive_sites.update(range(i, i + repeat_length))

        return list(repetitive_sites)

    def apply_substitution_mutation(self, sequence, mutation_rate, start_pos, end_pos):
        """
        Apply substitution mutations to the sequence within the specified range.
        :param sequence: The sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The start position of mutation.
        :param end_pos: The end position of mutation.
        :return: The mutated sequence.
        """
        mutated_sequence = list(sequence)
        for i in range(start_pos - 1, end_pos):  # Adjust for 0-indexed position
            if random.random() < mutation_rate / 100.0:
                mutated_sequence[i] = random.choice([b for b in "ATGC" if b != mutated_sequence[i]])
        return "".join(mutated_sequence)

    def mutate_base(self, base):
        """
        Mutates a single base (example for DNA: A -> T, G -> C, etc.).
        :param base: The base to mutate.
        :return: The mutated base.
        """
        mutation_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'
        }
        return mutation_map.get(base, base)

    def apply_codon_level_mutation(self, accession_id, mutation_rate, start_pos, end_pos):
        """
        Apply codon-level point mutations to the sequence within the specified range.
        :param accession_id: The ID of the sequence to mutate.
        :param mutation_rate: The rate of mutation (percentage).
        :param start_pos: The start position of mutation.
        :param end_pos: The end position of mutation.
        :return: The mutated sequence and a count of synonymous vs non-synonymous mutations.
        """
        sequence = self.seqs.get_sequence(accession_id)
        mutated_sequence = list(sequence)
        num_codons = (end_pos - start_pos + 1) // 3

        for codon_index in range(num_codons):
            codon_start = start_pos - 1 + codon_index * 3
            codon_end = codon_start + 3
            original_codon = sequence[codon_start:codon_end]

            # Apply mutation based on mutation rate
            if random.random() < (mutation_rate / 100.0):
                mutated_codon = self.mutate_codon(original_codon)
                mutated_sequence[codon_start:codon_end] = mutated_codon

                # Determine if the mutation is synonymous or non-synonymous
                original_aa = CODON_TABLE.get(original_codon, None)
                mutated_aa = CODON_TABLE.get(mutated_codon, None)

                if original_aa and mutated_aa:
                    if original_aa == mutated_aa:
                        self.synonymous_count += 1
                    else:
                        self.non_synonymous_count += 1

        return "".join(mutated_sequence)
    
    def set_transition_transversion_probs(self, transition_prob, transversion_prob):
        """
        Set the probabilities for transition and transversion mutations.
        :param transition_prob: Probability for transition mutations (0 to 1).
        :param transversion_prob: Probability for transversion mutations (0 to 1).
        """
        epsilon = 1e-6
        if abs(transition_prob + transversion_prob - 1.0) > epsilon:
            raise ValueError("Transition and transversion probabilities must sum to 1.")
        self.transition_prob = transition_prob
        self.transversion_prob = transversion_prob
    
    def mutate_codon(self, codon):
        """
        Mutate a given codon by randomly changing one of its bases.
        :param codon: The original codon (3 bases).
        :return: The mutated codon.
        """
        codon_list = list(codon)
        base_index = random.randint(0, 2)  # Choose one of the three bases to mutate
        original_base = codon_list[base_index]

        # Determine whether it's a transition or transversion
        if random.random() < self.transition_prob:
            # Transition mutation
            codon_list[base_index] = self.get_transition_base(original_base)
        else:
            # Transversion mutation
            codon_list[base_index] = self.get_transversion_base(original_base)

        return "".join(codon_list)

    def _log_mutation(self, accession_id, mutated_sequence, mutation_rate, start_pos, end_pos, mutation_type):
        """
        Log mutation details.
        :param accession_id: The ID of the sequence.
        :param mutated_sequence: The sequence after mutation.
        :param mutation_rate: The rate at which mutation was applied.
        :param start_pos: Starting position of mutation.
        :param end_pos: Ending position of mutation.
        :param mutation_type: Type of mutation applied.
        """
        mutation_entry = {
            "accession_id": accession_id,
            "mutation_rate": mutation_rate,
            "mutation_type": mutation_type,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "mutated_sequence": mutated_sequence,
        }
        self.mutation_log.append(mutation_entry)

    def get_mutation_log(self):
        """
        Retrieve mutation logs.
        """
        return self.mutation_log