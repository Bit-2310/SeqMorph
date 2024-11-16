import random
from protein_library import ProteinUtils
from gene_library import GeneUtils
from input_module import InputHandler
class MutationUtils:
    """
    A utility class for handling mutations on sequences (DNA, RNA, Protein).
    """
    def __init__(self, sequence):
        """
        Initialize the MutationUtils class with a sequence.
        Automatically determines the sequence type and initializes the mutation log.
        """
        self.sequence = sequence
        self.sequence_type = GeneUtils.detect_sequence_type(sequence)
        self.mutation_log = []  # Log of all mutations applied to this sequence

    def validate_sequence(self):
        """
        Validate the sequence based on its type (DNA, RNA, or Protein).
        :return: Boolean indicating whether the sequence is valid.
        """
        if self.sequence_type in ["DNA", "RNA"]:
            return GeneUtils.validate(self.sequence)
        elif self.sequence_type == "Protein":
            return ProteinUtils.validate(self.sequence)
        else:
            raise ValueError(f"Unknown sequence type: {self.sequence_type}")

    def count_ambiguities(self):
        """
        Count the number of ambiguous bases in the sequence.
        :return: Integer count of ambiguous bases.
        """
        if self.sequence_type in ["DNA", "RNA"]:
            return GeneUtils.count_ambiguities(self.sequence)
        else:
            raise ValueError("Ambiguities are not applicable for protein sequences.")

    def apply_random_mutation(self, mutation_rate):
        """
        Apply random mutations to the sequence and log them.
        :param mutation_rate: Percentage of sequence to mutate.
        :return: Mutated sequence.
        """
        if self.sequence_type in ["DNA", "RNA"]:
            mutated_sequence = self._random_mutation_gene(self.sequence, mutation_rate)
        elif self.sequence_type == "Protein":
            mutated_sequence = ProteinUtils.random_mutation(self.sequence, mutation_rate)
        else:
            raise ValueError(f"Unknown sequence type: {self.sequence_type}")
        
        # Log the mutation
        self._log_mutation(
            mutation_type="Random Mutation",
            subtype=None,
            original_sequence=self.sequence,
            mutated_sequence=mutated_sequence,
            mutation_rate=mutation_rate
        )
        
        # Update the sequence with the mutated version
        self.sequence = mutated_sequence
        return mutated_sequence

    def translate(self, handle_ambiguity=True):
        """
        Translate the sequence if it is DNA or RNA.
        :param handle_ambiguity: Whether to handle ambiguous bases during translation.
        :return: Translated protein sequence.
        """
        if self.sequence_type in ["DNA", "RNA"]:
            return GeneUtils.translate(self.sequence, handle_ambiguity)
        else:
            raise ValueError("Translation is only applicable to DNA or RNA sequences.")

    def _random_mutation_gene(self, sequence, mutation_rate):
        """
        Apply random mutations to a DNA or RNA sequence.
        Handles ambiguous bases during mutation.
        :param sequence: DNA or RNA sequence.
        :param mutation_rate: Percentage of sequence to mutate.
        :return: Mutated sequence.
        """
        sequence = list(sequence.upper())
        num_mutations = max(1, int(len(sequence) * mutation_rate / 100))

        for _ in range(num_mutations):
            position = random.randint(0, len(sequence) - 1)
            original_base = sequence[position]
            
            if original_base in GeneUtils.AMBIGUITY_CODES:
                possible_bases = GeneUtils.resolve_ambiguity(original_base)
            else:
                possible_bases = list(GeneUtils.TRANSITION_MAP.keys()) if random.random() < 0.5 else list(GeneUtils.TRANSVERSION_MAP.keys())
            
            possible_bases.remove(original_base)
            sequence[position] = random.choice(possible_bases)
        
        return ''.join(sequence)

    def _log_mutation(self, mutation_type, subtype, original_sequence, mutated_sequence, mutation_rate=None, position=None):
        """
        Log the details of a mutation.
        :param mutation_type: The high-level mutation type (e.g., 'Point Mutation').
        :param subtype: The specific mutation subtype (e.g., 'Missense').
        :param original_sequence: The sequence before mutation.
        :param mutated_sequence: The sequence after mutation.
        :param mutation_rate: Mutation rate (for random mutations).
        :param position: Position of the mutation (if applicable).
        """
        mutation_entry = {
            "mutation_type": mutation_type,
            "subtype": subtype,
            "original_sequence": original_sequence,
            "mutated_sequence": mutated_sequence,
            "mutation_rate": mutation_rate,
            "position": position
        }
        self.mutation_log.append(mutation_entry)

    def get_mutation_log(self):
        """
        Retrieve the mutation log.
        :return: List of mutation log entries.
        """
        return self.mutation_log
    
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-Mutation Engine-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MutationEngine:
    mutation_log = []  # Log of all mutations performed

    class GeneMutations:
        class PointMutation:
            @staticmethod
            def silent(codon, transition_weight=0.7, transversion_weight=0.3):
                amino_acid = GeneUtils.CODON_TABLE_DNA.get(codon)
                if not amino_acid:
                    raise ValueError(f"Invalid codon: {codon}")
                synonymous_codons = [
                    key for key, value in GeneUtils.CODON_TABLE_DNA.items()
                    if value == amino_acid and key != codon
                ]
                if not synonymous_codons:
                    return codon  # No synonymous mutation possible
                return MutationEngine.GeneMutations.PointMutation.weighted_mutation(
                    codon, synonymous_codons, transition_weight, transversion_weight
                )

            @staticmethod
            def missense(codon, transition_weight=0.7, transversion_weight=0.3):
                original_amino_acid = GeneUtils.CODON_TABLE_DNA.get(codon)
                if not original_amino_acid:
                    raise ValueError(f"Invalid codon: {codon}")
                nonsynonymous_codons = [
                    key for key, value in GeneUtils.CODON_TABLE_DNA.items()
                    if value != original_amino_acid and value != 'STOP'
                ]
                if not nonsynonymous_codons:
                    return codon  # No missense mutation possible
                return MutationEngine.GeneMutations.PointMutation.weighted_mutation(
                    codon, nonsynonymous_codons, transition_weight, transversion_weight
                )

            @staticmethod
            def nonsense(codon, transition_weight=0.7, transversion_weight=0.3):
                if codon in GeneUtils.STOP_CODONS_DNA:
                    return codon  # Codon is already a stop codon
                possible_stop_codons = list(GeneUtils.STOP_CODONS_DNA - {codon})
                if not possible_stop_codons:
                    return codon  # No nonsense mutation possible
                return MutationEngine.GeneMutations.PointMutation.weighted_mutation(
                    codon, possible_stop_codons, transition_weight, transversion_weight
                )

            @staticmethod
            def weighted_mutation(original_codon, candidate_codons, transition_weight, transversion_weight):
                weights = []
                for candidate in candidate_codons:
                    if GeneUtils.is_transition(original_codon, candidate):
                        weights.append(transition_weight)
                    else:
                        weights.append(transversion_weight)
                return random.choices(candidate_codons, weights=weights, k=1)[0]

        class InDels:
            @staticmethod
            def insertion(sequence, num_insertions=1, position=None, is_codon=False, sequence_type="DNA"):
                valid_bases = {"DNA": "ATGC", "RNA": "AUGC"}[sequence_type]
                sequence = list(sequence)
                for _ in range(num_insertions):
                    insert_base = "".join(random.choices(valid_bases, k=3 if is_codon else 1))
                    insert_pos = position if position is not None else random.randint(0, len(sequence))
                    sequence.insert(insert_pos, insert_base)
                    MutationEngine.mutation_log.append({
                        "type": "insertion",
                        "position": insert_pos,
                        "inserted": insert_base,
                        "is_codon": is_codon
                    })
                return "".join(sequence)

            @staticmethod
            def deletion(sequence, num_deletions=1, start_position=None, is_codon=False):
                sequence = list(sequence)
                deletion_length = 3 if is_codon else 1
                for _ in range(num_deletions):
                    del_start = start_position if start_position is not None else random.randint(0, len(sequence) - deletion_length)
                    del_end = del_start + deletion_length
                    deleted_segment = sequence[del_start:del_end]
                    del sequence[del_start:del_end]
                    MutationEngine.mutation_log.append({
                        "type": "deletion",
                        "position": del_start,
                        "deleted": "".join(deleted_segment),
                        "is_codon": is_codon
                    })
                return "".join(sequence)

            @staticmethod
            def frameshift(sequence, shift_type="insertion", position=None, sequence_type="DNA"):
                valid_bases = {"DNA": "ATGC", "RNA": "AUGC"}[sequence_type]
                sequence = list(sequence)
                shift_pos = position if position is not None else random.randint(0, len(sequence))
                if shift_type == "insertion":
                    insert_base = random.choice(valid_bases)
                    sequence.insert(shift_pos, insert_base)
                    MutationEngine.mutation_log.append({
                        "type": "frameshift_insertion",
                        "position": shift_pos,
                        "inserted": insert_base
                    })
                elif shift_type == "deletion":
                    if len(sequence) > shift_pos:
                        deleted_base = sequence[shift_pos]
                        del sequence[shift_pos]
                        MutationEngine.mutation_log.append({
                            "type": "frameshift_deletion",
                            "position": shift_pos,
                            "deleted": deleted_base
                        })
                return "".join(sequence)

#     class ProteinMutations:
#         class Missense:
#             @staticmethod
#             def conservative():
#                 return None

#             @staticmethod
#             def non_conservative():
#                 return None

#         class Frameshift:
#             @staticmethod
#             def test():
#                 return None
