import random
from protein_library import ProteinUtils
from gene_library import GeneUtils
from analysis import SequenceAnalysisReport

class MutationUtils:
    """
    A utility class for handling mutations on sequences (DNA, RNA, Protein).
    """
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.sequence_type = GeneUtils.detect_sequence_type(sequence)
        self.mutation_log = []  # Log of all mutations applied to this sequence

    def validate_sequence(self):
        """
        Validate the sequence using utility libraries.
        """
        if self.sequence_type in ["DNA", "RNA"]:
            return GeneUtils.validate(self.sequence)
        elif self.sequence_type == "Protein":
            return ProteinUtils.validate(self.sequence)
        raise ValueError(f"Unknown sequence type: {self.sequence_type}")

    def apply_random_mutation(self, mutation_rate):
        """
        Apply random mutations to the sequence.
        """
        if self.sequence_type in ["DNA", "RNA"]:
            mutated_sequence = GeneUtils.random_mutation(self.sequence, mutation_rate)
        elif self.sequence_type == "Protein":
            mutated_sequence = ProteinUtils.random_mutation(self.sequence, mutation_rate)
        else:
            raise ValueError(f"Unknown sequence type: {self.sequence_type}")

        # Log mutation details along with sequence analysis
        self._log_mutation("Random Mutation", None, self.sequence, mutated_sequence, mutation_rate)
        self.sequence = mutated_sequence
        return mutated_sequence

    def _log_mutation(self, mutation_type, subtype, original_sequence, mutated_sequence, mutation_rate=None, position=None):
        """
        Log details of a mutation and include analysis metrics for both sequences.
        """
        # Create analysis report for both original and mutated sequences
        analysis_report = self._generate_analysis_report(original_sequence, mutated_sequence)

        # Log mutation and analysis results
        mutation_entry = {
            "mutation_type": mutation_type,
            "subtype": subtype,
            "original_sequence": original_sequence,
            "mutated_sequence": mutated_sequence,
            "mutation_rate": mutation_rate,
            "position": position,
            "analysis_report": analysis_report  # Include analysis results
        }

        self.mutation_log.append(mutation_entry)

    def _generate_analysis_report(self, original_sequence, mutated_sequence):
        """
        Generate an analysis report comparing the original and mutated sequences.
        """
        original_analysis = SequenceAnalysisReport(original_sequence, self.sequence_type)
        mutated_analysis = SequenceAnalysisReport(mutated_sequence, self.sequence_type)

        # Compare the sequences and generate the comparison report
        comparison_report = original_analysis.compare_with(mutated_analysis)
        
        return comparison_report

    def get_mutation_log(self):
        """
        Retrieve mutation logs.
        """
        return self.mutation_log
class MutationEngine:
    mutation_log = []  # Global mutation log

    class GeneMutations:
        @staticmethod
        def weighted_mutation(original_codon, candidate_codons, transition_weight, transversion_weight):
            """
            Perform weighted mutation using transition and transversion probabilities.
            """
            weights = [
                transition_weight if GeneUtils.is_transition(original_codon, candidate) else transversion_weight
                for candidate in candidate_codons
            ]
            return random.choices(candidate_codons, weights=weights, k=1)[0]

        class PointMutation:
            @staticmethod
            def silent(codon, transition_weight=0.7, transversion_weight=0.3):
                """
                Perform silent mutation.
                """
                amino_acid = GeneUtils.CODON_TABLE_DNA.get(codon)
                if not amino_acid:
                    raise ValueError(f"Invalid codon: {codon}")

                synonymous_codons = [
                    key for key, aa in GeneUtils.CODON_TABLE_DNA.items()
                    if aa == amino_acid and key != codon
                ]
                if not synonymous_codons:
                    return codon  # No synonymous mutation possible

                return MutationEngine.GeneMutations.weighted_mutation(
                    codon, synonymous_codons, transition_weight, transversion_weight
                )

            @staticmethod
            def missense(codon, transition_weight=0.7, transversion_weight=0.3):
                """
                Perform missense mutation.
                """
                amino_acid = GeneUtils.CODON_TABLE_DNA.get(codon)
                if not amino_acid:
                    raise ValueError(f"Invalid codon: {codon}")

                nonsynonymous_codons = [
                    key for key, aa in GeneUtils.CODON_TABLE_DNA.items()
                    if aa != amino_acid and aa != 'STOP'
                ]
                if not nonsynonymous_codons:
                    return codon  # No missense mutation possible

                return MutationEngine.GeneMutations.weighted_mutation(
                    codon, nonsynonymous_codons, transition_weight, transversion_weight
                )

            @staticmethod
            def nonsense(codon, transition_weight=0.7, transversion_weight=0.3):
                """
                Perform nonsense mutation to introduce a stop codon.
                """
                if codon in GeneUtils.STOP_CODONS_DNA:
                    return codon  # Codon is already a stop codon

                possible_stop_codons = list(GeneUtils.STOP_CODONS_DNA - {codon})
                if not possible_stop_codons:
                    return codon  # No nonsense mutation possible

                return MutationEngine.GeneMutations.weighted_mutation(
                    codon, possible_stop_codons, transition_weight, transversion_weight
                )

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

    class ProteinMutations:
        class Missense:
            @staticmethod
            def conservative(protein_sequence, position):
                """
                Perform a conservative missense mutation.
                """
                amino_acid = protein_sequence[position]
                conservative_aa = ProteinUtils.get_conservative_group(amino_acid)
                if not conservative_aa:
                    return protein_sequence  # No conservative group available

                new_amino_acid = random.choice([aa for aa in conservative_aa if aa != amino_acid])
                return protein_sequence[:position] + new_amino_acid + protein_sequence[position + 1:]

            @staticmethod
            def non_conservative(protein_sequence, position):
                """
                Perform a non-conservative missense mutation.
                """
                amino_acid = protein_sequence[position]
                all_amino_acids = ProteinUtils.ALL_AMINO_ACIDS
                new_amino_acid = random.choice([aa for aa in all_amino_acids if aa != amino_acid])
                return protein_sequence[:position] + new_amino_acid + protein_sequence[position + 1:]