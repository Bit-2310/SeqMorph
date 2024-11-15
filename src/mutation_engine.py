import random
from rna_library import RNAUtils
from protein_library import ProteinUtils
from dna_library import DNAUtils


class MutationEngine:
    mutation_log = []  # Global log for tracking mutations

    @staticmethod
    def mutate(sequence, mutation_type, mutation_rate=10, start=None, end=None):
        """
        Mutate a given sequence based on mutation type and rate.

        Args:
            sequence (str): The input sequence (DNA, RNA, or Protein).
            mutation_type (str): The type of mutation ('substitute', 'insert', 'delete', 'frameshift').
            mutation_rate (float): Mutation rate as a percentage (0-100).
            start (int | list[int], optional): Start position(s) for mutation range.
            end (int | list[int], optional): End position(s) for mutation range.

        Returns:
            str: Mutated sequence.
        """
        if not sequence:
            raise ValueError("Input sequence is empty.")

        if mutation_rate < 0 or mutation_rate > 100:
            raise ValueError("Mutation rate must be between 0 and 100.")

        range_sequence = MutationEngine._extract_range(sequence, start, end)

        mutation_count = int(len(range_sequence) * (mutation_rate / 100))
        sequence_type = MutationEngine._detect_sequence_type(sequence)

        mutated_sequence = MutationEngine._apply_mutation(
            range_sequence, mutation_type, mutation_count, mutation_rate, sequence_type, start, end
        )

        # Integrate mutated region back into the full sequence if a range was specified
        if start is not None and end is not None:
            return MutationEngine._integrate_mutated_range(sequence, mutated_sequence, start, end)
        return mutated_sequence

    @staticmethod
    def highlight_mutations(original_sequence, mutated_sequence):
        """
        Highlight mutations in the mutated sequence compared to the original sequence.

        Args:
            original_sequence (str): The original sequence before mutation.
            mutated_sequence (str): The sequence after mutation.

        Returns:
            str: Highlighted sequence with mutations.
        """
        highlighted = []
        for o, m in zip(original_sequence, mutated_sequence):
            if o != m:
                highlighted.append(f"[{m}]")  # Highlight mutation
            else:
                highlighted.append(o)
        return ''.join(highlighted)

    @staticmethod
    def _extract_range(sequence, start, end):
        """Extracts a range or multiple ranges from a sequence."""
        if isinstance(start, list) and isinstance(end, list):
            if len(start) != len(end):
                raise ValueError("Start and end lists must have the same length.")
            return "".join(sequence[s:e] for s, e in zip(start, end) if 0 <= s < e <= len(sequence))
        elif start is not None and end is not None:
            if not (0 <= start < len(sequence)) or not (0 < end <= len(sequence)):
                raise ValueError("Invalid range specified. Ensure 0 <= start < end <= sequence length.")
            return sequence[start:end]
        return sequence

    @staticmethod
    def _integrate_mutated_range(sequence, mutated_sequence, start, end):
        """Integrates a mutated range back into the full sequence."""
        if isinstance(start, list) and isinstance(end, list):
            for s, e, sub_mutation in zip(start, end, mutated_sequence.split('|')):
                sequence = sequence[:s] + sub_mutation + sequence[e:]
            return sequence
        return sequence[:start] + mutated_sequence + sequence[end:]

    @staticmethod
    def _detect_sequence_type(sequence):
        """
        Detect the type of a sequence (DNA, RNA, Protein).
        """
        if set(sequence.upper()).issubset("ATGC"):
            return "DNA"
        elif set(sequence.upper()).issubset("AUGC"):
            return "RNA"
        elif set(sequence.upper()).issubset("ACDEFGHIKLMNPQRSTVWY"):
            return "Protein"
        else:
            raise ValueError("Invalid sequence type.")

    @staticmethod
    def _apply_mutation(sequence, mutation_type, mutation_count, mutation_rate, sequence_type, start, end):
        """
        Apply mutations to a sequence.

        Args:
            sequence (str): Sequence to mutate.
            mutation_type (str): Type of mutation.
            mutation_count (int): Number of mutations to apply.
            mutation_rate (float): Mutation rate as a percentage.
            sequence_type (str): Type of the sequence (DNA, RNA, Protein).
            start (int | list[int], optional): Start position(s) for mutation range.
            end (int | list[int], optional): End position(s) for mutation range.

        Returns:
            str: Mutated sequence.
        """
        sequence = list(sequence)
        if not sequence:
            raise ValueError("Sequence is empty or invalid.")

        positions = random.sample(range(len(sequence)), min(mutation_count, len(sequence)))

        for pos in positions:
            if mutation_type == 'substitute':
                original = sequence[pos]
                sequence[pos] = MutationEngine._substitute_mutation(sequence[pos], sequence_type)
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": mutation_type,
                    "original": original,
                    "mutated": sequence[pos],
                    "rate": mutation_rate,
                    "range": f"{start}-{end}" if start is not None and end is not None else "Full sequence",
                    "sequence_type": sequence_type
                })
            elif mutation_type == 'insert':
                insertion = MutationEngine._random_insert(sequence_type)
                sequence.insert(pos, insertion)
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": mutation_type,
                    "inserted": insertion
                })
            elif mutation_type == 'delete':
                original = sequence[pos]
                sequence[pos] = ''
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": mutation_type,
                    "original": original
                })
            elif mutation_type == 'frameshift':
                original = sequence[pos]
                sequence[pos] = ''  # Deleting to cause a shift
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": mutation_type,
                    "original": original
                })
        return ''.join(sequence)

    @staticmethod
    def _substitute_mutation(base, sequence_type):
        """Generate a substitution mutation."""
        if sequence_type == "DNA":
            return random.choice("ATGC".replace(base, ""))
        elif sequence_type == "RNA":
            return random.choice("AUGC".replace(base, ""))
        elif sequence_type == "Protein":
            return random.choice("ACDEFGHIKLMNPQRSTVWY".replace(base, ""))
        raise ValueError("Unknown sequence type.")

    @staticmethod
    def _random_insert(sequence_type):
        """Generate a random insertion for the given sequence type."""
        if sequence_type == "DNA":
            return random.choice("ATGC")
        elif sequence_type == "RNA":
            return random.choice("AUGC")
        elif sequence_type == "Protein":
            return random.choice("ACDEFGHIKLMNPQRSTVWY")
        raise ValueError("Unknown sequence type.")
