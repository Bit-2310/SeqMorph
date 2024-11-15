import random
from rna_library import RNAUtils
from protein_library import ProteinUtils

class MutationEngine:
    mutation_log = []  # Global log for tracking mutations

    @staticmethod
    def mutate(sequence, mutation_type, mutation_rate=10, start=None, end=None):
        if not sequence:
            raise ValueError("Input sequence is empty.")

        if mutation_rate < 0 or mutation_rate > 100:
            raise ValueError("Mutation rate must be between 0 and 100.")

        if start is not None and end is not None:
            if not (0 <= start < len(sequence)) or not (0 < end <= len(sequence)):
                raise ValueError("Invalid range specified. Ensure 0 <= start < end <= sequence length.")
            range_sequence = sequence[start:end]
        else:
            range_sequence = sequence

        mutation_count = int(len(range_sequence) * (mutation_rate / 100))
        mutated_sequence = MutationEngine._apply_mutation(range_sequence, mutation_type, mutation_count)
        
        # Integrate mutated region back into the full sequence if a range was specified
        if start is not None and end is not None:
            return sequence[:start] + mutated_sequence + sequence[end:]
        return mutated_sequence

    @staticmethod
    def mutate_protein_to_rna_and_back(protein_sequence, mutation_type, mutation_rate=10, start=None, end=None):
       # Step 1: Convert protein to RNA sequence using codon table from ProteinUtils
        rna_sequence = ProteinUtils.protein_to_rna(protein_sequence)

        # Step 2: Mutate the RNA sequence
        mutated_rna = MutationEngine.mutate(rna_sequence, mutation_type, mutation_rate, start, end)

        # Step 3: Convert back to protein sequence using RNAUtils translation
        mutated_protein = RNAUtils.translate(mutated_rna)

        return mutated_protein

    @staticmethod
    def _apply_mutation(sequence, mutation_type, mutation_count):
        sequence = list(sequence)
        if not sequence:
            raise ValueError("Sequence is empty or invalid.")

        positions = random.sample(range(len(sequence)), min(mutation_count, len(sequence)))

        for pos in positions:
            if mutation_type == 'substitute':
                original = sequence[pos]
                if sequence[pos] in "ATGC":
                    sequence[pos] = random.choice("ATGC".replace(sequence[pos], ""))
                elif sequence[pos] in "AUGC":
                    sequence[pos] = random.choice("AUGC".replace(sequence[pos], ""))
                elif sequence[pos] in "ACDEFGHIKLMNPQRSTVWY":
                    sequence[pos] = random.choice("ACDEFGHIKLMNPQRSTVWY".replace(sequence[pos], ""))
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": "substitute",
                    "original": original,
                    "mutated": sequence[pos]
                })
            elif mutation_type == 'insert':
                insertion = random.choice(sequence)
                sequence.insert(pos, insertion)
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": "insert",
                    "inserted": insertion
                })
            elif mutation_type == 'delete':
                original = sequence[pos]
                sequence[pos] = ''
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": "delete",
                    "original": original
                })
            elif mutation_type == 'frameshift':
                original = sequence[pos]
                sequence[pos] = ''  # Deleting to cause a shift
                MutationEngine.mutation_log.append({
                    "position": pos,
                    "type": "frameshift",
                    "original": original
                })
        return ''.join(sequence)

    @staticmethod
    def export_mutation(sequence, mutated_sequence, output_file):
        try:
            with open(output_file, "w") as f:
                f.write(">Original Sequence\n")
                f.write(sequence + "\n\n")
                f.write(">Mutated Sequence\n")
                f.write(mutated_sequence + "\n")
            print(f"Mutated sequence exported to {output_file}")
        except Exception as e:
            raise IOError(f"Failed to write to file {output_file}: {e}")

    @staticmethod
    def highlight_mutations(original_sequence, mutated_sequence):
        highlighted = []
        for o, m in zip(original_sequence, mutated_sequence):
            if o != m:
                highlighted.append(f"[{m}]")  # Highlight mutation
            else:
                highlighted.append(o)
        return ''.join(highlighted)

if __name__ == "__main__":
    sequence = "ATGCGTACGTTAGC"

    # Random substitutions across the entire sequence (10% mutated)
    print("Original Sequence:", sequence)
    mutated_sequence = MutationEngine.mutate(sequence, "substitute", mutation_rate=10)
    print("Substitution Mutation:", mutated_sequence)

    # Substitutions in a specific range (50% of the range mutated)
    range_mutated_sequence = MutationEngine.mutate(sequence, "substitute", mutation_rate=50, start=3, end=10)
    print("Range-Specific Substitution:", range_mutated_sequence)

    # Insertions across the entire sequence (5% mutated)
    inserted_sequence = MutationEngine.mutate(sequence, "insert", mutation_rate=5)
    print("Insertion Mutation:", inserted_sequence)

    # Deletions in a specific range (20% of the range mutated)
    deleted_sequence = MutationEngine.mutate(sequence, "delete", mutation_rate=20, start=5, end=12)
    print("Range-Specific Deletion:", deleted_sequence)

    # Frameshift mutation example
    frameshift_sequence = MutationEngine.mutate(sequence, "frameshift", mutation_rate=10)
    print("Frameshift Mutation:", frameshift_sequence)

    # Highlighting mutations
    highlighted_sequence = MutationEngine.highlight_mutations(sequence, mutated_sequence)
    print("Highlighted Mutations:", highlighted_sequence)

    # Export the mutated sequence to a file
    MutationEngine.export_mutation(sequence, mutated_sequence, "mutated_output.txt")

    # Protein to RNA and back mutation example
    protein_sequence = "ACDEFGHIK"
    mutated_protein = MutationEngine.mutate_protein_to_rna_and_back(protein_sequence, "substitute", mutation_rate=20)
    print("Original Protein Sequence:", protein_sequence)
    print("Mutated Protein Sequence:", mutated_protein)

    # Print mutation log
    print("Mutation Log:", MutationEngine.mutation_log)
