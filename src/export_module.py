from gene_library import GeneUtils
from protein_library import ProteinUtils
from input_module import InputHandler
class ExportHandler:
    @staticmethod
    def export_sequence(sequence, mutated_sequence, output_file, sequence_name="Sequence"):
        """
        Export the original and mutated sequences along with their attributes to a file.

        Args:
            sequence (str): Original sequence.
            mutated_sequence (str): Mutated sequence.
            output_file (str): File name to save the output.
            sequence_name (str): Name or identifier for the sequence.
        """
        try:
            # Detect sequence type
            sequence_type = InputHandler._detect_sequence_type(sequence)

            # Attributes
            gc_content = None
            nucleotide_frequency = None
            hydrophobicity = None

            if sequence_type in ["DNA", "RNA"]:
                gc_content = GeneUtils.gc_content(sequence) if sequence_type == "RNA" else GeneUtils.gc_content(sequence)
                nucleotide_frequency = GeneUtils.nucleotide_frequency(sequence) if sequence_type == "DNA" else None
            elif sequence_type == "Protein":
                hydrophobicity = ProteinUtils.hydrophobicity(sequence)

            # Write to file
            with open(output_file, "w") as f:
                f.write(f"> Sequence Name: {sequence_name}\n")
                if gc_content is not None:
                    f.write(f"GC Content: {gc_content:.2f}%\n")
                if nucleotide_frequency is not None:
                    f.write(f"Nucleotide Frequency: {nucleotide_frequency}\n")
                if hydrophobicity is not None:
                    f.write(f"Hydrophobicity: {hydrophobicity:.2f}\n")
                f.write("\nOriginal Sequence:\n")
                f.write(sequence + "\n")
                f.write("\nMutated Sequence:\n")
                f.write(mutated_sequence + "\n")

            print(f"Mutated sequence exported to {output_file}")
        except Exception as e:
            raise IOError(f"Failed to write to file {output_file}: {e}")