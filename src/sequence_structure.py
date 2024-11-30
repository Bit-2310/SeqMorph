from collections import defaultdict
from input_module import InputHandler

class Sequence_Structure:
    """
    A custom data structure to store and manage biological sequences.
    Supports storing sequences, metadata (including sequence length), sequence type count, and mutation operations.
    """
    def __init__(self, sequence_structure):
        """
        Initializes the Sequence_Structure object.
        :param sequence_structure: The sequence data to initialize the object.
        """
        self.seqs = sequence_structure
        self.Sequences = {}  # Hash map for storing sequences by chunk keys (Seq_ID_segment_1, Seq_ID_segment_2, ...)
        self.Sequence_type = defaultdict(int)  # A summary count of sequence types (DNA, RNA, Protein)
        self.MetaData = {}  # A dictionary to store metadata (e.g., GC content, Hydrophobicity, Length)

    def add_sequence(self, sequence_id, sequence, seq_type, chunk_size=999):
        """
        Adds a sequence to the data structure and divides it into smaller chunks for storage.
        :param sequence_id: The unique identifier of the sequence (e.g., NCBI or UniProt ID).
        :param sequence: The sequence string (DNA, RNA, or Protein).
        :param seq_type: The type of sequence ('DNA', 'RNA', or 'Protein').
        :param chunk_size: The size of each chunk (default is 999 base pairs for better alignment with k-mers).
        """
        self.Seq_ID = sequence_id
        # Divide the sequence into smaller chunks of size chunk_size (999)
        sequence_chunks = [sequence[i:i+chunk_size] for i in range(0, len(sequence), chunk_size)]
        
        # Store each chunk in the hash map with a unique key
        for index, chunk in enumerate(sequence_chunks):
            chunk_key = f"{sequence_id}_segment_{index+1}"
            self.Sequences[chunk_key] = chunk  # Store the chunk in the hash map
        
        # Update sequence type count
        self.Sequence_type[seq_type] += 1
        
        # Calculate and store metadata based on the sequence type
        self.update_metadata(sequence_id, sequence, seq_type)

    def update_metadata(self, sequence_id, sequence, seq_type=None):
        """
        Updates the metadata for a given sequence.
        :param sequence_id: The unique identifier of the sequence.
        :param sequence: The sequence string.
        :param seq_type: The type of sequence ('DNA', 'RNA', or 'Protein').
        """
        if seq_type is None:
            seq_type = self.MetaData.get(sequence_id, {}).get('Type')
        
        sequence_length = len(sequence)
        if seq_type in ['DNA', 'RNA']:
            self.MetaData[sequence_id] = {
                'GC_content': self.calculate_gc_content(sequence),
                'Length': sequence_length,  # Store the length of the sequence
                'Type': seq_type  # Store the type of the sequence
            }
        elif seq_type == 'Protein':
            self.MetaData[sequence_id] = {
                'Hydrophobicity': self.calculate_hydrophobicity(sequence),
                'Length': sequence_length,  # Store the length of the sequence
                'Type': seq_type  # Store the type of the sequence
            }

    def get_sequence(self, accession_id):
        """
        Retrieves the full sequence for a given accession_id by combining all chunks.
        :param accession_id: The accession ID of the sequence to retrieve.
        :return: The full sequence string if found, else None.
        """
        # Retrieve all chunks for the given accession_id
        full_sequence = ''.join([chunk for chunk_key, chunk in self.Sequences.items() if chunk_key.startswith(accession_id)])

        if not full_sequence:
            print(f"Error: Sequence with accession ID {accession_id} not found in Sequences.")
            return None
        
        return full_sequence

    def get_sequences(self):
        """
        Returns all stored sequences (chunks).
        :return: Dictionary with chunk keys and their corresponding sequences.
        """
        return self.Sequences

    def get_summary(self):
        """
        Returns the summary of the sequence types in the data structure.
        :return: Dictionary with counts of each sequence type (DNA, RNA, Protein).
        """
        return dict(self.Sequence_type)

    def get_sequence_length(self, sequence_id):
        """
        Returns the length of the sequence.
        :param sequence_id: The ID of the sequence to retrieve length for.
        :return: Length of the sequence.
        """
        return len(self.Sequences.get(sequence_id, ""))

    def calculate_gc_content(self, sequence):
        """
        Calculate the GC content of a DNA/RNA sequence.
        :param sequence: The DNA or RNA sequence string.
        :return: The GC content as a percentage.
        """
        gc_count = sum(1 for base in sequence.upper() if base in 'GC')
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    def get_full_sequence_length(self, sequence_id):
        """
        Retrieve the full sequence length from metadata.
        :param sequence_id: The ID of the sequence.
        :return: The length of the sequence.
        """
        return self.MetaData.get(sequence_id, {}).get('Length', 0)  # Use metadata length if available

    def calculate_hydrophobicity(self, sequence):
        """
        Calculate the hydrophobicity of a protein sequence.
        Simple hydrophobicity scale based on amino acid properties.
        :param sequence: The protein sequence string.
        :return: Hydrophobicity score.
        """
        hydrophobicity_scale = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
            'H': -0.5, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
            'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
            'W': -0.9, 'Y': -1.3
        }
        hydrophobicity_score = sum(hydrophobicity_scale.get(aa, 0) for aa in sequence.upper())
        return hydrophobicity_score / len(sequence) if sequence else 0