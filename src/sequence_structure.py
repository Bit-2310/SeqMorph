import os
from typing import Dict, Generator, Optional, List  # Added List import
from gene_library import SequenceValidation
from Bio import SeqIO  # Ensure SeqIO is imported for FASTA parsing

class SequenceStructure:
    """
    Manages sequences and their metadata, ensuring efficient storage and retrieval.
    """
    def __init__(self):
        self.sequences = {}  # Stores sequences by accession ID
        self.metadata = {}   # Stores metadata by accession ID

    def add_sequence(self, accession_id: str, sequence: str, 
                     seq_type: Optional[str] = None, **metadata):
        """
        Add a sequence with optional metadata.
        :param accession_id: Unique identifier for the sequence.
        :param sequence: The sequence data (DNA, RNA, or protein).
        :param seq_type: Sequence type (DNA, RNA, or protein).
        :param metadata: Additional metadata (e.g., source, description).
        """
        if not accession_id:
            raise ValueError("Accession ID cannot be empty.")
        if accession_id in self.sequences:
            raise ValueError(f"Sequence with ID {accession_id} already exists.")

        # Validate sequence type
        if not seq_type:
            seq_type = SequenceValidation.detect_seq(sequence)
        if not SequenceValidation.validate(sequence, seq_type):
            raise ValueError(f"Invalid {seq_type} sequence.")

        # Store sequence and metadata
        self.sequences[accession_id] = sequence
        self.metadata[accession_id] = {
            "type": seq_type,
            **metadata
        }

    def get_sequence(self, accession_id: str) -> str:
        """
        Retrieve a sequence by its accession ID.
        :param accession_id: The ID of the sequence to retrieve.
        :return: The sequence string.
        """
        if accession_id not in self.sequences:
            raise KeyError(f"Sequence {accession_id} not found.")
        return self.sequences[accession_id]

    def get_metadata(self, accession_id: str) -> Dict:
        """
        Retrieve metadata for a sequence by its accession ID.
        :param accession_id: The ID of the sequence.
        :return: A dictionary of metadata.
        """
        if accession_id not in self.metadata:
            raise KeyError(f"Metadata for {accession_id} not found.")
        return self.metadata[accession_id]

    def update_metadata(self, accession_id: str, **metadata):
        """
        Update metadata for a sequence.
        :param accession_id: The ID of the sequence.
        :param metadata: Key-value pairs to update.
        """
        if accession_id not in self.metadata:
            raise KeyError(f"Sequence {accession_id} not found.")
        self.metadata[accession_id].update(metadata)

    def remove_sequence(self, accession_id: str):
        """
        Remove a sequence and its metadata.
        :param accession_id: The ID of the sequence to remove.
        """
        if accession_id not in self.sequences:
            raise KeyError(f"Sequence {accession_id} not found.")
        del self.sequences[accession_id]
        del self.metadata[accession_id]

    def load_from_file(self, filepath: str, seq_type: Optional[str] = None):
        """
        Load sequences from a file (FASTA or TXT).
        :param filepath: Path to the file.
        :param seq_type: Optional sequence type (DNA, RNA, or protein).
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        with open(filepath, "r") as file:
            if filepath.endswith(".fasta"):
                self._load_fasta(file, seq_type)
            elif filepath.endswith(".txt"):
                self._load_txt(file, seq_type)
            else:
                raise ValueError("Unsupported file format.")

    def _load_fasta(self, file, seq_type: Optional[str] = None):
        """
        Load sequences from a FASTA file.
        """
        for record in SeqIO.parse(file, "fasta"):
            self.add_sequence(record.id, str(record.seq), seq_type)

    def _load_txt(self, file, seq_type: Optional[str] = None):
        """
        Load sequences from a plain text file.
        """
        sequence = file.read().strip()
        if not sequence:
            raise ValueError("No sequence found in the file.")
        self.add_sequence("default_id", sequence, seq_type)

    def export_to_file(self, filepath: str, accession_ids: Optional[List[str]] = None):
        """
        Export sequences to a file (FASTA format).
        :param filepath: Path to the output file.
        :param accession_ids: Optional list of IDs to export.
        """
        if not accession_ids:
            accession_ids = list(self.sequences.keys())

        with open(filepath, "w") as file:
            for acc_id in accession_ids:
                seq = self.sequences[acc_id]
                metadata = self.metadata[acc_id]
                file.write(f">{acc_id} {metadata.get('description', '')}\n")
                file.write(f"{seq}\n")

    def __len__(self):
        """Return the number of sequences stored."""
        return len(self.sequences)

    def __contains__(self, accession_id: str):
        """Check if a sequence exists by its accession ID."""
        return accession_id in self.sequences

    def __iter__(self):
        """Iterate over sequences and metadata."""
        for acc_id in self.sequences:
            yield acc_id, self.sequences[acc_id], self.metadata[acc_id]