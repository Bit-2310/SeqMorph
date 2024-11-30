# input_module.py : Handles importing and storing sequences from files or online sources.
import os
import re
import requests
from Bio import SeqIO
from tkinter import Tk, filedialog
from gene_library import SequenceValidation

class InputHandler:
    """
    Handles user input for sequences through file uploads or direct entry, 
    using a generator to handle large sequences and a hashmap for storage.
    """
    SUPPORTED_FILETYPES = [
        ("FASTA Files", "*.fasta"),
        ("Text Files", "*.txt"),
        ("All Files", "*.*")
    ]

    def select_file(self):
        """
        Open a file dialog to select a single file.
        :return: Path of the selected file.
        """
        root = Tk()
        root.withdraw()  # Hide the root window
        return filedialog.askopenfilename(
            title="Select a Sequence File",
            filetypes=self.SUPPORTED_FILETYPES
        )

    def select_multiple_files(self):
        """
        Open a file dialog to select multiple files.
        :return: List of paths for the selected files.
        """
        root = Tk()
        root.withdraw()  # Hide the root window
        return filedialog.askopenfilenames(
            title="Select Sequence Files",
            filetypes=self.SUPPORTED_FILETYPES
        )

    def load_sequence_file(self, filepath):
        """
        Load and parse a sequence file using a generator.
        Supports `.fasta` and `.txt` file formats.
        :param filepath: Path to the sequence file.
        :return: A generator yielding sequences from the file.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        file_extension = os.path.splitext(filepath)[1].lower()
        if file_extension not in [".fasta", ".txt"]:
            raise ValueError(f"Unsupported file format: {file_extension}")

        try:
            if file_extension == ".fasta":
                return self._parse_fasta(filepath)
            elif file_extension == ".txt":
                return self._parse_text(filepath)
        except Exception as e:
            raise ValueError(f"Error reading file {filepath}: {str(e)}")

    def _parse_fasta(self, filepath):
        """
        Generator to parse a FASTA file to yield sequences.
        :param filepath: Path to the FASTA file.
        :return: A generator yielding sequence strings.
        """
        for record in SeqIO.parse(filepath, "fasta"):
            yield str(record.seq)

    def _parse_text(self, filepath):
        """
        Parse a plain text file to extract a single sequence.
        :param filepath: Path to the text file.
        :return: A list containing a single sequence.
        """
        with open(filepath, "r") as file:
            sequence = file.read().strip()
            if not sequence:
                raise ValueError("No content found in the TXT file.")
            yield sequence  # Returning a generator

    def sequence_type(self, sequence):
        """
        Determine the type of sequence (DNA, RNA, or Protein) and validate it.
        :param sequence: Input sequence string.
        :return: The type of sequence ('DNA', 'RNA', 'Protein').
        """
        sequence = sequence.upper()
        try:
            # Use SequenceValidation methods for detection and validation
            seq_type = SequenceValidation.detect_seq(sequence)
            if not SequenceValidation.validate(sequence, seq_type):
                raise ValueError(f"Invalid {seq_type} sequence.")
            return seq_type
        except ValueError as e:
            raise ValueError(f"Sequence validation failed: {e}")
class SequenceFetcher:
    """
    Fetch sequences from online resources like NCBI and UniProt.
    Sequences are stored in a hash map for efficient retrieval.
    """
    BASE_URLS = {
        "NCBI": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id={id}&report=fasta",
        "UniProt": "https://www.uniprot.org/uniprot/{id}.fasta"
    }
    HEADERS = {"User-Agent": "SeqMorphTool/1.0 (+https://github.com/Bit_2310/SeqMorph)"}

    def __init__(self):
        # Initialize an empty hash map (dictionary) to store sequences by their accession ID
        self.sequence_map = {}

    def fetch_sequence(self, accession_id):
        """
        Fetch a sequence from the appropriate source based on the accession ID.
        Automatically detects whether it's from NCBI or UniProt.
        :param accession_id: The accession ID (NCBI or UniProt).
        :return: The sequence string.
        """
        if accession_id in self.sequence_map:
            # If the sequence is already stored in the map, return it directly
            return self.sequence_map[accession_id]
        else:
            # Otherwise, fetch the sequence from the source
            if is_ncbi_accession(accession_id):
                seq = self.fetch_from_ncbi(accession_id)
            elif is_uniprot_accession(accession_id):
                seq = self.fetch_from_uniprot(accession_id)
            else:
                raise ValueError(f"Invalid accession ID format: {accession_id}")
            
            # Store the fetched sequence in the hash map for future use
            self.sequence_map[accession_id] = seq
            return seq

    @staticmethod
    def fetch_from_ncbi(accession_id):
        """
        Fetch a sequence from NCBI using the accession ID.
        :param accession_id: NCBI accession ID.
        :return: The sequence string.
        """
        url = SequenceFetcher.BASE_URLS["NCBI"].format(id=accession_id)
        return SequenceFetcher._fetch_sequence(url, "NCBI")

    @staticmethod
    def fetch_from_uniprot(accession_id):
        """
        Fetch a sequence from UniProt using the accession ID.
        :param accession_id: UniProt ID.
        :return: The sequence string.
        """
        url = SequenceFetcher.BASE_URLS["UniProt"].format(id=accession_id)
        return SequenceFetcher._fetch_sequence(url, "UniProt")

    @staticmethod
    def _fetch_sequence(url, source):
        """
        Internal method to fetch and validate a sequence from a given URL.
        :param url: URL to fetch the sequence.
        :param source: Source name (e.g., 'NCBI', 'UniProt').
        :return: The sequence string.
        """
        try:
            response = requests.get(url, headers=SequenceFetcher.HEADERS, timeout=10)
            response.raise_for_status()

            lines = response.text.splitlines()
            if not lines or not lines[0].startswith(">"):
                raise ValueError(f"Invalid FASTA format received from {source}.")
            
            sequence_data = "".join(lines[1:]).strip()
            if not sequence_data:
                raise ValueError(f"No sequence data found from {source}.")
            return sequence_data
        except requests.exceptions.Timeout:
            raise ValueError(f"Request to {source} timed out.")
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Error fetching sequence from {source}: {str(e)}")

# Helper functions for accession ID detection
def is_ncbi_accession(accession_id):
    """Check if an accession ID matches NCBI's pattern."""
    return bool(re.match(r"^(NC_|NM_|NW_)[0-9]+", accession_id))

def is_uniprot_accession(accession_id):
    """Check if an accession ID matches UniProt's pattern."""
    return bool(re.match(r"^[A-Za-z0-9]{6,}", accession_id))  # Adjust as per UniProt ID length