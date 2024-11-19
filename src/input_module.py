import os
from Bio import SeqIO
from tkinter import Tk, filedialog
import requests
from gene_library import SequenceValidation


class InputHandler:
    """
    Handles user input for sequences through file uploads or direct entry.
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
        Load and parse a sequence file.
        Supports `.fasta` and `.txt` file formats.
        :param filepath: Path to the sequence file.
        :return: List of sequences from the file.
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

    @staticmethod
    def _parse_fasta(filepath):
        """
        Parse a FASTA file to extract sequences.
        :param filepath: Path to the FASTA file.
        :return: List of sequences from the file.
        """
        sequences = [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]
        if not sequences:
            raise ValueError("No sequences found in the FASTA file.")
        return sequences

    @staticmethod
    def _parse_text(filepath):
        """
        Parse a plain text file to extract a sequence.
        :param filepath: Path to the text file.
        :return: List containing a single sequence.
        """
        with open(filepath, "r") as file:
            sequence = file.read().strip()
            if not sequence:
                raise ValueError("No content found in the TXT file.")
            return [sequence]
        
    @staticmethod
    def sequence_type(sequence):
        """
        Determine the type of sequence (DNA, RNA, or Protein) and validate it.
        :param sequence: Input sequence string.
        :return: The type of sequence ('DNA', 'RNA', 'Protein').
        """
        sequence = sequence.upper()
        try:
            # Use SequenceValidation methods for detection and validation
            sequence_type = SequenceValidation.detect_sequence_type(sequence)
            if not SequenceValidation.validate(sequence, sequence_type):
                raise ValueError(f"Invalid {sequence_type} sequence.")
            return sequence_type
        except ValueError as e:
            raise ValueError(f"Sequence validation failed: {e}")
        
class SequenceFetcher:
    """
    Fetch sequences from online resources like NCBI and UniProt.
    """
    BASE_URLS = {
        "NCBI": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id={id}&report=fasta",
        "UniProt": "https://www.uniprot.org/uniprot/{id}.fasta"
    }
    HEADERS = {"User-Agent": "SeqMorphTool/1.0 (+https://github.com/your-repo)"}

    @staticmethod
    def fetch_sequence_ncbi(accession_id):
        """
        Fetch a sequence from NCBI using the accession ID.
        :param accession_id: NCBI accession ID.
        :return: The sequence string.
        """
        url = SequenceFetcher.BASE_URLS["NCBI"].format(id=accession_id)
        return SequenceFetcher._fetch_sequence(url, "NCBI")

    @staticmethod
    def fetch_sequence_uniprot(uniprot_id):
        """
        Fetch a sequence from UniProt using the UniProt ID.
        :param uniprot_id: UniProt ID.
        :return: The sequence string.
        """
        url = SequenceFetcher.BASE_URLS["UniProt"].format(id=uniprot_id)
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