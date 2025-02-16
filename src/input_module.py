import os
import re
import json
import requests
from Bio import SeqIO
from tkinter import Tk, filedialog
from gene_library import SequenceValidation
import logging
from typing import Dict

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
    BASE_URLS = {
        "NCBI": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id={id}&report=fasta",
        "UniProt": "https://www.uniprot.org/uniprot/{id}.fasta"
    }
    HEADERS = {"User-Agent": "SeqMorphTool/1.0 (+https://github.com/Bit_2310/SeqMorph)"}

    def __init__(self, cache_file="sequence_cache.json", timeout=10):
        self.cache_file = cache_file
        self.timeout = timeout
        self.sequence_map = self._load_cache()

        # Initialize cache file if it doesn't exist or is invalid
        if not os.path.exists(self.cache_file) or not self.sequence_map:
            self._save_cache()

    def _load_cache(self) -> Dict:
        if not os.path.exists(self.cache_file):
            logging.info(f"Cache file {self.cache_file} not found. Initializing empty cache.")
            return {}
        try:
            with open(self.cache_file, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logging.warning(f"Invalid or empty cache file {self.cache_file}. Initializing empty cache. Error: {e}")
            return {}

    def _save_cache(self):
        with open(self.cache_file, "w") as f:
            json.dump(self.sequence_map, f)

    def fetch_sequence(self, accession_id):
        if accession_id in self.sequence_map:
            return self.sequence_map[accession_id]
        else:
            seq = self._fetch_from_source(accession_id)
            self.sequence_map[accession_id] = seq
            self._save_cache()
            return seq

    def _fetch_from_source(self, accession_id):
        if is_ncbi_accession(accession_id):
            return self.fetch_from_ncbi(accession_id)
        elif is_uniprot_accession(accession_id):
            return self.fetch_from_uniprot(accession_id)
        else:
            raise ValueError(f"Invalid accession ID format: {accession_id}")

    def fetch_from_ncbi(self, accession_id):
        url = self.BASE_URLS["NCBI"].format(id=accession_id)
        return self._fetch_sequence(url, "NCBI")

    def fetch_from_uniprot(self, accession_id):
        url = self.BASE_URLS["UniProt"].format(id=accession_id)
        return self._fetch_sequence(url, "UniProt")

    def _fetch_sequence(self, url, source):
        try:
            response = requests.get(url, headers=self.HEADERS, timeout=self.timeout)
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

def is_ncbi_accession(accession_id):
    return bool(re.match(r"^(AC_|NC_|NG_|NM_|NP_|NR_|NT_|NW_|XM_|XP_|XR_|YP_|ZP_)[0-9]+", accession_id))

def is_uniprot_accession(accession_id):
    return bool(re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", accession_id))