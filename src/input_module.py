import os
import re
import json
import logging
import requests
from Bio import SeqIO
from typing import Dict, Generator

# Optional better GUI alternative
try:
    import tkinter as tk
    from tkinter import filedialog
except ImportError:
    tk = None  # fallback if GUI is not supported

try:
    import pyperclip
except ImportError:
    pyperclip = None

from library import SequenceValidation


class InputHandler:
    """
    Handles sequence input from files, clipboard, or text, using generator parsing and validation.
    """
    SUPPORTED_EXTENSIONS = {".fasta", ".fa", ".txt"}

    def select_file(self):
        if not tk:
            raise RuntimeError("Tkinter GUI is not supported in this environment.")
        root = tk.Tk()
        root.withdraw()
        return filedialog.askopenfilename(filetypes=[("Sequence Files", "*.fasta *.fa *.txt")])

    def select_multiple_files(self):
        if not tk:
            raise RuntimeError("Tkinter GUI is not supported in this environment.")
        root = tk.Tk()
        root.withdraw()
        return filedialog.askopenfilenames(filetypes=[("Sequence Files", "*.fasta *.fa *.txt")])

    def load_sequence_file(self, filepath: str) -> Generator[str, None, None]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        if not os.access(filepath, os.R_OK):
            raise PermissionError(f"Cannot read file: {filepath}")

        ext = os.path.splitext(filepath)[1].lower()
        if ext not in self.SUPPORTED_EXTENSIONS:
            raise ValueError(f"Unsupported file format: {ext}")

        try:
            if ext in [".fasta", ".fa"]:
                return self._parse_fasta(filepath)
            return self._parse_text(filepath)
        except Exception as e:
            raise ValueError(f"Error while reading file: {e}")

    def _parse_fasta(self, filepath: str) -> Generator[str, None, None]:
        records = list(SeqIO.parse(filepath, "fasta"))
        if not records:
            raise ValueError("No sequences found in FASTA file.")
        for record in records:
            sequence = str(record.seq).replace("\n", "").replace(" ", "")
            if len(sequence) < 10:
                raise ValueError(f"Sequence in {record.id} is too short to be useful.")
            yield sequence

    def _parse_text(self, filepath: str) -> Generator[str, None, None]:
        with open(filepath, "r") as f:
            content = f.read().strip()
            if not content:
                raise ValueError("No content found in the text file.")
            cleaned = content.replace("\n", "").replace(" ", "")
            if len(cleaned) < 10:
                raise ValueError("Sequence from text is too short.")
            yield cleaned

    def read_from_clipboard(self) -> str:
        if not pyperclip:
            raise RuntimeError("pyperclip not installed. Install it via pip to use clipboard support.")
        text = pyperclip.paste().strip()
        if not text:
            raise ValueError("Clipboard is empty.")
        cleaned = text.replace("\n", "").replace(" ", "")
        if not re.match(r"^[ACGTURYKMSWBDHVNacgturykmswbdhvnXx]+$", cleaned):
            raise ValueError("Clipboard content does not appear to be a valid biological sequence.")
        return cleaned

    def sequence_type(self, sequence: str) -> str:
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        if len(sequence) < 10:
            raise ValueError("Sequence is too short to classify.")
        seq_type = SequenceValidation.detect_seq(sequence)
        if not SequenceValidation.validate(sequence, seq_type):
            raise ValueError(f"Invalid {seq_type} sequence.")
        return seq_type


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
        if not self.sequence_map:
            self._save_cache()

    def _load_cache(self) -> Dict:
        if not os.path.exists(self.cache_file):
            return {}
        try:
            with open(self.cache_file, "r") as f:
                return json.load(f)
        except Exception:
            return {}

    def _save_cache(self):
        with open(self.cache_file, "w") as f:
            json.dump(self.sequence_map, f, indent=2)

    def fetch_sequence(self, accession_id: str) -> str:
        if accession_id in self.sequence_map:
            return self.sequence_map[accession_id]
        sequence = self._fetch_from_source(accession_id)
        if len(sequence) < 10:
            raise ValueError("Fetched sequence is unusually short. Please verify accession ID.")
        self.sequence_map[accession_id] = sequence
        self._save_cache()
        return sequence

    def _fetch_from_source(self, accession_id: str) -> str:
        if is_ncbi_accession(accession_id):
            return self.fetch_from_ncbi(accession_id)
        elif is_uniprot_accession(accession_id):
            return self.fetch_from_uniprot(accession_id)
        raise ValueError(f"Unrecognized accession ID format: {accession_id}")

    def fetch_from_ncbi(self, accession_id: str) -> str:
        url = self.BASE_URLS["NCBI"].format(id=accession_id)
        return self._fetch_sequence(url, "NCBI")

    def fetch_from_uniprot(self, accession_id: str) -> str:
        url = self.BASE_URLS["UniProt"].format(id=accession_id)
        return self._fetch_sequence(url, "UniProt")

    def _fetch_sequence(self, url: str, source: str) -> str:
        try:
            response = requests.get(url, headers=self.HEADERS, timeout=self.timeout)
            response.raise_for_status()
            lines = response.text.splitlines()
            if not lines or not lines[0].startswith(">"):
                raise ValueError(f"Invalid FASTA format from {source}.")
            return "".join(lines[1:]).strip()
        except requests.exceptions.Timeout:
            raise TimeoutError(f"Request to {source} timed out. Try again later.")
        except requests.exceptions.HTTPError as e:
            raise ValueError(f"{source} returned an error: {e.response.status_code} - {e.response.reason}")
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Failed to fetch sequence from {source}: {e}")
def is_ncbi_accession(accession_id: str) -> bool:
    return bool(re.match(r"^(AC_|NC_|NG_|NM_|NP_|NR_|NT_|NW_|XM_|XP_|XR_|YP_|ZP_)[0-9]+", accession_id))

def is_uniprot_accession(accession_id: str) -> bool:
    return bool(re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]{5}$", accession_id))