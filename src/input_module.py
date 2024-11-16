# input_module.py: Handles file uploads and sequence retrieval for SeqMorph

import os
from Bio import SeqIO
from tkinter import Tk, filedialog
import requests

class InputHandler:
    def __init__(self):
        pass

    def select_file(self):
        root = Tk()  # Create the root window
        root.withdraw()  # Hide the root window
        file_path = filedialog.askopenfilename(
            title="Select a Sequence File",
            filetypes=(
                ("FASTA Files", "*.fasta"),
                ("Text Files", "*.txt"),
                ("All Files", "*.*")
            )
        )
        return file_path

    def select_multiple_files(self):
        root = Tk()  # Create the root window
        root.withdraw()  # Hide the root window
        file_paths = filedialog.askopenfilenames(
            title="Select Sequence Files",
            filetypes=(
                ("FASTA Files", "*.fasta"),
                ("Text Files", "*.txt"),
                ("All Files", "*.*")
            )
        )
        return file_paths

    def load_sequence_file(self, filepath):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        file_extension = os.path.splitext(filepath)[1].lower()
        if file_extension not in [".fasta", ".txt"]:
            raise ValueError(f"Unsupported file format: {file_extension}")

        try:
            if file_extension == ".fasta":
                with open(filepath, "r") as file:
                    sequences = [str(record.seq) for record in SeqIO.parse(file, "fasta")]
                    if not sequences:
                        raise ValueError("No sequences found in the FASTA file.")
                    return sequences[0]  # Return the first sequence
            elif file_extension == ".txt":
                with open(filepath, "r") as file:
                    sequence = file.read().strip()
                    if not sequence:
                        raise ValueError("No content found in the TXT file.")
                    return sequence
        except Exception as e:
            raise ValueError(f"Error reading file: {str(e)}")

    @staticmethod
    def sequence_type(sequence):
        """
        Detect the type of a sequence (DNA, RNA, Protein).

        Args:
            sequence (str): The input sequence to analyze.

        Returns:
            str: The type of sequence ('DNA', 'RNA', 'Protein').

        Raises:
            ValueError: If the sequence type cannot be determined.
        """
        if set(sequence.upper()).issubset("ATGC"):
            return "DNA"
        elif set(sequence.upper()).issubset("AUGC"):
            return "RNA"
        elif set(sequence.upper()).issubset("ACDEFGHIKLMNPQRSTVWY"):
            return "Protein"
        else:
            raise ValueError("Cannot detect sequence type.")
        
    def fetch_sequence_ncbi(self, accession_id):
        url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id={accession_id}&report=fasta"
        try:
            response = requests.get(url)
            response.raise_for_status()
            sequence_data = "".join(response.text.splitlines()[1:])  # Remove header line
            if not sequence_data:
                raise ValueError("No sequence data found for the given NCBI accession ID.")
            return sequence_data
        except Exception as e:
            raise ValueError(f"Error fetching sequence from NCBI: {str(e)}")

    def fetch_sequence_uniprot(self, uniprot_id):
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        try:
            response = requests.get(url)
            response.raise_for_status()
            sequence_data = "".join(response.text.splitlines()[1:])  # Remove header line
            if not sequence_data:
                raise ValueError("No sequence data found for the given UniProt ID.")
            return sequence_data
        except Exception as e:
            raise ValueError(f"Error fetching sequence from UniProt: {str(e)}")