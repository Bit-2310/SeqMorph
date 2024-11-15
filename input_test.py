import unittest
from unittest.mock import patch, mock_open
from input_module import InputHandler

class TestInputHandler(unittest.TestCase):
    @patch("input_module.Entrez.efetch")
    def test_fetch_sequence_ncbi(self, mock_efetch):
        """
        Test fetching a sequence from NCBI using a mock response.
        """
        # Mock response for NCBI efetch
        mock_efetch.return_value.__enter__.return_value.read.return_value = (
            ">MockSequence\nATGCATGCATGC"
        )

        handler = InputHandler(email="test@example.com")
        accession_id = "NM_001200025"  # Example accession ID
        sequence = handler.fetch_sequence_ncbi(accession_id)

        self.assertEqual(sequence, "ATGCATGCATGC")
        mock_efetch.assert_called_once_with(
            db="nucleotide", id=accession_id, rettype="fasta", retmode="text"
        )

    @patch("requests.get")
    def test_fetch_sequence_uniprot(self, mock_get):
        """
        Test fetching a sequence from UniProt using a mock response.
        """
        # Mock response for UniProt
        mock_get.return_value.status_code = 200
        mock_get.return_value.text = ">MockProtein\nMAMAPRTEINSTRING"

        handler = InputHandler()
        uniprot_id = "P69905"  # Example UniProt ID
        sequence = handler.fetch_sequence_uniprot(uniprot_id)

        self.assertEqual(sequence, "MAMAPRTEINSTRING")
        mock_get.assert_called_once_with(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta")

if __name__ == "__main__":
    unittest.main()
