#!/usr/bin/env python3

import unittest
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from clavata_homolog_ncbi_search import fetch_sequence, fetch_protein_sequence, search_online_blast, print_blast_results

class TestBlastSearch(unittest.TestCase):

    def test_fetch_sequence(self):
        seq = fetch_sequence("XM_019238096.1", "test.fna")
        self.assertTrue(seq.startswith("ATG"))

    def test_fetch_protein_sequence(self):
        seq = fetch_protein_sequence("XP_019093641.1", "test.faa")
        self.assertTrue(seq.startswith("M"))

    def test_search_online_blast(self):
        result = search_online_blast("ATGCGTACG", "blastn", "refseq", "Trifolium repens", "XM_019238096.1")
        self.assertIsNotNone(result)

    def test_print_blast_results(self):
        mock_alignment = MagicMock(title="Test Sequence", length=100, hsps=[MagicMock(expect=0.01)])
        mock_blast_record = MagicMock(alignments=[mock_alignment])
        with patch('clavata_homolog_ncbi_search.logging.info') as mock_logging:
            print_blast_results(mock_blast_record, "XM_019238096.1", "NT")
            mock_logging.assert_any_call("\nXM_019238096.1 Homologs in NT database:")

if __name__ == '__main__':
    unittest.main()
