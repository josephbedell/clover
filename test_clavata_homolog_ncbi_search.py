#!/usr/bin/env python3

import unittest
from unittest.mock import patch, MagicMock
from io import StringIO
from Bio.Blast import NCBIXML
from clavata_homolog_ncbi_search import fetch_sequence, fetch_protein_sequence, search_online_blast, print_blast_results

class TestBlastSearch(unittest.TestCase):

    @patch('clavata_homolog_ncbi_search.Entrez.efetch')
    def test_fetch_sequence(self, mock_efetch):
        mock_efetch.return_value = StringIO(">XM_019238096.1\nATGCGTACG")
        seq = fetch_sequence("XM_019238096.1", "test.fna")
        self.assertEqual(str(seq), "ATGCGTACG")

    @patch('clavata_homolog_ncbi_search.Entrez.efetch')
    def test_fetch_protein_sequence(self, mock_efetch):
        mock_efetch.return_value = StringIO(">XP_019093641.1\nMKSTFF")
        seq = fetch_protein_sequence("XP_019093641.1", "test.faa")
        self.assertEqual(str(seq), "MKSTFF")

    @patch('clavata_homolog_ncbi_search.NCBIWWW.qblast')
    @patch('clavata_homolog_ncbi_search.NCBIXML.read')
    def test_search_online_blast(self, mock_read, mock_qblast):
        mock_qblast.return_value = StringIO("<BlastOutput></BlastOutput>")
        mock_read.return_value = MagicMock(alignments=[MagicMock(hsps=[MagicMock(expect=0.01)])])
        result = search_online_blast("ATGCGTACG", "blastn", "nt", "Trifolium repens", "XM_019238096.1")
        self.assertIsNotNone(result)

    @patch('clavata_homolog_ncbi_search.logging.info')
    def test_print_blast_results(self, mock_logging):
        mock_alignment = MagicMock(title="Test Sequence", length=100, hsps=[MagicMock(expect=0.01)])
        mock_blast_record = MagicMock(alignments=[mock_alignment])
        print_blast_results(mock_blast_record, "XM_019238096.1", "NT")
        mock_logging.assert_any_call("\nXM_019238096.1 Homologs in NT database:")

if __name__ == '__main__':
    unittest.main()
