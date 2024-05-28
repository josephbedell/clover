#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for CLV1 and BAM3
clv1_accession = "NM_001198022"
bam3_accession = "NM_127982"

# Define the SRA accession numbers for Trifolium repens samples
sra_accessions = ["ERS12548248"]

def fetch_sequence(accession, output_file):
    """Fetch the DNA sequence for a given GenBank accession number and save it to a file."""
    print(f"Fetching sequence for {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print(f"Saving sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return record.seq

def search_online_blast(sequence, db, accession):
    """Search for homologs of a given sequence using NCBI online BLAST."""
    print(f"Performing BLAST search for {accession} against {db}...")
    result_handle = NCBIWWW.qblast("blastn", db, str(sequence))
    blast_record
