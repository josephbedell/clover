#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for Trifolium repens samples
accessions = ["SAMEA110451297", "SAMEA110451296", "SAMEA110450218"]

def fetch_sequence(accession):
    """Fetch the DNA sequence for a given GenBank accession number."""
    print(f"Fetching sequence for {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    output_file = f"{accession}.fna"
    print(f"Saving sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return record.seq

def search_online_blast(sequence, accession):
    """Search for homologs of a given sequence using NCBI online BLAST."""
    print(f"Performing BLAST search for {accession}...")
    result_handle = NCBIWWW.qblast("blastn", "nt", str(sequence))
    blast_records = NCBIXML.read(result_handle)
    return blast_records

def print_blast_results(blast_records, accession):
    """Print BLAST search results."""
    print(f"\n{accession} Homologs:")
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.05:
                print(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")

def main():
    for accession in accessions:
        # Fetch sequence and save to file
        sequence = fetch_sequence(accession)
        
        # Perform BLAST search
        blast_results = search_online_blast(sequence, accession)
        
        # Print BLAST results
        print_blast_results(blast_results, accession)

if __name__ == "__main__":
    print("Starting BLAST search script for Trifolium repens accessions...")
    main()
    print("Script completed.")
