#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the BioSample accession numbers for Trifolium repens samples
biosample_accessions = ["SAMEA110451297", "SAMEA110451296", "SAMEA110450218"]

def fetch_nucleotide_accessions(biosample_accession):
    """Fetch nucleotide sequence accessions for a given BioSample accession number."""
    print(f"Fetching nucleotide accessions for BioSample {biosample_accession}...")
    handle = Entrez.esearch(db="nucleotide", term=biosample_accession)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def fetch_sequence(nucleotide_accession):
    """Fetch the DNA sequence for a given nucleotide accession number and save it to a file."""
    print(f"Fetching sequence for {nucleotide_accession}...")
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    output_file = f"{nucleotide_accession}.fna"
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
    for biosample_accession in biosample_accessions:
        nucleotide_accessions = fetch_nucleotide_accessions(biosample_accession)
        for nucleotide_accession in nucleotide_accessions:
            # Fetch sequence and save to file
            sequence = fetch_sequence(nucleotide_accession)
            
            # Perform BLAST search
            blast_results = search_online_blast(sequence, nucleotide_accession)
            
            # Print BLAST results
            print_blast_results(blast_results, nucleotide_accession)

if __name__ == "__main__":
    print("Starting BLAST search script for Trifolium repens BioSamples...")
    main()
    print("Script completed.")
