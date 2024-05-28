#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for CLV1 and BAM3
clv1_accession = "NM_001198022"
bam3_accession = "NM_127982"

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

def search_online_blast(sequence, gene_name):
    """Search for homologs of a given sequence using NCBI online BLAST."""
    print(f"Performing BLAST search for {gene_name}...")
    result_handle = NCBIWWW.qblast("blastn", "nt", str(sequence))
    blast_records = NCBIXML.read(result_handle)
    return blast_records

def print_blast_results(blast_records, gene_name):
    """Print BLAST search results."""
    print(f"\n{gene_name} Homologs:")
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.05:
                print(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")

def main():
    # Fetch sequences for CLV1 and BAM3 and save them to files
    clv1_sequence = fetch_sequence(clv1_accession, "clv1.NM_001198022.fna")
    bam3_sequence = fetch_sequence(bam3_accession, "bam3.NM_127982.fna")

    # Search for CLV1 homologs
    clv1_blast_results = search_online_blast(clv1_sequence, "CLV1")
    print_blast_results(clv1_blast_results, "CLV1")

    # Search for BAM3 homologs
    bam3_blast_results = search_online_blast(bam3_sequence, "BAM3")
    print_blast_results(bam3_blast_results, "BAM3")

if __name__ == "__main__":
    print("Starting CLAVATA homolog search script...")
    main()
    print("Script completed.")
