#!/usr/bin/env python3

import logging
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for CLV1 and BAM3
clv1_accession = "NM_001198022"
bam3_accession = "NM_127982"

def fetch_sequence(accession, output_file):
    """Fetch the DNA sequence for a given GenBank accession number and save it to a file."""
    logging.info(f"Fetching sequence for {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    logging.info(f"Saving sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return record.seq

def fetch_protein_sequence(accession, output_file):
    """Fetch the protein sequence for a given GenBank accession number and save it to a file."""
    logging.info(f"Fetching protein sequence for {accession}...")
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    logging.info(f"Saving protein sequence to {output_file}...")
    with open(output_file, "w") as file:
        SeqIO.write(record, file, "fasta")
    return record.seq

def search_online_blast(sequence, db, organism, accession):
    """Search for homologs of a given sequence using NCBI online BLAST."""
    logging.info(f"Performing BLAST search for {accession} against {db} database...")
    result_handle = NCBIWWW.qblast("blastn" if db == "nt" else "blastp", db, str(sequence), entrez_query=f"{organism}[Organism]")
    blast_records = NCBIXML.read(result_handle)
    return blast_records

def print_blast_results(blast_records, accession, db):
    """Print BLAST search results."""
    logging.info(f"\n{accession} Homologs in {db} database:")
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.05:
                logging.info(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")

def main():
    # Fetch DNA sequences for CLV1 and BAM3
    clv1_sequence = fetch_sequence(clv1_accession, "clv1.NM_001198022.fna")
    bam3_sequence = fetch_sequence(bam3_accession, "bam3.NM_127982.fna")

    # Fetch protein sequences for CLV1 and BAM3
    clv1_protein_sequence = fetch_protein_sequence(clv1_accession, "clv1.NM_001198022.faa")
    bam3_protein_sequence = fetch_protein_sequence(bam3_accession, "bam3.NM_127982.faa")

    # Perform BLAST searches against NT and NR databases
    organism = "Trifolium repens"
    
    # BLAST searches for DNA sequences against NT
    clv1_blast_results_nt = search_online_blast(clv1_sequence, "nt", organism, clv1_accession)
    print_blast_results(clv1_blast_results_nt, clv1_accession, "NT")

    bam3_blast_results_nt = search_online_blast(bam3_sequence, "nt", organism, bam3_accession)
    print_blast_results(bam3_blast_results_nt, bam3_accession, "NT")
    
    # BLAST searches for protein sequences against NR
    clv1_blast_results_nr = search_online_blast(clv1_protein_sequence, "nr", organism, clv1_accession)
    print_blast_results(clv1_blast_results_nr, clv1_accession, "NR")

    bam3_blast_results_nr = search_online_blast(bam3_protein_sequence, "nr", organism, bam3_accession)
    print_blast_results(bam3_blast_results_nr, bam3_accession, "NR")

if __name__ == "__main__":
    logging.info("Starting BLAST search script for CLV1 and BAM3 against NT and NR databases...")
    main()
    logging.info("Script completed.")
