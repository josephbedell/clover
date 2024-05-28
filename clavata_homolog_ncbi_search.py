#!/usr/bin/env python3

import logging
import time
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import requests

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define your email for NCBI
Entrez.email = "joseph.bedell@gmail.com"

# Define the GenBank accession numbers for CLV1 mRNA and protein
clv1_mrna_accessions = ["XM_019238096.1", "XM_010473414.2"]
clv1_protein_accessions = ["XP_019093641.1", "XP_010471716.1"]

# Define the GenBank accession numbers for BAM3 mRNA and protein
bam3_mrna_accession = "NM_118146.5"
bam3_protein_accession = "NP_193760.1"

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

def search_online_blast(sequence, program, db, organism, accession):
    """Search for homologs of a given sequence using NCBI online BLAST with a timeout."""
    logging.info(f"Performing BLAST search for {accession} against {db} database with {program}...")
    try:
        start_time = time.time()
        result_handle = NCBIWWW.qblast(program, db, str(sequence), entrez_query=f"{organism}[Organism]")
        blast_records = NCBIXML.read(result_handle)
        end_time = time.time()
        logging.info(f"BLAST search completed in {end_time - start_time:.2f} seconds")
        return blast_records
    except requests.exceptions.Timeout:
        logging.error("BLAST search timed out.")
        return None

def print_blast_results(blast_records, accession, db):
    """Print BLAST search results."""
    if blast_records:
        logging.info(f"\n{accession} Homologs in {db} database:")
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.05:
                    logging.info(f"Sequence: {alignment.title}\nLength: {alignment.length}\nE-value: {hsp.expect}\n")
    else:
        logging.error(f"No results for {accession} in {db} database.")

def main():
    organism = "Trifolium repens"
    
    # Fetch and BLAST search for CLV1 mRNA sequences
    for accession in clv1_mrna_accessions:
        clv1_sequence = fetch_sequence(accession, f"{accession}.fna")
        clv1_blast_results_nt = search_online_blast(clv1_sequence, "blastn", "nt", organism, accession)
        print_blast_results(clv1_blast_results_nt, accession, "NT")

    # Fetch and BLAST search for CLV1 protein sequences
    for accession in clv1_protein_accessions:
        clv1_protein_sequence = fetch_protein_sequence(accession, f"{accession}.faa")
        clv1_blast_results_nr = search_online_blast(clv1_protein_sequence, "blastp", "nr", organism, accession)
        print_blast_results(clv1_blast_results_nr, accession, "NR")

    # Fetch and BLAST search for BAM3 mRNA sequence
    bam3_sequence = fetch_sequence(bam3_mrna_accession, f"{bam3_mrna_accession}.fna")
    bam3_blast_results_nt = search_online_blast(bam3_sequence, "blastn", "nt", organism, bam3_mrna_accession)
    print_blast_results(bam3_blast_results_nt, bam3_mrna_accession, "NT")

    # Fetch and BLAST search for BAM3 protein sequence
    bam3_protein_sequence = fetch_protein_sequence(bam3_protein_accession, f"{bam3_protein_accession}.faa")
    bam3_blast_results_nr = search_online_blast(bam3_protein_sequence, "blastp", "nr", organism, bam3_protein_accession)
    print_blast_results(bam3_blast_results_nr, bam3_protein_accession, "NR")

if __name__ == "__main__":
    logging.info("Starting BLAST search script for CLV1 and BAM3 mRNA and protein against NT and NR databases...")
    main()
    logging.info("Script completed.")
